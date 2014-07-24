#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

# Script to dump Ensembl triples. It was lashed together rapidly, and could stand to use a library for the writing of triples, with accommodation for the scale of the data involved.
# Requires installation of Ensembl Core and Compara APIs, as well as dependencies such as BioPerl


use strict;

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $species = '';
my $version = '';
my $path;
my $virtgraph;
my $man;
my $limit;
my $help;

GetOptions (
    'species=s' => \$species,
    'version=s' => \$version,
    'virtgraph' => \$virtgraph,
    'out=s' => \$path,
    'limit=i' => \$limit,
    'help|?' => \$help, 
    man => \$man
    ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

if (!$species || !$version) {
    pod2usage(1);
}
print STDOUT "Generating mappings using Ensembl core DB $species version $version\n";

my $outfile = ${species}.'_'.${version}.'_xrefs.ttl';
if ($path) {
    $path = $path."/";
    $outfile = $path.$outfile;
}

# create the output file
open OUT, ">$outfile" || die "Can't open out file $outfile\n";

my %prefix = (
  ensembl => 'http://rdf.ebi.ac.uk/resource/ensembl/',
  term => 'http://rdf.ebi.ac.uk/terms/ensembl/',
  rdfs => 'http://www.w3.org/2000/01/rdf-schema#',
  sio => 'http://semanticscience.org/resource/',
  dc => 'http://purl.org/dc/terms/',
  rdf => 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
  faldo => 'http://biohackathon.org/resource/faldo#',
  obo => 'http://purl.obolibrary.org/obo/',
  skos => 'http://www.w3.org/2004/02/skos/core#',
  taxon => 'http://identifiers.org/taxonomy/',
);

foreach (keys %prefix) {
    triple('@prefix',$_.':',u($prefix{$_}) );
}

# Choice of database host is a factor in how fast the script runs. Try to find your nearest mirror, and check the database version before running.
Bio::EnsEMBL::Registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
  -db_version => 75,
  -no_cache => 1,
);

my $ga = Bio::EnsEMBL::Registry->get_adaptor($species,'Core','Gene');

my $meta = Bio::EnsEMBL::Registry->get_adaptor($species,'Core','MetaContainer');

# create a map of taxon id to species name, we will create some triples for these at the end

# get the taxon id for this species 
my $taxon_id = $meta->get_taxonomy_id;
my $scientific_name = $meta->get_scientific_name;

my %dbname2type;
my %dbname2base;
my %dbname2short;
my %dbid2regex;
my %ignore;

# check if we want to create a virtuoso graph file for this output file
if ($virtgraph) {
    my $versionGraphUri = "http://rdf.ebi.ac.uk/dataset/ensembl/".$version;
    my $graphUri = $versionGraphUri."/".$taxon_id;
    my $graphFile = $outfile.'.graph';
    open GRAPH, ">$graphFile" || die "Can't create the virtuoso graph file\n";
    print GRAPH $graphUri;
    close GRAPH;
    # make the species graph a subgraph of the version graph
    triple (u($graphUri), '<http://www.w3.org/2004/03/trix/rdfg-1/subGraphOf>', u($versionGraphUri)); 
}


# read config
open DBNAME , '<xref_config.txt' || die "can't open mapping config file\n";
while (<DBNAME>) {
    chomp;
    my $line = $_;
    my @split = split /\t/, $line;
    my $dbname = $split[0];
    my $baseUri = $split[1];
    my $typeUri = $split[2];
    my $shortname = $split[3];
    my $regex = $split[4];
    my $na = $split[5];
    
    if ($na) {
	$ignore{$dbname} =1;
    }
    if ($shortname) {
	$dbname2short{$dbname} = $shortname;
    }
    if ($regex) {
	$dbid2regex{$dbname} = $regex;
    }
    if ($typeUri) {
	$dbname2type{$dbname} = $typeUri;
    }
    if ($baseUri) {
	$dbname2base{$dbname} = $baseUri;
    }
}
close DBNAME;

my %relations;
my %ensemblSources;
my $genes = $ga->fetch_all();

#print "GENE ", $gene->stable_id(), "\n";
my $count = 0;
while (my $gene = shift @$genes) {
    $count++;
    print_DBEntries( $gene->get_all_DBEntries() , $gene->stable_id());
    
    foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
	print_DBEntries( $transcript->get_all_DBEntries() , $transcript->stable_id());
	
	if ( defined $transcript->translation() ) {
	    my $translation = $transcript->translation();
	    print_DBEntries( $translation->get_all_DBEntries(), $translation->stable_id() );
	}
    }
  last if ($limit && $count == $limit);
}

# print relation assertions as sub property of skos:related
open RELOUT, ">${path}${species}_xref_relations.txt" || die "can't open ${species} xref relations\n";
open SRCOUT, ">${path}${species}_xref_sources.txt" || die "can't open ${species} xref relations\n";

while ( my ($key, $value) = each(%relations) ) {
    print RELOUT "$key\n";
}

while ( my ($key, $relhash) = each(%ensemblSources) ) {
    
    while ( my ($rel, $value) = each(%$relhash) ) {	
	print SRCOUT "${key}\t${rel}\t${value}\n";
    }
}

close RELOUT;
close SRCOUT;

sub print_DBEntries
{
    my $db_entries = shift;
    my $ensId = shift;
    foreach my $dbe ( @{$db_entries} ) {

	#printf "\tXREF %s (db:%s) (id:%s) (desc:%s) (enstype:%s) (exttype:%s) (evidence:%s) (linkage:%s)\n",
	my $label = $dbe->display_id();
	my $name = $dbe->dbname();
	next if ($ignore{$name});
	my $id = $dbe->primary_id();
	my $desc = $dbe->description();

	my $ensemblUri = $prefix{ensembl}.$ensId; 
	my $relation = $prefix{term}.$dbe->info_type(); 
	$relations{$relation}=1;

	if (!$ensemblSources{$name}) {
	    $ensemblSources{$name}{$relation} = $id;
	    triple(u($relation), 'rdfs:subPropertyOf', 'skos:related');
	}
	elsif (!$ensemblSources{$name}{$relation}) {
	    $ensemblSources{$name}{$relation} = $id;
	    triple(u($relation), 'rdfs:subPropertyOf', 'skos:related');
	}

	my $xrefUri;
	my $xrefTypeUri;
	
	if ($dbid2regex{$name}) {
	    eval "\$id =~ $dbid2regex{$name}";
	}

	if ($dbname2base{$name}) {
	    $xrefUri= $dbname2base{$name}.$id;
	}
	elsif ($dbname2short{$name}) {
	    $xrefUri = "http://identifiers.org/".$dbname2short{$name}."/".$id;
	}
	else {
	    $xrefUri = "http://identifiers.org/".lc($name)."/".$id;
	}
	
	if ($dbname2type{$name}) {
	    $xrefTypeUri = $dbname2type{$name};
	    triple(u($xrefTypeUri), 'rdfs:subClassOf', 'ensembl:EnsemblDBEntry');
	    triple(u($xrefTypeUri), 'rdfs:label', '"'.$name.'"');
	}
	else {
	    $xrefTypeUri = $prefix{term}.$name;	    
	    triple(u($xrefTypeUri), 'rdfs:subClassOf', 'ensembl:EnsemblDBEntry');
	    triple(u($xrefTypeUri), 'rdfs:label', '"'.$name.'"');
	}
	
	triple(u($ensemblUri), u($relation), u($xrefUri));
	triple(u($xrefUri), 'rdfs:label', '"'.$label.'"');
	triple(u($xrefUri), 'dc:identifier', '"'.$dbe->primary_id().'"');
	triple(u($xrefUri), 'dc:description', '"'.$desc.'"');

	# type the xref
	triple(u($xrefUri), 'rdfs:subClassOf', u($xrefTypeUri));
	
    }
}


sub u {
    my $stuff= shift;
    return '<'.$stuff.'>';
}

sub triple {
    my ($subject,$predicate,$object) = @_;
    
    printf OUT "%s %s %s .\n",$subject,$predicate,$object;
}

