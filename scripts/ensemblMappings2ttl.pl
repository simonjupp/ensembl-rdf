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
use Digest::MD5 qw(md5_hex);

my $species = '';
my $version = '';
my $path;
my $virtgraph;
my $man;
my $limit;
my $help;
#"ensembldb.ensembl.org";
my $host = "mysql-ensembl-mirror.ebi.ac.uk"; 
my $port = "4240";
my $genome = "";

GetOptions (
    'species=s' => \$species,
    'version=s' => \$version,
    'virtgraph' => \$virtgraph,
    'out=s' => \$path,
    'limit=i' => \$limit,
    'host=s' => \$host,
    'port=i' => \$port,
    'genome=s' => \$genome,
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
  dataset => 'http://rdf.ebi.ac.uk/dataset/ensembl/',
  term => 'http://rdf.ebi.ac.uk/terms/ensembl/',
  rdfs => 'http://www.w3.org/2000/01/rdf-schema#',
  void => 'http://rdfs.org/ns/void#',
  sio => 'http://semanticscience.org/resource/',
  dc => 'http://purl.org/dc/terms/',
  rdf => 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
  prov => 'http://www.w3.org/ns/prov#',
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
  -host => $host,
#  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
  -port => $port,
  -db_version => $version,
#  -no_cache => 1,
);

my $ga = Bio::EnsEMBL::Registry->get_adaptor($species,'Core','Gene');

my $meta = Bio::EnsEMBL::Registry->get_adaptor($species,'Core','MetaContainer');

my $pfa = Bio::EnsEMBL::Registry->get_adaptor($species,"funcgen","ProbeFeature");
my $aa = Bio::EnsEMBL::Registry->get_adaptor($species,"funcgen","Array");
my $pba = Bio::EnsEMBL::Registry->get_adaptor($species,"funcgen","ProbeSet");

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
my $versionGraphUri = "http://rdf.ebi.ac.uk/dataset/ensembl/".$version;
my $graphUri = $versionGraphUri."/".$taxon_id;
if ($virtgraph) {
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
    next if $. < 2; # Skip first line

    chomp;
    my $line = $_;
    my @split = split /\t/, $line;
    my $dbname = $split[0];
    my $shortname = $split[2];
    my $baseUri = $split[3];
    my $typeUri = $split[4];
    my $na = $split[5];
    my $regex = $split[6];
    
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
#my $gene = $ga->fetch_by_stable_id ('ENSG00000105393');
#print "GENE ", $gene->stable_id(), "\n";
my $count = 0;
while (my $gene = shift @$genes) {
    $count++;
    print_DBEntries( $gene->get_all_DBEntries() , $gene, undef, $gene->biotype());
    
    foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
	print_DBEntries( $transcript->get_all_DBEntries() , $gene, 'INFERRED_FROM_TRANSCRIPT', $gene->biotype());
	print_DBEntries( $transcript->get_all_DBEntries() , $transcript,  undef, "transcript");
#	print_Probe_Features ($transcript);
	
	if ( defined $transcript->translation() ) {
	    my $translation = $transcript->translation();
	    print_DBEntries( $translation->get_all_DBEntries() , $gene, 'INFERRED_FROM_TRANSLATION', $gene->biotype());
	    print_DBEntries( $translation->get_all_DBEntries(), $translation, undef,  "protein");
	}
    }
  last if ($limit && $count == $limit);
}

# print relation assertions as sub property of skos:related
#open RELOUT, ">${path}${species}_xref_relations.txt" || die "can't open ${species} xref relations\n";
open SRCOUT, ">${path}${species}_xref_sources.txt" || die "can't open ${species} xref relations\n";
open LINKSETOUT, ">${path}${species}_linksets_void.ttl" || die "can't open ${species} xref relations\n";

#while ( my ($key, $value) = each(%relations) ) {
#    print RELOUT "$key\n";
#}

my %linksets;
while ( my ($type, $nameHash) = each(%ensemblSources) ) {
    
    while ( my ($objectName, $relhash) = each(%$nameHash) ) {
	next if ($objectName ne "Uniprot/SPTREMBL");
	while ( my ($rel, $subjectNameHash) = each(%$relhash) ) {	
	    while ( my ($subjectName, $count) = each(%$subjectNameHash) ) {	
		print SRCOUT "${objectName}\t${type}\t${rel}\t${subjectName}\t${count}\n";
		
		my $linksetid = "linkset-" . md5_hex( ($objectName, $rel, $subjectName));
		my $linksetLabel = $subjectName . " " . $rel . " " . $objectName . " linkset";
		$linksets{$linksetid} = $linksetLabel;
		my $linksetUri = $prefix{dataset}.$version."/".$linksetid;
		my $subjectPartitionUri = $prefix{dataset}.$version."/".$subjectName."-dataset-partition";
		my $objectPartitionUri .=  $prefix{dataset}.$version."/".$objectName."-dataset-partition";
		
		# define the partitions
		printf LINKSETOUT "%s %s %s .\n",u( $prefix{dataset}.$version), u($prefix{void}."classPartition"), u($subjectPartitionUri);
		printf LINKSETOUT "%s %s %s .\n",u( $prefix{dataset}.$version), u($prefix{void}."classPartition"), u($objectPartitionUri);
		printf LINKSETOUT "%s %s %s .\n",u($subjectPartitionUri), u($prefix{void}."class"), u($prefix{term}.$subjectName);
		printf LINKSETOUT "%s %s %s .\n",u($objectPartitionUri), u($prefix{void}."class"), u($type);
		
		# define the linksets 
		printf LINKSETOUT "%s %s %s .\n",u($linksetUri), u($prefix{void}."linkPredicate"), u($prefix{term}.$rel);
		printf LINKSETOUT "%s %s %s .\n",u($linksetUri), u($prefix{void}."subjectTarget"), u($subjectPartitionUri);
		printf LINKSETOUT "%s %s %s .\n",u($linksetUri), u($prefix{void}."objectTarget"), u($objectPartitionUri);
		printf LINKSETOUT "%s %s %s .\n",u($linksetUri), u($prefix{void}."triples"), $count;

		# if inferred link, link to provenance
		if ($rel =~/^INFERRED_FROM/) {
		    foreach my $inferredFromRel  (keys %{$relations{$objectName}} ) {
			if ($inferredFromRel != $rel) {
			    my $inferredLinksetid = "linkset-" . md5_hex( ($objectName, $inferredFromRel, "transcript"));
			    if ($rel =~/TRANSLATION/) {
				$inferredLinksetid = "linkset-" . md5_hex( ($objectName, $inferredFromRel, "protein"));
			    }
			    my $inferredLinksetUri = $prefix{dataset}.$version."/".$inferredLinksetid;
#			    my $inferredLinksetLabel = $subjectName . " " . $rel . " " . $objectName . " linkset";
#			    $linksets{$inferredLinksetid} = $inferredLinksetLabel;
			    
			    # derived from subject partition, transcript partition and protein partition
			    printf LINKSETOUT "%s %s %s .\n",u($linksetUri), u($prefix{prov}."wasDerivedFrom"), u($inferredLinksetUri);
			} 
		    }
	#		print "\n";
		    my $geneTranscriptPartitionUri = $prefix{dataset}.$version."/linkset-transcript-SO_transcribed_from-".$subjectName;
		    printf LINKSETOUT "%s %s %s .\n",u($linksetUri), u($prefix{prov}."wasDerivedFrom"), u($geneTranscriptPartitionUri);
		    if ($rel =~/TRANSLATION/) {
			my $transcriptProteinPartitionUri = $prefix{dataset}.$version."/linkset-transcript-SO_translates_to-protein";
			printf LINKSETOUT "%s %s %s .\n",u($linksetUri), u($prefix{prov}."wasDerivedFrom"), u($transcriptProteinPartitionUri);
		    }
		}
	    }
#	print SRCOUT "${key}\t${rel}\t${value}\n";
	}
    }
}

while (my ($key, $value) = each (%linksets)) {
    printf LINKSETOUT "%s %s %s .\n", u($prefix{dataset}.$version."/".$key), "a", u($prefix{void}."Linkset");
    printf LINKSETOUT "%s %s \"%s\" .\n", u($prefix{dataset}.$version."/".$key), u($prefix{rdfs}."label"),  $value;
}

printf LINKSETOUT "%s %s %s .\n", u($prefix{dataset}.$version."/".$taxon_id), "a", u($prefix{void}."Dataset");
printf LINKSETOUT "%s %s \"%s Ensembl %s RDF\" .\n", u($prefix{dataset}.$version."/".$taxon_id), u($prefix{dc}."title"), $scientific_name, $version;
printf LINKSETOUT "%s %s <ftp://ftp.ebi.ac.uk/pub/databases/RDF/ensembl/%s/%s> .\n", u($prefix{dataset}.$version."/".$taxon_id), u($prefix{void}."datadump"), $version, ${species}.'_'.${version}.'.ttl';


sub print_Probe_Features {
    my $transcript = shift;
    
    my @probesets = @{$pba->fetch_all_by_external_name($transcript->stable_id)};

    foreach my $probeset (@probesets){
	
	my $arrays_string = "";
	foreach my $array (@{$probeset->get_all_Arrays}) {
	    
	    print "probeset name: " . $probeset->name . "\n";
	    print "array desc: " . $array->description . "\n";
	    print "array name: " . $array->name . "\n";
	    print "array vendor: " . $array->vendor . "\n";
	    print "array type: " . $array->type . "\n";
	    $arrays_string .= $array->name;
	    
	}
    }
}

#close RELOUT;
close SRCOUT;
close LINKSETOUT;

sub print_DBEntries
{
    my $db_entries = shift;
    my $feature = shift;
    my $rel = shift;
    my $biotype = shift;

    foreach my $dbe ( @{$db_entries} ) {

	#printf "\tXREF %s (db:%s) (id:%s) (desc:%s) (enstype:%s) (exttype:%s) (evidence:%s) (linkage:%s)\n",
	my $label = $dbe->display_id();
	my $name = $dbe->dbname();
	next if ($ignore{$name});
	my $id = $dbe->primary_id();
	my $desc = $dbe->description();

	my $ensemblUri = $prefix{ensembl}.$feature->stable_id(); 
	if (!$rel) {
	    $rel = $dbe->info_type();
	}
	my $relation = $prefix{term}.$rel; 
	$relations{$name}{$rel}=1;


	my $xrefUri;
	my $xrefTypeUri;
	

	if ($dbname2short{$name}) {
	    $xrefUri = "http://identifiers.org/".$dbname2short{$name}."/".$id;
	}
	

	if ($dbname2type{$name}) {
	    $xrefTypeUri = $dbname2type{$name};
	    triple(u($xrefTypeUri), 'rdfs:subClassOf', 'term:EnsemblExternalDBEntry');
	    triple(u($xrefTypeUri), 'rdfs:label', '"'.$name.'"');
	}
	else {
	    $xrefTypeUri = $prefix{term}.$name;	    
	    triple(u($xrefTypeUri), 'rdfs:subClassOf', 'term:EnsemblExternalDBEntry');
	    triple(u($xrefTypeUri), 'rdfs:label', '"'.$name.'"');
	}

	if (!$ensemblSources{$xrefTypeUri}) {
	    $ensemblSources{$xrefTypeUri}{$name}{$rel}{$biotype} = 1;
	    triple(u($relation), 'rdfs:subPropertyOf', 'skos:related');
	}
	elsif (!$ensemblSources{$xrefTypeUri}{$name}{$rel}{$biotype}) {
	    $ensemblSources{$xrefTypeUri}{$name}{$rel}{$biotype} = 1;
	    triple(u($relation), 'rdfs:subPropertyOf', 'skos:related');
	}
	$ensemblSources{$xrefTypeUri}{$name}{$rel}{$biotype}++;
	
	triple(u($ensemblUri), u($relation), u($xrefUri));
	triple(u($xrefUri), 'rdfs:label', '"'.$label.'"');
	triple(u($xrefUri), 'dc:identifier', '"'.$dbe->primary_id().'"');
	if ($desc) {
	    triple(u($xrefUri), 'dc:description', '"'.escape($desc).'"');
	}
	# type the xref
	triple(u($xrefUri), 'a', u($xrefTypeUri));

	# add a canoncial LOD URI if provided
	if ($dbname2base{$name}) {
	    $xrefUri= $dbname2base{$name}.$id;
	    if ($dbid2regex{$name}) {
		eval "\$id =~ $dbid2regex{$name}";
	    }
	    $xrefUri= $dbname2base{$name}.$id;
	    triple(u($ensemblUri), u($relation), u($xrefUri));
	    triple(u($xrefUri), 'rdfs:label', '"'.$label.'"');
	    triple(u($xrefUri), 'dc:identifier', '"'.$dbe->primary_id().'"');
	    if ($desc) {
		triple(u($xrefUri), 'dc:description', '"'.escape($desc).'"');
	    }
	    # type the xref
	    triple(u($xrefUri), 'a', u($xrefTypeUri));
	}
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

sub escape {
    my $string = shift;
    $string =~s/(["])/\\$1/g;
    return $string;
}
