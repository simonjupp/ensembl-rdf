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

=head1 NAME

    ensembl2ttl - Convert Ensembl core data to RDF turtle

=head1 AUTHOR

Kieron Taylor, Simon Jupp, European Bioinformatics Institute

=head1 SYNOPSIS
    ensembl2ttl [options]
     Options:
       -species         The common name for species to conver e.g. Human
       -version         The Ensembl version e.g. 75 
       -help            brief help message
       -man             full documentation

=head1 OPTIONS

=over 8

=item B<-species>

       The common name for species to conver e.g. Human

=item B<-version>

       The Ensembl version e.g. 75

=item B<-help>

    Print a brief help message and exits.

=item B<-man>

    Prints the manual page and exits.

=back

=head1 DESCRIPTION

    Script to dump Ensembl triples. It was lashed together rapidly, and could stand to use a library for the writing of triples, with accommodation for the scale of the data involved.
    Requires installation of Ensembl Core and Compara APIs, as well as dependencies such as BioPerl

=cut


use strict;

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $species = '';
my $version = '';
my $path;
my $split;
my $virtgraph;
my $man = 0;
my $limit;
my $help = 0;

GetOptions (
    'species=s' => \$species,
    'version=s' => \$version,
    'virtgraph' => \$virtgraph,
    'out=s' => \$path,
    'split=i' => \$split,
    'limit=i' => \$limit,
    'help|?' => \$help, 
    man => \$man
    ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

if (!$species || !$version) {
    pod2usage(1);
}
print STDOUT "Using Ensembl core DB $species version $version\n";

my $outfile = ${species}.'_'.${version}.'.ttl';
if ($path) {
    $outfile = $path."/".$outfile;
}

# create the output file
open OUT, ">$outfile" || die "Can't open out file $outfile\n";

# common prefixes used
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

# create a map of prefixes for easy lookup
foreach (keys %prefix) {
  triple('@prefix',$_.':',u($prefix{$_}) );
}

# Choice of database host is a factor in how fast the script runs. Try to find your nearest mirror, and check the database version before running.
Bio::EnsEMBL::Registry->load_registry_from_db(
  -host => 'mysql-ensembl-mirror.ebi.ac.uk',
#  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
  -port => 4240,
  -db_version => $version,
  -no_cache => 1,
);


my $ga = Bio::EnsEMBL::Registry->get_adaptor($species,'Core','Gene');

my $ontoa =
    Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

# get all genes
my $genes = $ga->fetch_all();

# create a lookup table for term to ontology id
my %term2ontologyId;

# simple count for testing
my $count = 0;

my $meta = Bio::EnsEMBL::Registry->get_adaptor($species,'Core','MetaContainer');

# create a map of taxon id to species name, we will create some triples for these at the end

# get the taxon id for this species 
my $taxon_id = $meta->get_taxonomy_id;
my $scientific_name = $meta->get_scientific_name;

# print out global triples about the organism  
triple('taxon:'.$taxon_id, 'rdfs:subClassOf', 'obo:OBI_0100026');
triple('taxon:'.$taxon_id, 'rdfs:label', '"'.$scientific_name.'"');

# get the current ensembl release number 
my $schemaVersion = $meta->get_schema_version;

# check if we want to create a virtuoso graph file for this output file
if ($virtgraph) {
    my $versionGraphUri = "http://rdf.ebi.ac.uk/dataset/ensembl/".$version;
    my $graphUri = $versionGraphUri."/".$taxon_id;
    my $graphFile = $outfile.'.graph';
    open GRAPH, ">$graphFile" || die "Can't create the virtuoso graph file\n";
    print GRAPH $graphUri;
    close GRAPH;
    # make the species graph a subgraph of the version graph
    triple (u($graphUri), '<http://rdfs.org/ns/void#subset>', u($versionGraphUri)); 
}

# start to process all genes
#my $gene = $ga->fetch_by_stable_id('ENSG00000189167');
# dump all features#
while (my $gene = shift @$genes) {
  $count++;

  # get all the trancripts for this gene
  my @trans = @{$gene->get_all_Transcripts};

  # get the bio type and try and match to a SO ontology term, default to creating a new term
  # in the ensmebl namespace if no SO term is found todo: fix this to map to something in SO
  # most are protein coding, but we also have nonsence_mediated_decay and miscmRNA that don't map 
  my $ontoTypeId = getSOOntologyId($gene->biotype());

  # create the gene as a sublclass of the type coming back from SO, usually protein coding. 
  triple('ensembl:'.$gene->stable_id, 'rdfs:subClassOf', $ontoTypeId );

  dump_feature($gene);

  # add some useful meta data
  triple('ensembl:'.$gene->stable_id, 'rdfs:label', ($gene->external_name)? '"'.$gene->external_name.'"':'"'.$gene->display_id.'"' );
  triple('ensembl:'.$gene->stable_id, 'dc:description', '"'.$gene->description.'"');
  dump_synonyms($gene);
  
  # relate gene to its taxon
  taxonTriple('ensembl:'.$gene->stable_id);

  # loop through the transcripts
  foreach my $transcript (@trans) {
      # transcripts are transcribed from a gene 
      triple( 'ensembl:'.$transcript->stable_id, 'obo:SO_transcribed_from', 'ensembl:'.$gene->stable_id);
      taxonTriple('ensembl:'.$transcript->stable_id);
      
      # type the transcript using SO again
      triple( 'ensembl:'.$transcript->stable_id, 'rdfs:subClassOf', 'obo:SO_0000673');
      triple( 'ensembl:'.$transcript->stable_id, 'dc:identifier', '"'.$transcript->stable_id.'"');
      triple( 'ensembl:'.$transcript->stable_id, 'rdfs:label', ($transcript->external_name)? '"'.$transcript->external_name.'"':'"'.$transcript->display_id.'"');

      # dump the features
      dump_feature($transcript);
      
      # if the transcript encodes a protein via translation 
      if ($transcript->translation) {
	  my $trans = $transcript->translation;
	  triple('ensembl:'.$transcript->stable_id, 'obo:SO_translates_to', 'ensembl:'.$trans->stable_id);
	  triple('ensembl:'.$trans->stable_id, 'rdfs:subClassOf', 'obo:SO_0000104');
	  triple('ensembl:'.$trans->stable_id, 'dc:identifier', '"'.$trans->stable_id.'"' );
	  triple('ensembl:'.$trans->stable_id, 'rdfs:label', '"'.$trans->display_id.'"');
	  taxonTriple('ensembl:'.$trans->stable_id);
      }
      
      # get the exons for the transcript
      my @exons = @{$transcript->get_all_Exons};
      # make all in the same strand
      #if ($transcript->strand == -1) {
	#  @exons = reverse @exons;
      #}
      my $position = 1;
      # Assert Exon bag for a given transcript, exons are ordered by position on the transcript.
      foreach my $exon (@exons) {
	  
	  # exon type of SO exon, both gene and transcript are linked via has part
	  triple('ensembl:'.$exon->stable_id,'rdfs:subClassOf','obo:SO_0000147');
	  triple('ensembl:'.$exon->stable_id,'dc:identifier', '"'.$exon->stable_id.'"');
	  triple('ensembl:'.$exon->stable_id, 'rdfs:label', '"'.$exon->display_id.'"');
	  taxonTriple('ensembl:'.$exon->stable_id);

	  # we assert that both the gene and the transcript have a part that is the exon
	  # The exon can refer to both the part on the gene and the part on the transcript 
	  triple('ensembl:'.$transcript->stable_id, 'obo:SO_has_part', 'ensembl:'.$exon->stable_id);
	  triple('ensembl:'.$gene->stable_id, 'obo:SO_has_part', 'ensembl:'.$exon->stable_id);
	  
	  # we add some triple to support querying for order of exons
	  triple('ensembl:'.$transcript->stable_id, 'sio:SIO_000974',  u($prefix{ensembl}.$transcript->stable_id.'#Exon_'.$position));
	  triple('ensembl:'.$transcript->stable_id.'#Exon_'.$position,  'rdfs:subClassOf', 'sio:SIO_001261');
	  triple('ensembl:'.$transcript->stable_id.'#Exon_'.$position, 'sio:SIO_000628', 'ensembl:'.$exon->stable_id);
	  triple('ensembl:'.$transcript->stable_id.'#Exon_'.$position, 'sio:SIO_000300', $position);
	  dump_feature($exon);
	  $position++;
      }
  }

  # Homology
  my $similar_genes = $gene->get_all_homologous_Genes;
  foreach my $alt_gene (map {$_->[0]} @$similar_genes) {
    triple('ensembl:'.$gene->stable_id, 'sio:SIO_000558', 'ensembl:'.$alt_gene->stable_id);
  }
  print STDERR ".";
  last if ($limit && $count == $limit);
}
print "Dumped triples for $count genes \n";

my %reference_hash;

sub dump_feature {
    my $feature = shift;
    
    my $slice = $feature->slice;
    my $region_name = $slice->seq_region_name;
    my $cs = $slice->coord_system;
    
    # generate a version specific portion of a URL that includes the database version, species, assembly version and region name
    # e.g. The URI for human chromosme 1 in ensmebl 75 on assembly GRCh37 would be http://rdf.ebi.ac.uk/resource/ensembl/75/homo_sapiens/GRCh37/1
    my $version_url = $prefix{ensembl}.$schemaVersion."/".$taxon_id. '/'.$cs->version.'/'.$feature->seq_region_name; 
    
    # we also create a non versioned URI that is a super class e.g. 
    # http://rdf.ebi.ac.uk/resource/ensembl/homo_sapiens/1
    my $non_version_url = $prefix{ensembl}.$taxon_id.'/'.$feature->seq_region_name; 

    # these are typed as chromosome or patches e.g. HG991_PATCH
    my $reference = u($version_url);
    if (! $reference_hash{$reference}) {
	triple($reference, 'rdfs:subClassOf', u($non_version_url));
	if ($cs->name eq 'chromosome') {	
	    triple(u($non_version_url), 'rdfs:subClassOf', 'obo:SO_0000340');
	}
	else {
	    triple(u($non_version_url), 'rdfs:subClassOf', 'term:'.$cs->name);
	    triple('term:'.$cs->name, 'rdfs:subClassOf', 'term:EnsemblRegion');
	}
	triple(u($non_version_url), 'rdfs:label', '"'.${species}.' '.$cs->name.' '.$feature->seq_region_name.'"');	
	triple($reference, 'rdfs:label', '"'.${species}.' '.$cs->name.' '.$feature->seq_region_name.' ('.$cs->version.')"');	
	triple($reference, 'dc:identifier', '"'.$feature->seq_region_name.'"');	
	triple($reference, 'term:inEnsemblSchemaNumber', '"'.$schemaVersion.'"');	
	triple($reference, 'term:inEnsemblAssembly', '"'.$cs->version.'"');	
	taxonTriple($reference);
	taxonTriple(u($non_version_url));

	$reference_hash{$reference} = 1;
    }

    # implement the FALDO model:  A semantic standard for describing the location of nucleotide and protein feature annotation
    # dx.doi.org/10.1101/002121
    my $begin = ($feature->strand == 1) ? $feature->start : $feature->end;
    my $end = ($feature->strand == 1) ? $feature->end : $feature->start;
    my $location = u($version_url.':'.$begin.'-'.$end.':'.$feature->strand);
    my $beginUri = u($version_url.':'.$begin.':'.$feature->strand);
    my $endUri = u($version_url.':'.$end.':'.$feature->strand);
    triple('ensembl:'.$feature->stable_id, 'faldo:location', $location);
    triple($location, 'rdfs:label', '"'.$cs->name.' '.$feature->seq_region_name.':'.$begin.'-'.$end.':'.$feature->strand.'"');
    triple($location, 'rdf:type', 'faldo:Region');
    triple($location, 'faldo:begin', $beginUri);
    triple($location, 'faldo:end', $endUri);
    triple($location, 'faldo:reference', $reference);
    triple($beginUri, 'rdf:type', 'faldo:ExactPosition');
    triple($beginUri, 'rdf:type', ($feature->strand == 1)? 'faldo:ForwardStrandPosition':'faldo:ReverseStrandPosition');
#    triple($beginUri, 'faldo:position', ($feature->strand == 1) ? $feature->start : $feature->end);
    triple($beginUri, 'faldo:position', $begin);
    triple($beginUri, 'faldo:reference', $reference);

    triple($endUri, 'rdf:type', 'faldo:ExactPosition');
    triple($endUri, 'rdf:type', ($feature->strand == 1)? 'faldo:ForwardStrandPosition':'faldo:ReverseStrandPosition');
#    triple($endUri, 'faldo:position', ($feature->strand == 1) ? $feature->end : $feature->start);
    triple($endUri, 'faldo:position', $end);
    triple($endUri, 'faldo:reference', $reference);

    triple('ensembl:'.$feature->stable_id, 'dc:identifier', '"'.$feature->stable_id.'"' );
}

sub dump_synonyms {

    my $feature = shift;
    
    my $db_entries = $feature->get_all_DBEntries();

    foreach my $dbe ( @{$db_entries} ) {
	foreach my $syn ( @{ $dbe->get_all_synonyms }) {	    
	    triple('ensembl:'.$feature->stable_id, 'skos:altLabel', '"'.$syn.'"');
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

sub taxonTriple {
    my $subject= shift;
    triple($subject, 'obo:RO_0002162', 'taxon:'.$taxon_id);
}

sub getSOOntologyId {
    my $term = shift;
    if ($term2ontologyId{$term}) {
	return $term2ontologyId{$term};
    }
    
    my ($typeterm) =
	@{ $ontoa->fetch_all_by_name( $term, 'SO' ) };
    
    if (!$typeterm) {
	print STDERR "WARN: Can't find SO term for $term\n";
	$term2ontologyId{$term} = "term:" . $term; 
	# make the biotype a child of an ensmebl biotype
	triple($term2ontologyId{$term}, 'rdfs:subClassOf', 'term:EnsemblFeature');
 	return $term2ontologyId{$term};
    }
    
    my $id = $typeterm->accession;
    $id=~s/SO:/obo:SO_/;
    $term2ontologyId{$term} = $id;

    # make the biotype a child of an ensmebl biotype
    triple($term2ontologyId{$term}, 'rdfs:subClassOf', 'term:EnsemblFeature');
 
    return $term2ontologyId{$term};
    
}

sub dumpVirtuoso {


}
close OUT;
