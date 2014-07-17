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

# Author: Kieron Taylor, Simon Jupp, European Bioinformatics Institute

# Script to dump Ensembl triples. It was lashed together rapidly, and could stand to use a library for the writing of triples, with accommodation for the scale of the data involved.
# Requires installation of Ensembl Core and Compara APIs, as well as dependencies such as BioPerl

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


# Choice of database host is a factor in how fast the script runs. Try to find your nearest mirror, and check the database version before running.
Bio::EnsEMBL::Registry->load_registry_from_db(
  -host => 'mysql-ensembl-mirror.ebi.ac.uk',
  -user => 'anonymous',
  -port => 4240,
  -db_version => 75,
  -no_cache => 1,
);

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
my $ga = Bio::EnsEMBL::Registry->get_adaptor('Human','Core','Gene');

my $ontoa =
    Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

# get all genes
my $genes = $ga->fetch_all;

# create a lookup table for term to ontology id
my %term2ontologyId;

# simple count for testing
my $count = 0;

my $meta = Bio::EnsEMBL::Registry->get_adaptor('Human','Core','MetaContainer');

# create a map of taxon id to species name, we will create some triples for these at the end

# get the taxon id for this species 
my $taxon_id = $meta->get_taxonomy_id;
my $scientific_name = $meta->get_scientific_name;

# print out global triples about the organism  
triple('taxon:'.$taxon_id, 'rdfs:subClassOf', 'obo:OBI_0100026');
triple('taxon:'.$taxon_id, 'rdfs:label', '"'.$scientific_name.'"');

# get the current ensembl release number 
my $schemaVersion = $meta->get_schema_version;

# start to process all genes
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

  # dump all features
  dump_feature($gene);

  # add some useful meta data
  triple('ensembl:'.$gene->stable_id, 'rdfs:label', ($gene->external_name)? '"'.$gene->external_name.'"':'"'.$gene->display_id.'"' );

  triple('ensembl:'.$gene->stable_id, 'dc:description', '"'.$gene->description.'"');
# we will get all synonyms as mappings using the mapings script
# triple('ensembl:'.$gene->stable_id, 'skos:altLabel', '"'.$gene->external_name.'"');
  
  # relate gene to its taxon
  taxonTriple('ensembl:'.$gene->stable_id);

  # loop through the transcripts
  foreach my $transcript (@trans) {
      # transcripts are transcribed from a gene 
      triple( 'ensembl:'.$transcript->stable_id, 'obo:SO_transcribed_from', 'ensembl:'.$gene->stable_id);
      taxonTriple('ensembl:'.$transcript->stable_id);
      
      # type the transcript using SO again
      my $transcriptTypeId = getSOOntologyId($transcript->biotype());
      triple( 'ensembl:'.$transcript->stable_id, 'rdfs:subClassOf', $transcriptTypeId);
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
      if ($transcript->strand == -1) {
	  @exons = reverse @exons;
      }
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
  print STDERR ".\n";
  last if ($count == 100);
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
    my $version_url = $prefix{ensembl}.$schemaVersion."/".$cs->species. '/'.$cs->version.'/'.$feature->seq_region_name; 
    
    # we also create a non versioned URI that is a super class e.g. 
    # http://rdf.ebi.ac.uk/resource/ensembl/homo_sapiens/1
    my $non_version_url = $prefix{ensembl}.$cs->species.'/'.$feature->seq_region_name; 

    # these are typed as chromosome or patches e.g. HG991_PATCH
    my $reference = u($version_url);
    if (! $reference_hash{$reference}) {
	triple($reference, 'rdfs:subClassOf', u($non_version_url));
	if ($cs->name eq 'chromosome') {	
	    triple(u($non_version_url), 'rdfs:subClassOf', 'obo:SO_0000340');
	}
	else {
	    triple(u($non_version_url), 'rdfs:subClassOf', 'ensembl:'.$cs->name);
	}
	triple(u($non_version_url), 'rdfs:label', '"'.$cs->name.' '.$feature->seq_region_name.'"');	
	triple($reference, 'rdfs:label', '"'.$cs->name.' '.$feature->seq_region_name.'"');	
	triple($reference, 'dc:identifier', '"'.$feature->seq_region_name.'"');	
	triple($reference, 'ensembl:inEnsemblSchemaNumber', '"'.$schemaVersion.'"');	
	triple($reference, 'ensembl:inEnsemblAssembly', '"'.$cs->version.'"');	
	taxonTriple($reference);
	taxonTriple(u($non_version_url));

	$reference_hash{$reference} = 1;
    }

    # implement the FALDO model:  A semantic standard for describing the location of nucleotide and protein feature annotation
    # dx.doi.org/10.1101/002121
    my $location = u($version_url.':'.$feature->start.'-'.$feature->end.':'.$feature->strand);
    my $begin = u($version_url.':'.$feature->start.':'.$feature->strand);
    my $end = u($version_url.':'.$feature->end.':'.$feature->strand);
    triple('ensembl:'.$feature->stable_id, 'faldo:location', $location);
    triple($location, 'rdf:type', 'faldo:Region');
    triple($location, 'faldo:begin', $begin);
    triple($location, 'faldo:end', $end);
    triple($begin, 'rdf:type', 'faldo:ExactPosition');
    triple($begin, 'rdf:type', ($feature->strand == 1)? 'faldo:ForwardStrandPosition':'faldo:ReverseStrandPosition');
    triple($begin, 'faldo:position', ($feature->strand == 1) ? $feature->start : $feature->end);
    triple($begin, 'faldo:reference', $reference);

    triple($end, 'rdf:type', 'faldo:ExactPosition');
    triple($end, 'rdf:type', ($feature->strand == 1)? 'faldo:ForwardStrandPosition':'faldo:ReverseStrandPosition');
    triple($end, 'faldo:position', ($feature->strand == 1) ? $feature->end : $feature->start);
    triple($end, 'faldo:reference', $reference);

    triple('ensembl:'.$feature->stable_id, 'dc:identifier', '"'.$feature->stable_id.'"' );
}

sub u {
    my $stuff= shift;
    return '<'.$stuff.'>';
}
sub triple {
    my ($subject,$predicate,$object) = @_;
    
    printf "%s %s %s .\n",$subject,$predicate,$object;
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
	$term2ontologyId{$term} = "obo:" . $term; 
	return $term2ontologyId{$term};
    }
    
    my $id = $typeterm->accession;
    $id=~s/SO:/obo:SO_/;
    $term2ontologyId{$term} = $id; 
    return $term2ontologyId{$term};
    
}
