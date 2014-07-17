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
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


# Choice of database host is a factor in how fast the script runs. Try to find your nearest mirror, and check the database version before running.
Bio::EnsEMBL::Registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
  -db_version => 75,
  -no_cache => 1,
);


my %prefix = (
  base => 'http://rdf.ebi.ac.uk/resource/ensembl/',
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
my $ga = Bio::EnsEMBL::Registry->get_adaptor('Human','Core','Gene');

sub print_DBEntries
{
    my $db_entries = shift;
    
    foreach my $dbe ( @{$db_entries} ) {
        printf "\tXREF %s (db:%s) (id:%s) (desc:%s) (enstype:%s) (exttype:%s) (evidence:%s) (linkage:%s)\n",
	$dbe->display_id(), 
	$dbe->dbname(),
	$dbe->primary_id(),
	$dbe->description(),
	$dbe->ensembl_object_type(),
	$dbe->type(),
	$dbe->info_type(),
	$dbe->linkage_annotation()
	;
    }
}

my $gene = $ga->fetch_by_stable_id('ENSG00000139618');

print "GENE ", $gene->stable_id(), "\n";
print_DBEntries( $gene->get_all_DBEntries() );

foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
    print "TRANSCRIPT ", $transcript->stable_id(), "\n";
    print_DBEntries( $transcript->get_all_DBEntries() );
    
# Watch out: pseudogenes have no translation
    if ( defined $transcript->translation() ) {
	my $translation = $transcript->translation();
	
	print "TRANSLATION ", $translation->stable_id(), "\n";
	print_DBEntries( $translation->get_all_DBEntries() );
    }
}

sub u {
    my $stuff= shift;
    return '<'.$stuff.'>';
}
sub triple {
    my ($subject,$predicate,$object) = @_;
    
    printf "%s %s %s .\n",$subject,$predicate,$object;
}

