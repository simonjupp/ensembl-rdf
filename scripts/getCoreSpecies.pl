use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-ensembl-mirror.ebi.ac.uk', # ensembldb.ensembl.orgalternatively 'useastdb.ensembl.org'
    -port => 4240,
    -user => 'anonymous'
    );

my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

foreach my $db_adaptor (@db_adaptors) {
    my $db_connection = $db_adaptor->dbc();

    if ( $db_adaptor->group() =~/core/) {
	printf("%s\n",$db_adaptor->species());
    }
}
