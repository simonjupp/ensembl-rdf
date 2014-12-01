use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';
# non verts  -host mysql-eg-publicsql.ebi.ac.uk -port 4157"
$registry->load_registry_from_multiple_dbs(
    #{
#	-host => 'mysql-ensembl-mirror.ebi.ac.uk', # ensembldb.ensembl.org alternativel 'useastdb.ensembl.org'
#	-port => 4240,
#	-user => 'anonymous'
#    },
    {
	-host => 'mysql-eg-publicsql.ebi.ac.uk',
	-port => 4157, 
	-user => 'anonymous'
    }
);

my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

foreach my $db_adaptor (@db_adaptors) {
    print "hello";
    printf("%s\n",$db_adaptor->species());
    if ( $db_adaptor->group() =~/core/) {
	
    }
}


