

/*
 * Copyright (c) 2013 EMBL - European Bioinformatics Institute
 * Licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either express or implied. See the License for the specific language governing permissions
 * and limitations under the License.
 */

var exampleQueries = [

    {
        shortname : "Query 1",
        description: "Show all transcripts for human BRCA2 gene including location",
        query:
            "SELECT DISTINCT ?transcript ?id ?typeLabel ?referenceLabel ?begin ?end ?location {\n\n" +
            " # query human data only\n" +
            " GRAPH <http://rdf.ebi.ac.uk/dataset/ensembl/75/9606> {\n" +
            "   ?transcript obo:SO_transcribed_from ensembl:ENSG00000139618 ;\n" +
            "               rdfs:subClassOf ?type;\n" +
            "               dcterms:identifier ?id .\n" +
            " OPTIONAL {\n" +
            "  ?transcript faldo:location ?location .\n" +
            "  ?location faldo:begin [faldo:position ?begin] .\n" +
            "  ?location faldo:end [faldo:position ?end ] .\n" +
            "  ?location faldo:reference ?reference .\n" +
            "  ?reference rdfs:label ?referenceLabel\n" +
            "  } \n" +
            " }\n" +
            "OPTIONAL {?type rdfs:label ?typeLabel}\n" +
            "}"
    },
    {
        shortname : "Query 2",
        description: "Show ordered exons with their length for transcript ENST00000380152",
        query:
            "SELECT DISTINCT ?exon ?id ?referenceLabel ?begin ?end ?strand {\n\n" +
            " # query human data only\n" +
            " GRAPH <http://rdf.ebi.ac.uk/dataset/ensembl/75/9606> {\n" +
            "   ensembl:ENST00000380152 obo:SO_has_part ?exon;\n" +
            "                           sio:SIO_000974 ?orderedPart .\n" +
            "   ?exon dcterms:identifier ?id .\n" +
            "   # we include an explicit exon order\n" +
            "   # so that we can order correctly in both + and - strand \n" +
            "   ?orderedPart sio:SIO_000628 ?exon .\n" +
            "   ?orderedPart sio:SIO_000300 ?order .\n\n" +
            "   OPTIONAL {\n" +
            "     ?exon faldo:location ?location .\n" +
            "     ?location faldo:begin\n" +
            "        [a ?strand ;\n" +
            "           faldo:position ?begin] .\n" +
            "     ?location faldo:end\n" +
            "        [a ?strand ;\n" +
            "           faldo:position ?end] .\n" +
            "     ?location faldo:reference ?reference .\n" +
            "     ?reference rdfs:label ?referenceLabel\n" +
            "   }\n" +
            " }\n" +
            "FILTER (?strand != faldo:ExactPosition)\n" +
            "}\n" +
            "ORDER BY ASC(?order)"
     },
    {
        shortname : "Query 3",
        description: "Show proteins translated from ENST00000380152 that map to reviewed (SwissProt) and unreviewed (TrEMBL) UniProt entries",
        query:   "PREFIX core: <http://purl.uniprot.org/core/>\n" +
                 "SELECT DISTINCT ?id ?xrefRelationType ?xrefLabel ?xrefUri ?xrefType  {\n" +
                 " # query human data only\n" +
                 " GRAPH <http://rdf.ebi.ac.uk/dataset/ensembl/75/9606> {\n" +
                 "  ensembl:ENST00000380152 obo:SO_translates_to ?peptide .\n" +
                 "  ?peptide rdfs:subClassOf ?type ;\n" +
                 "           dcterms:identifier ?id ;\n" +
                 "           ?xrefRelationType ?xrefUri .\n" +
                 "  ?xrefUri rdfs:subClassOf ?xrefType ;\n" +
                 "           rdfs:label ?xrefLabel .\n" +
                 "  VALUES ?xrefType {core:Reviewed_Protein core:Protein}\n" +
                 " }\n" +
                 "}\n"
    },
    {
        shortname : "Query 4",
        description: "Show all external database cross references for gene ENSG00000139618",

        query:           "PREFIX skos: <http://www.w3.org/2004/02/skos/core#>\n" +
                         "SELECT DISTINCT ?id ?xrefRelationType ?xrefLabel ?xrefUri ?xrefType ?xrefLabel {\n" +
                         " # query human data only\n" +
                         " GRAPH <http://rdf.ebi.ac.uk/dataset/ensembl/75/9606> {\n" +
                         "  {\n" +
                         "    {?subject obo:SO_transcribed_from ensembl:ENSG00000139618}" +
                         "      UNION\n" +
                         "    {?transcript obo:SO_transcribed_from ensembl:ENSG00000139618 .\n" +
                         "     ?transcript obo:SO_translates_to  ?subject .}\n" +
                         "  }\n" +
                         "  ?subject rdfs:subClassOf ?type ;\n" +
                         "           dcterms:identifier ?id ;\n" +
                         "           ?xrefRelationType ?xrefUri .\n" +
                         "  ?xrefUri rdfs:subClassOf ?xrefType ;\n" +
                         "           rdfs:label ?xrefLabel .\n" +
                         "  ?xrefRelationType rdfs:subPropertyOf skos:related .\n" +
                         " }\n" +
                         "}\n"
    },
    {
        shortname : "Query 5",
        description: "Get all mouse genes on chromosome 11 between location 101,100,523 and 101,190,725 forward strand",
        query:
            "SELECT DISTINCT ?gene ?id ?label ?typelabel ?desc ?reference ?begin ?end {\n" +
            " GRAPH <http://rdf.ebi.ac.uk/dataset/ensembl/75/10090> {\n" +
            "  ?gene rdfs:subClassOf ?type ;\n" +
            "        rdfs:label ?label ;\n" +
            "        dcterms:description ?desc ;\n" +
           	"        dcterms:identifier ?id ;\n" +
            "        faldo:location ?location .\n" +
            "  ?location faldo:begin\n" +
            "       [a faldo:ForwardStrandPosition ;\n" +
            "        faldo:position ?begin] .\n" +
            "  ?location faldo:end\n" +
            "       [a faldo:ForwardStrandPosition ;\n" +
            "        faldo:position ?end] .\n" +
            "  ?location faldo:reference ?reference .\n" +
            "\n" +
            " # ensembl mouse chromosome 2 stable URI\n" +
            " ?reference rdfs:subClassOf <http://rdf.ebi.ac.uk/resource/ensembl/10090/11>\n" +
            " FILTER (?begin >= 101100523 && ?end <= 101190725 )\n" +
            " }\n" +
            " ?type rdfs:label ?typelabel .\n" +
            "}"
    },
    {
        shortname : "Query 6",
        description: "Get orthologs for human gene ENSG00000139618",
        query:
            "SELECT DISTINCT ?ortholog ?orthologLabel ?orthologSpecies {\n" +
            " ?gene sio:SIO_000558 ?ortholog .\n" +
            " ?gene obo:RO_0002162 [rdfs:label ?species] .\n" +
            " ?gene rdfs:label ?geneLabel .\n" +
            " ?ortholog rdfs:label ?orthologLabel .\n" +
            " ?ortholog obo:RO_0002162 [rdfs:label ?orthologSpecies] .\n" +
            " VALUES ?gene {ensembl:ENSG00000139618} \n" +
            " FILTER (?species != ?orthologSpecies) \n" +
            "}"
    },
    {
        shortname : "Query 7",
        description: "Get natural variant annotations and citation from UniProt for proteins encoded by ENSG00000139618 using a federated query",
        query:
                 "PREFIX core: <http://purl.uniprot.org/core/>\n" +
                 "PREFIX taxon:<http://purl.uniprot.org/taxonomy/>\n" +
                 "SELECT DISTINCT ?ensemblprotein ?xrefUri ?xrefLabel ?substitution ?text ?citation {\n" +
                 " # query human data only\n" +
                 " GRAPH <http://rdf.ebi.ac.uk/dataset/ensembl/75/9606> {\n" +
                 "  ?transcript obo:SO_transcribed_from ensembl:ENSG00000172936 .\n" +
                 "  ?transcript obo:SO_translates_to ?peptide .\n" +
                 "  ?peptide rdfs:subClassOf ?type ;\n" +
                 "           dcterms:identifier ?ensemblprotein ;\n" +
                 "           ?xrefRelationType ?xrefUri .\n" +
                 "  ?xrefUri rdfs:subClassOf core:Reviewed_Protein ;\n" +
                 "           rdfs:label ?xrefLabel .\n" +
                 "  }\n" +
                 "  SERVICE <http://beta.sparql.uniprot.org/sparql> {\n" +
                 "  	?xrefUri core:annotation ?annotation .\n" +
                 "  	?annotation a core:Natural_Variant_Annotation .\n" +
                 "  	?annotation rdfs:comment ?text .\n" +
                 "  	?annotation core:substitution ?substitution .\n" +
                 "  	?annotation core:range [faldo:begin [faldo:position ?location]] .\n" +
                 "  	?statement rdf:object ?annotation .\n" +
                 "  	?statement core:attribution ?ref .\n" +
                 "  	?ref core:source ?citation .\n" +
                 "  }\n" +
                 "}\n"
    },
    {
        shortname : "Query 8",
        description: "Show all species graphs loaded",
        query: "# all species data are loaded into different named graphs\n" +
               "# ontologies are also place in their own graph\n" +
               "# this query shows all the graphs available\n\n" +
               "SELECT ?g ?title WHERE {\n" +
               "?g <http://www.w3.org/2004/03/trix/rdfg-1/subGraphOf> <http://rdf.ebi.ac.uk/dataset/ensembl/75> . \n" +
               "?g dcterms:title ?title . \n" +
               "}"
    }


]