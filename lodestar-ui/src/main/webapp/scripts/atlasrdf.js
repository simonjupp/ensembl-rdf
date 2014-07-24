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

(function($) {

    $.fn.atlasmapping = function(options) {

        var mappingQuery = "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n" +
        "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>  \n" +
        "PREFIX owl: <http://www.w3.org/2002/07/owl#>  \n" +
        "PREFIX dcterms: <http://purl.org/dc/terms/>   \n" +
        "PREFIX obo: <http://purl.obolibrary.org/obo/>  \n" +
        "PREFIX efo: <http://www.ebi.ac.uk/efo/>     \n" +
        "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/>   \n" +
        "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/>  \n" +
        "PREFIX sio: <http://semanticscience.org/resource/>    \n" +

        "SELECT DISTINCT ?term ?label ?description ?relation ?relatedterm WHERE   \n" +
        "{           \n" +

        "?term a owl:Class .   \n" +
        "?term rdfs:label ?label .     \n" +
        "?term dcterms:description ?description .   \n" +
        " ?term ?relation ?relatedterm .   \n" +
        " ?relatedterm a owl:Class .     \n" +
        " FILTER regex(str(?term), \"^http://rdf.ebi.ac.uk/terms/atlas/\") .    \n" +
        " FILTER (! regex( str(?relatedterm), \"^http://rdf.ebi.ac.uk/terms/atlas/\")) . \n" +
        "}    \n" +
        "ORDER BY ?label";

    }
})

