package uk.ac.ebi.pride.archive.indexer.services.ws;

/**
 * This code is licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * ==Overview==
 * <p>
 * This class
 * <p>
 * Created by ypriverol (ypriverol@gmail.com) on 02/07/2018.
 */
public interface PrideArchiveAPIField {
    String PRIDE_PROJECT_ACCESSION = "accession";
    String PRIDE_PROJECT_TITLE     = "title";
    String PRIDE_PROJECT_ADDITIONAL_ATTRIBUTES = "additionalAttributes";
    String PRIDE_PROJECT_DESCRIPTION =  "projectDescription";
    String PRIDE_PROJECT_SAMPLE_DESCRIPTION  =  "sampleProcessingProtocol";
    String PRIDE_PROJECT_DATA_PROCESSING =      "dataProcessingProtocol";
    String PRIDE_PROJECT_TAGS = "projectTags";
    String PRIDE_PROJECT_KEYWORDS = "keywords";
    String PRIDE_PROJECT_DOI = "doi";
    String PRIDE_PROJECT_SUBMISSION_DATE = "submissionDate";
    String PRIDE_PROJECT_PUBLICATION_DATE = "publicationDate";
    String PRIDE_PROJECT_UPDATED_DATE = "updatedDate";
    String PRIDE_PROJECT_SUBMITTERS = "submitters";
    String PRIDE_PROJECT_LAB_HEADS = "labPIs";
    String PRIDE_PROJECTS_AFFILIATIONS = "affiliations";
    String PRIDE_PROJECT_INSTRUMENTS = "instruments";
    String PRIDE_PROJECT_SOFTWARE = "softwares";
    String PRIDE_PROJECT_QUANT_METHODS = "quantificationMethods";
    String PRIDE_PROJECT_COUNTRY = "countries";
    String PRIDE_PROJECT_SAMPLE_ATTRIBUTES = "sampleAttributes";
    String PRIDE_PROJECT_ORGANISM = "organisms";
    String PRIDE_PROJECT_ORGANISM_PART = "organismParts";
    String PRIDE_PROJECT_DISEASES = "diseases";
    String PRIDE_PROJECT_REFERENCES = "references";
    String PRIDE_PROJECT_PTMS = "identifiedPTMS";
    String QUERY_SCORE = "queryScore";

}
