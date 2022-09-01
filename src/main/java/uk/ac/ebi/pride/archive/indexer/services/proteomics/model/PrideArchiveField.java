package uk.ac.ebi.pride.archive.indexer.services.proteomics.model;

/**
 * All properties used by the PRIDE Mongo Project, this values are used to perform queries in each service.
 *
 * @author ypriverol
 */
public interface PrideArchiveField {

    String ID = "id";
    /** Additional Attributes Accessions **/
    String ADDITIONAL_ATTRIBUTES = "additionalAttributes";
    /** Sample metadata names **/
    String SAMPLE_ATTRIBUTES_NAMES = "sample_attributes";
    /** Identified PTMs in the Project**/
    String PROJECT_IDENTIFIED_PTM = "ptmList";
    String PSM_SPECTRUM_ACCESSIONS = "psmAccessions";
    String PEPTIDE_SEQUENCE = "peptideSequence";
    String MODIFIED_PEPTIDE_SEQUENCE = "modifiedPeptideSequence";
    String PEPTIDOFORM = "peptidoform";
    String EXTERNAL_PROJECT_ACCESSION = "projectAccession";
    String EXTERNAL_REANALYSIS_ACCESSION = "reanalysisAccession";
    String PROTEIN_ASSAY_ACCESSION = "assayAccession";
    String IDENTIFICATION_DATABASE = "database";
    String CHARGE = "charge";
    String PRECURSOR_MASS = "precursorMass";
    String START_POSITION = "startPosition";
    String END_POSITION = "endPosition";
    String SCORES = "scores";
    String MISSED_CLEAVAGES = "missedCleavages";
    /** Alias Protein Table **/
    String PRIDE_PROTEIN_COLLECTION_NAME = "pride_protein_evidences";
    String PROTEIN_SEQUENCE = "proteinSequence";
    String UNIPROT_MAPPED_PROTEIN_ACCESSION = "uniprotMappedProteinAccession";
    String ENSEMBL_MAPPED_PROTEIN_ACCESSION = "ensemblMappedProteinAccession";
    String PROTEIN_GROUP_MEMBERS = "proteinGroupMembers";
    String PROTEIN_DESCRIPTION = "proteinDescription";
    String PROTEIN_MODIFICATIONS = "ptms";
    String IS_DECOY = "isDecoy";
    String BEST_SEARCH_ENGINE = "bestSearchEngineScore";
    String PROTEIN_REPORTED_ACCESSION = "reportedAccession";
    String PEPTIDE_ACCESSION = "peptideAccession";
    String PROTEIN_ACCESSION = "proteinAccession";
    String QUALITY_ESTIMATION_METHOD = "qualityEstimationMethods";
    String IS_VALIDATED = "isValid";
    String VALUE = "value";
    String NUMBER_PEPTIDEEVIDENCES = "numberPeptides";
    String NUMBER_PSMS = "numberPSMs";
    String PROTEIN_COVERAGE = "sequenceCoverage";
    String USI = "usi";
    String SPECTRA_USI = "spectraUsi";
    String PSM_SUMMARY_FILE = "fileName";
}
