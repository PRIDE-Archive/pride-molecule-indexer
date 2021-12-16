package uk.ac.ebi.pride.archive.indexer.services.proteomics.model;

/**
 * All properties used by the PRIDE Mongo Project, this values are used to perform queries in each service.
 *
 * @author ypriverol
 */
public interface PrideArchiveField {

    String ID = "id";

    /** Project Accession **/
    String ACCESSION = "accession";

    /** Assay **/
    String ASSAY = "assayType";

    /** Project Accession **/
    String PX_ACCESSION = "pxAccession";

    /** Project Title **/
    String PROJECT_TILE = "title";

    /** RELATIONS BETWEEN FILES **/
    String FILE_RELATIONS_IN_PROJECT = "submittedFileRelations";

    String ACCESSION_SUBMISSION_FILE = "accessionSubmissionFile";

    /** Additional Attributes Accessions **/
    String ADDITIONAL_ATTRIBUTES = "additionalAttributes";

    /** Project Description **/
    String PROJECT_DESCRIPTION = "description";

    String SAMPLE = "sample";
    String SAMPLE_FILE_CHECKSUM = "filechecksum";

    /** Sample Protocol **/
    String PROJECT_SAMPLE_PROTOCOL = "sampleProtocol";

    /** Data Processing Protocol **/
    String PROJECT_DATA_PROTOCOL = "dataProtocol";

    /** Project Tags **/
    String PROJECT_TAGS = "tags";

    /** Keywords **/
    String PROJECT_KEYWORDS = "keywords";

    /** PROJECT DOI **/
    String PROJECT_DOI = "projectDoi";

    /** PROJECT OMICS **/
    String PROJECT_OMICS_LINKS = "project_other_omics";

    /** Submission Type **/
    String PROJECT_SUBMISSION_TYPE = "submissionType";

    /** Submission Date **/
    String SUBMISSION_DATE = "submissionDate";

    /** Publication Date **/
    String PUBLICATION_DATE = "publicationDate";

    /** Update Date **/
    String UPDATED_DATE = "updatedDate";

    /** Submitter FirstName **/
    String PROJECT_SUBMITTER = "submitters";

    /** List of Lab Head Names **/
    String PROJECT_PI_NAMES = "lab_heads";

    /** List of instruments Ids*/
    String INSTRUMENTS = "instruments";

    String SOFTWARES = "softwares";

    String QUANTIFICATION_METHODS = "quantificationMethods";

    /** This field store all the countries associated with the experiment **/
    String COUNTRIES = "countries";

    /** Sample metadata names **/
    String SAMPLE_ATTRIBUTES_NAMES = "sample_attributes";

    /** References related with the project **/
    String PROJECT_REFERENCES = "project_references";

    /** Identified PTMs in the Project**/
    String PROJECT_IDENTIFIED_PTM = "ptmList";

    /** Collections Name **/
    String PRIDE_PROJECTS_COLLECTION_NAME = "pride_projects";
    String PRIDE_EXPERIMENTAL_DESIGN_COLLECTION_NAME = "pride_experimental_design";
    String PRIDE_FILE_COLLECTION_NAME = "pride_files";
    String PRIDE_MSRUN_COLLECTION_NAME = "pride_msruns";
    String PRIDE_REANALYSIS_COLLECTION_NAME = "pride_reanalysis_collection";
    String PRIDE_MOLECULES_COLLECTION_NAME = "pride_molecules";
    String PRIDE_SDRF_COLLECTION_NAME = "pride_sdrf";

    String PUBLIC_PROJECT = "public_project";

    /** Experimental Factors **/
    String EXPERIMENTAL_FACTORS = "experimentalFactors";

    /** External Project accessions that use this following file **/
    String EXTERNAL_PROJECT_ACCESSIONS = "projectAccessions";

    /** External Project Analysis Accessions that use the file **/
    String EXTERNAL_ANALYSIS_ACCESSIONS = "analysisAccessions";

    //** File Fields **/
    String FILE_CATEGORY = "fileCategory";
    String FILE_SOURCE_TYPE = "fileSourceType";
    String FILE_SOURCE_FOLDER = "fileSourceFolder";
    String FILE_CHECKSUM = "fileChecksum";
    String CHECKSUM = "checksum";
    String SUBMITTER_FILE_CHECKSUM = "submitterFileChecksum";
    String FILE_PUBLIC_LOCATIONS = "filePublicLocations";
    String FILE_SIZE_MB = "fileSizeMB";
    String FILE_EXTENSION = "fileExtension";
    String FILE_NAME = "fileName";
    String FILE_IS_COMPRESS = "fileCompress";

    /** PSM Collections **/
    String PRIDE_PEPTIDE_COLLECTION_NAME = "pride_peptide_evidences";
    String PEPTIDE_SUMMARY_COLLECTION_NAME = "peptide_summary";
    String PRIDE_PSM_COLLECTION_NAME = "pride_psm_evidences";
    String PSM_SPECTRUM_ACCESSIONS = "psmAccessions";
    String PEPTIDE_SEQUENCE = "peptideSequence";
    String MODIFIED_PEPTIDE_SEQUENCE = "modifiedPeptideSequence";
    String REPORTED_PROTEIN_ACCESSION = "reportedProteinAccession";
    String REPORTED_FILE_ID = "reportedFileID";
    String EXTERNAL_PROJECT_ACCESSION = "projectAccession";
    String EXTERNAL_ANALYSIS_ACCESSION = "analysisAccession";
    String PROTEIN_ASSAY_ACCESSION = "assayAccession";
    String IDENTIFICATION_DATABASE = "database";
    String PEPTIDE_UNIQUE = "peptideUnique";
    String BEST_PSM_SCORE = "bestPSMScore";
    String RETENTION_TIME = "retentionTime";
    String CHARGE = "charge";
    String PRECURSOR_MASS = "precursorMass";
    String EXPERIMENTAL_MASS_TO_CHARGE = "expMassToCharge";
    String CALCULATED_MASS_TO_CHARGE = "calculatedMassToCharge";
    String PRE_AMINO_ACID = "preAminoAcid";
    String POST_AMINO_ACID  = "postAminoAcid";
    String START_POSITION = "startPosition";
    String END_POSITION = "endPosition";
    String SEARCH_ENGINE_SCORES = "searchEngineScores";
    String SCORES = "scores";
    String MISSED_CLEAVAGES = "missedCleavages";

    /** PRIDE Analysis Collection **/
    String PRIDE_ANALYSIS_COLLECTION = "pride_analysis_collection";

    /** PRIDE Stats Fields*/
    String PRIDE_STATS_COLLECTION = "pride_stats_collection";
    String PEPTIDOME_COLLECTION = "peptidome_stats";
    String STATS_SUBMISSION_COUNTS = "pride_submission_counts";
    String STATS_ESTIMATION_DATE = "estimationDate";
    String STATS_COMPLEX_COUNTS = "pride_complex_counts";

    /** MSRun file Properties **/
    String MS_RUN_FILE_PROPERTIES = "FileProperties";
    String MS_RUN_INSTRUMENT_PROPERTIES = "InstrumentProperties";
    String MS_RUN_MS_DATA = "MsData";
    String MS_RUN_SCAN_SETTINGS = "ScanSettings";
    String MS_RUN_ID_SETTINGS = "IdSettings";

    /*MSRun ID Settigns fields*/
    String MS_RUN_ID_SETTINGS_FIXED_MODIFICATIONS = "FixedModifications";
    String MS_RUN_ID_SETTINGS_VARIABLE_MODIFICATIONS = "VariableModifications";
    String MS_RUN_ID_SETTINGS_ENZYMES = "Enzymes";
    String MS_RUN_ID_SETTINGS_FRAGMENT_TOLERANCE = "FragmentTolerance";
    String MS_RUN_ID_SETTINGS_PARENT_TOLERANCE = "ParentTolerance";



    /** Alias **/
    String MONGO_MSRUN_DOCUMENT_ALIAS = "MongoPrideMSRun";
    String MONGO_FILE_DOCUMENT_ALIAS = "MongoPrideFile";

    String SAMPLES = "samples";
    String SAMPLES_MSRUN = "samples_msrun";

    /** ASSAY Properties **/
    String PRIDE_ASSAY_COLLECTION_NAME = "pride_assays";
    String ASSAY_ACCESSION = "accession";
    String ASSAY_FILE_NAME = "fileName";
    String ASSAY_TITLE = "assayTitle";
    String ASSAY_DESCRIPTION = "assayDescription";
    String ASSAY_DATA_ANALYSIS_SOFTWARE = "dataAnalysisSoftwares";
    String ASSAY_DATA_ANALYSIS_DATABASE = "dataAnalysisDatabase";
    String ASSAY_DATA_ANALYSIS_RESULTS = "summaryResults";
    String ASSAY_DATA_ANALYSIS_PROTOCOL = "dataAnalysisProperties";
    String ASSAY_DATA_ANALYSIS_PTMS = "ptmsResults";
    String ASSAY_FILES = "assayFiles";
    String VALID_ASSAY = "validAssay";

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

    String MSRUN_PROPERTIES = "MSRunProperties";

    String PEPTIDE_ACCESSION = "peptideAccession";
    String PROTEIN_ACCESSION = "proteinAccession";
    String PEPTIDE_ACCESSIONS = "peptideAccessions";
    String PROTEIN_ACCESSIONS = "proteinAccessions";
    String QUALITY_ESTIMATION_METHOD = "qualityEstimationMethods";
    String IS_VALIDATED = "isValid";
    String VALUE = "value";

    String NUMBER_PEPTIDEEVIDENCES = "numberPeptides";
    String NUMBER_PSMS = "numberPSMs";
    String PROTEIN_COVERAGE = "sequenceCoverage";

    String USI = "usi";
    String SPECTRA_USI = "spectraUsi";
    String PSM_SUMMARY_FILE = "fileName";

    String BEST_SEARCH_ENGINE_SCORE = "best_search_engine_score";
    String PSMS_COUNT = "psms_count";
    String BEST_USIS = "best_usis";
    String PTMS_MAP = "ptms_map";
    String IS_UNIPROT = "is_uniprot_accession";
    String IS_UNIQUE  = "is_uniq_peptide_within_organism";
    String IS_MULTI_ORGANISM = "is_multiorganism_peptides";
    String TAX_ID = "TaxId";
    String ORGANISM = "Organism";
    String GENE = "Gene";
    String PROTEIN_NAME = "ProteinName";
    String EVIDENCE_LEVEL = "EvidenceLevel";

}
