package uk.ac.ebi.pride.archive.indexer.services.proteomics.model;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import lombok.Builder;
import lombok.Data;
import uk.ac.ebi.pride.archive.dataprovider.param.Param;

import java.util.List;
import java.util.Set;


@Data
@Builder
@JsonIgnoreProperties(ignoreUnknown = true)
public class PridePsmSummaryEvidence implements PrideArchiveField{

    /** Generated accession **/
    @JsonProperty(PrideArchiveField.ID)
    @JsonIgnore
    private String id;

    /** Accession Provided by PRIDE Pipelines **/
    @JsonProperty(PrideArchiveField.USI)
    String usi;

    /** Accession Provided by PRIDE Pipelines **/
    @JsonProperty(PrideArchiveField.SPECTRA_USI)
    String spectraUsi;

    /** Accession in Reported File **/
    @JsonProperty(PrideArchiveField.PROTEIN_ASSAY_ACCESSION)
    private String assayAccession;

    /** External Project that contains the PSM **/
    @JsonProperty(PrideArchiveField.EXTERNAL_PROJECT_ACCESSION)
    private String projectAccession;

    /** External Project that contains the PSM **/
    @JsonProperty(PrideArchiveField.EXTERNAL_REANALYSIS_ACCESSION)
    private String reanalysisAccession;

    @JsonProperty("proteinAccessions")
    List<String> proteinAccessions;

    /** Peptide Sequence **/
    @JsonProperty(PrideArchiveField.PEPTIDE_SEQUENCE)
    private String peptideSequence;

    /** Modified  Sequence **/
    @JsonProperty(PrideArchiveField.MODIFIED_PEPTIDE_SEQUENCE)
    private String modifiedPeptideSequence;

    /** Additional Attributes **/
    @JsonProperty(PrideArchiveField.SCORES)
    private Set<Param> scores;

    /** Additional Attributes **/
    @JsonProperty(PrideArchiveField.SAMPLE_ATTRIBUTES_NAMES)
    private Set<Param> sampleProperties;

    @JsonProperty(PrideArchiveField.IS_DECOY)
    private Boolean isDecoy;

    @JsonProperty(PrideArchiveField.IS_VALIDATED)
    private Boolean isValid;

    @JsonProperty(PrideArchiveField.CHARGE)
    private Integer charge;

    @JsonProperty(PrideArchiveField.PRECURSOR_MASS)
    private Double precursorMass;

    @JsonProperty(PrideArchiveField.PSM_SUMMARY_FILE)
    private String fileName;

    @JsonProperty(PrideArchiveField.BEST_SEARCH_ENGINE)
    private Param bestSearchEngineScore;

}
