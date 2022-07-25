package uk.ac.ebi.pride.archive.indexer.services.proteomics.model;


import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import lombok.Builder;
import lombok.Data;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PeptideSpectrumOverview;
import uk.ac.ebi.pride.archive.dataprovider.data.protein.ProteinDetailProvider;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModificationProvider;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParam;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParamProvider;
import uk.ac.ebi.pride.archive.dataprovider.param.Param;

import java.util.*;
import java.util.stream.Collectors;

/**
 * The {@link PrideProteinEvidence} contains the information of a protein identification in an specific experiment/dataset.
 * It contains the
 */

@Data
@Builder
@JsonIgnoreProperties(ignoreUnknown = true)
public class PrideProteinEvidence implements PrideArchiveField, ProteinDetailProvider {

    /** Generated accession **/
    @JsonProperty(PrideArchiveField.ID)
    @JsonIgnore
    private String id;

    /** Reported File ID is the Identifier of the File mzTab in PRIDE **/
    @JsonProperty(PrideArchiveField.PROTEIN_REPORTED_ACCESSION)
    private String reportedAccession;

    /** Accession in Reported File **/
    @JsonProperty(PrideArchiveField.PROTEIN_ASSAY_ACCESSION)
    private String assayAccession;

    /** External Project that contains the PSM **/
    @JsonProperty(EXTERNAL_PROJECT_ACCESSION)
    private String projectAccession;

    /** Uniprot protein identifier mapper **/
    @JsonProperty(PrideArchiveField.UNIPROT_MAPPED_PROTEIN_ACCESSION)
    private String uniprotMappedProteinAccession;

    /** Ensembl protein identifier mapper **/
    @JsonProperty(PrideArchiveField.ENSEMBL_MAPPED_PROTEIN_ACCESSION)
    private String ensemblMappedProteinAccession;

    /** Ensembl protein identifier mapper **/
    @JsonProperty(PrideArchiveField.PROTEIN_GROUP_MEMBERS)
    private Set<String> proteinGroupMembers;

    /** Ensembl protein identifier mapper **/
    @JsonProperty(PrideArchiveField.PROTEIN_DESCRIPTION)
    private String proteinDescription;

    /** Additional Attributes, Contains also scores**/
    @JsonProperty(PrideArchiveField.ADDITIONAL_ATTRIBUTES)
    private Set<Param> additionalAttributes;

    @JsonProperty(PrideArchiveField.PROTEIN_MODIFICATIONS)
    private List<String> ptms;

    @JsonProperty(PrideArchiveField.BEST_SEARCH_ENGINE)
    private Param bestSearchEngineScore;

    @JsonProperty(PrideArchiveField.SCORES)
    private Set<Param> scores;

    @JsonProperty(PrideArchiveField.SAMPLE_ATTRIBUTES_NAMES)
    private Set<Param> sampleProperties;

    @JsonProperty(PrideArchiveField.QUALITY_ESTIMATION_METHOD)
    private Set<Param> qualityEstimationMethods;

    @JsonProperty(PrideArchiveField.IS_VALIDATED)
    private Boolean isValid;

    @JsonProperty(PrideArchiveField.IS_DECOY)
    private boolean isDecoy;

    @JsonProperty(PrideArchiveField.NUMBER_PEPTIDEEVIDENCES)
    private Integer numberPeptides;

    @JsonProperty(PrideArchiveField.NUMBER_PSMS)
    private Integer numberPSMs;

    @JsonProperty(PrideArchiveField.PROTEIN_COVERAGE)
    private Double sequenceCoverage;

    @JsonProperty(PSM_SPECTRUM_ACCESSIONS)
    private Set<PeptideSpectrumOverview> psmAccessions;

    @Override
    @JsonIgnore
    public String getUniprotMapping() {
        return uniprotMappedProteinAccession;
    }

    @Override
    @JsonIgnore
    public String getEnsemblMapping() {
        return ensemblMappedProteinAccession;
    }

    @Override
    @JsonIgnore
    public Set<String> getProteinGroupMembers() {
        return proteinGroupMembers;
    }

    @Override
    public String getSubmittedSequence() {
        return null;
    }

    @JsonIgnore
    public Collection<String> getIdentifiedModifications() {
        return ptms;
    }

    @Override
    @JsonIgnore
    public String getAccession() {
        return reportedAccession;
    }

    @Override
    @JsonIgnore
    public String getDescription() {
        return proteinDescription;
    }

    @Override
    public Comparable getId() {
        return id;
    }

    @Override
    @JsonIgnore
    public Collection<? extends String> getAdditionalAttributesStrings() {
        return (additionalAttributes != null)? additionalAttributes
                .stream()
                .map(Param::getValue)
                .collect(Collectors.toList()): Collections.emptyList() ;
    }

    /**
     * This method add an attribute as {@link CvParamProvider} to the list of attributes.
     * @param attribute Attribute in {@link CvParamProvider}
     */
    public void addAttribute(CvParam attribute){
        if(additionalAttributes == null)
            additionalAttributes = new HashSet<>();
        additionalAttributes.add(new Param(attribute.getName(),
                attribute.getValue()));
    }

    /**
     * Add a to the list of modifications of a Protein.
     * @param modification {@link IdentifiedModificationProvider}
     */
    public void addIdentifiedModification(String modification){
        if(ptms == null)
            ptms = new ArrayList<>();
        ptms.add(modification);
    }

    @Override
    public String toString() {
        return "PrideMongoProteinEvidence{" +
                "id=" + id +
                ", reportedAccession='" + reportedAccession + '\'' +
                ", assayAccession='" + assayAccession + '\'' +
                ", projectAccession='" + projectAccession + '\'' +
                ", uniprotMappedProteinAccession='" + uniprotMappedProteinAccession + '\'' +
                ", ensemblMappedProteinAccession='" + ensemblMappedProteinAccession + '\'' +
                ", proteinGroupMembers=" + proteinGroupMembers +
                ", proteinDescription='" + proteinDescription + '\'' +
                ", additionalAttributes=" + additionalAttributes +
                ", ptms=" + ptms +
                ", bestSearchEngineScore=" + bestSearchEngineScore +
                ", scores=" + scores +
                ", sampleProperties=" + sampleProperties +
                ", qualityEstimationMethods=" + qualityEstimationMethods +
                ", isValid=" + isValid +
                ", isDecoy=" + isDecoy +
                ", numberPeptides=" + numberPeptides +
                ", numberPSMs=" + numberPSMs +
                ", sequenceCoverage=" + sequenceCoverage +
                '}';
    }
}
