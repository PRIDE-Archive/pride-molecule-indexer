package uk.ac.ebi.pride.archive.indexer.services.proteomics.model;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import lombok.Builder;
import lombok.Data;
import uk.ac.ebi.pride.archive.dataprovider.data.database.DatabaseProvider;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PeptideSequenceProvider;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PeptideSpectrumOverview;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModificationProvider;
import uk.ac.ebi.pride.archive.dataprovider.param.Param;
import uk.ac.ebi.pride.archive.dataprovider.param.ParamProvider;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * {@link PrideMongoPeptideEvidence} stores the information of a peptide related to a protein within one particular project/experiment.
 * Aprt of multiple properties such as Sequence, Protein accession, Sample properties etc. It contains a list of all the PSMs that support the
 * particular peptide evidence using a {@link Set } of {@link PeptideSpectrumOverview}. Properties are store as key value pairs.
 *
 * @author ypriverol
 *
 **/

@Data
@Builder
@JsonIgnoreProperties(ignoreUnknown = true)
public class PrideMongoPeptideEvidence implements PrideArchiveField, PeptideSequenceProvider {

    /** Generated accession **/
    @JsonIgnore
    private String id;

    /** Accession Provided by PRIDE Pipelines **/
    @JsonProperty(PrideArchiveField.PEPTIDE_ACCESSION)
    String peptideAccession;

    /** Reported File ID is the Identifier of the File mzTab in PRIDE **/
    @JsonProperty(PrideArchiveField.PROTEIN_ACCESSION)
    private String proteinAccession;

    /** Accession in Reported File **/
    @JsonProperty(PROTEIN_ASSAY_ACCESSION)
    private String assayAccession;

    /** External Project that contains the PSM **/
    @JsonProperty(EXTERNAL_PROJECT_ACCESSION)
    private String projectAccession;

    /** Peptide Sequence **/
    @JsonProperty(PrideArchiveField.PEPTIDE_SEQUENCE)
    private String peptideSequence;

    /** Database information used to perform the identification/quantification **/
    @JsonProperty(IDENTIFICATION_DATABASE)
    private DatabaseProvider database;

    /** PTMs Identified in the PEptide Sequence **/
    @JsonProperty(PROJECT_IDENTIFIED_PTM)
    private Collection<? extends IdentifiedModificationProvider> ptmList;

    /** Best Search engine scores **/
    @JsonProperty(BEST_SEARCH_ENGINE)
    ParamProvider bestSearchEngineScore;

    /** Best Search engine scores **/
    @JsonProperty(SCORES)
    Set<Param> scores;

    /** Sample properties **/
    @JsonProperty(SAMPLE_ATTRIBUTES_NAMES)
    Set<Param> sampleProperties;

    /** Additional Attributes **/
    @JsonProperty(PrideArchiveField.ADDITIONAL_ATTRIBUTES)
    private Set<Param> additionalAttributes;

    @JsonProperty(PSM_SPECTRUM_ACCESSIONS)
    private Set<PeptideSpectrumOverview> psmAccessions;

    @JsonProperty(IS_DECOY)
    private Boolean isDecoy;

    @JsonProperty(START_POSITION)
    private Integer startPosition;

    @JsonProperty(END_POSITION)
    private Integer endPosition;

    @JsonProperty(MISSED_CLEAVAGES)
    Integer missedCleavages;

    @JsonProperty(PrideArchiveField.QUALITY_ESTIMATION_METHOD)
    private Set<Param> qualityEstimationMethods;

    @JsonProperty(PrideArchiveField.IS_VALIDATED)
    private Boolean isValid;

    @Override
    @JsonIgnore
    public Collection<? extends IdentifiedModificationProvider> getModifications() {
        return null;
    }

    @Override
    @JsonIgnore
    public Collection<String> getModificationNames() {
        List<String> ptms = Collections.EMPTY_LIST;
        if(this.ptmList != null && !this.ptmList.isEmpty())
            ptms = ptmList.stream().map(x-> x.getModificationCvTerm().getName()).collect(Collectors.toList());
        return ptms;
    }

    @Override
    @JsonIgnore
    public Integer getNumberModifiedSites() {
        final int[] sites = {0};
        if(this.ptmList != null && !this.ptmList.isEmpty()){
            this.ptmList.forEach(x -> sites[0] += x.getPositionMap().size());
            return sites[0];
        }
        return 0;
    }

    @Override
    public Integer getMissedCleavages() {
        return missedCleavages;
    }

    @Override
    @JsonIgnore
    public Boolean isDecoy() {
        return isDecoy;
    }

    @Override
    @JsonIgnore
    public Collection<? extends String> getAdditionalAttributesStrings() {
        List<String> attributes = Collections.EMPTY_LIST;
        if(this.additionalAttributes != null )
            return additionalAttributes
                    .stream()
                    .map(Param::getValue)
                    .collect(Collectors.toList());
        return attributes;
    }




}
