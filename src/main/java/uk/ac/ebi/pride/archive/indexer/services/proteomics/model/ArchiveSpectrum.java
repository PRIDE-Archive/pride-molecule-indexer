package uk.ac.ebi.pride.archive.indexer.services.proteomics.model;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonRootName;
import com.fasterxml.jackson.annotation.JsonTypeName;
import lombok.Builder;
import lombok.Data;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PSMProvider;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModification;
import uk.ac.ebi.pride.archive.dataprovider.param.Param;
import uk.ac.ebi.pride.archive.dataprovider.param.ParamProvider;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@JsonRootName("ArchiveSpectrum")
@JsonTypeName("ArchiveSpectrum")
@Data
@Builder
@JsonIgnoreProperties(ignoreUnknown = true)
public class ArchiveSpectrum implements PSMProvider {

    @JsonProperty("usi")
    String usi;

    @JsonProperty("projectAccession")
    String projectAccession;

    @JsonProperty("assayAccession")
    String assayAccession;

    @JsonProperty("spectrumFile")
    String spectrumFile;

    @JsonProperty("sourceID")
    String sourceID;

    @JsonProperty("spectrumTitle")
    String spectrumTitle;

    @JsonProperty("masses")
    Double[] masses;

    @JsonProperty("intensities")
    Double[] intensities;

    @JsonProperty("numPeaks")
    Integer numPeaks;

    @JsonProperty("msLevel")
    Integer msLevel;

    @JsonProperty("precursorCharge")
    Integer precursorCharge;

    @JsonProperty("precursorMz")
    Double precursorMz;

    @JsonProperty("retentionTime")
    Double retentionTime;

    @JsonProperty("properties")
    Set<Param> properties;

    /** Interpretation of the Spectra **/

    @JsonProperty("peptideSequence")
    String peptideSequence;

    @JsonProperty("proteinAccessions")
    List<String> proteinAccessions;

    @JsonProperty("missedCleavages")
    Integer missedCleavages;

    @JsonProperty("modifications")
    Collection<IdentifiedModification> modifications;

    @JsonProperty("annotations")
    List<String> annotations;

    @JsonProperty("isDecoy")
    Boolean isDecoy;

    @JsonProperty("qualityEstimationMethods")
    private Set<Param> qualityEstimationMethods;

    @JsonProperty("isValid")
    private Boolean isValid;

    @JsonProperty("scores")
    private Set<Param> scores;

    @JsonProperty("bestSearhScore")
    Param bestSearchEngineScore;

    @JsonProperty("sample_properties")
    private Set<Param> sampleProperties;

    public ArchiveSpectrum() { }

    public ArchiveSpectrum(String usi, String projectAccession, String assayAccession, String spectrumFile, String sourceID, String spectrumTitle, Double[] masses, Double[] intensities, Integer numPeaks, Integer msLevel, Integer precursorCharge, Double precursorMz, Double retentionTime, Set<Param> properties, String peptideSequence, List<String> proteinAccessions, Integer missedCleavages, Collection<IdentifiedModification> modifications, List<String> annotations, Boolean isDecoy, Set<Param> qualityEstimationMethods, Boolean isValid, Set<Param> scores, Param bestSearchEngineScore, Set<Param> sampleProperties) {
        this.usi = usi;
        this.projectAccession = projectAccession;
        this.assayAccession = assayAccession;
        this.spectrumFile = spectrumFile;
        this.sourceID = sourceID;
        this.spectrumTitle = spectrumTitle;
        this.masses = masses;
        this.intensities = intensities;
        this.numPeaks = numPeaks;
        this.msLevel = msLevel;
        this.precursorCharge = precursorCharge;
        this.precursorMz = precursorMz;
        this.retentionTime = retentionTime;
        this.properties = properties;
        this.peptideSequence = peptideSequence;
        this.proteinAccessions = proteinAccessions;
        this.missedCleavages = missedCleavages;
        this.modifications = modifications;
        this.annotations = annotations;
        this.isDecoy = isDecoy;
        this.qualityEstimationMethods = qualityEstimationMethods;
        this.isValid = isValid;
        this.scores = scores;
        this.bestSearchEngineScore = bestSearchEngineScore;
        this.sampleProperties = sampleProperties;
    }

    @Override
    @JsonProperty(access = JsonProperty.Access.WRITE_ONLY)
    public Collection<? extends String> getAdditionalAttributesStrings() {
        return properties.stream().map(Param::getName).collect(Collectors.toList());
    }

    @Override
    @JsonProperty(access = JsonProperty.Access.WRITE_ONLY)
    public Collection<String> getModificationNames() {
        return modifications.stream().map(x -> x.getModificationCvTerm().getName()).collect(Collectors.toList());
    }

    @Override
    @JsonProperty(access = JsonProperty.Access.WRITE_ONLY)
    public Integer getNumberModifiedSites() {
        return modifications.size();
    }

    @Override
    public Boolean isDecoy() {
        return isDecoy;
    }

    @JsonProperty(access = JsonProperty.Access.WRITE_ONLY)
    public Collection<? extends ParamProvider> getAttributes() {
        return properties;
    }

    @Override
    public Boolean isValid() {
        return isValid;
    }
}
