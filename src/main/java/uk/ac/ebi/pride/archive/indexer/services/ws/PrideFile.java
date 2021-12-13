package uk.ac.ebi.pride.archive.indexer.services.ws;

import com.fasterxml.jackson.annotation.JsonFormat;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParamProvider;

import java.util.Date;
import java.util.List;

@JsonIgnoreProperties
public class PrideFile {

    @JsonProperty
    List<String> projectAccessions;

    @JsonProperty
    String accession;
    
    @JsonProperty
    String fileName;

    @JsonProperty
    List<CvParamProvider> publicFileLocations;

    @JsonFormat(pattern="yyyy-MM-dd")
    @JsonProperty
    private Date publicationDate;

    @JsonProperty
    private CvParamProvider fileCategory;

    public List<String> getProjectAccessions() {
        return projectAccessions;
    }

    public String getAccession() {
        return accession;
    }

    public String getFileName() {
        return fileName;
    }

    public List<CvParamProvider> getPublicFileLocations() {
        return publicFileLocations;
    }

    public Date getPublicationDate() {
        return publicationDate;
    }

    public CvParamProvider getFileCategory() {
        return fileCategory;
    }

    @Override
    public String toString() {
        return "PrideFile{" +
                "projectAccessions=" + projectAccessions +
                ", accession='" + accession + '\'' +
                ", fileName='" + fileName + '\'' +
                ", publicFileLocations=" + publicFileLocations +
                ", publicationDate=" + publicationDate +
                ", fileCategory=" + fileCategory +
                '}';
    }
}
