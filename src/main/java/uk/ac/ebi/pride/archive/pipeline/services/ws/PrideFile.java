package uk.ac.ebi.pride.archive.pipeline.services.ws;

import com.fasterxml.jackson.annotation.JsonFormat;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import lombok.Data;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParamProvider;

import java.util.Date;
import java.util.List;

@JsonIgnoreProperties
@Data
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
