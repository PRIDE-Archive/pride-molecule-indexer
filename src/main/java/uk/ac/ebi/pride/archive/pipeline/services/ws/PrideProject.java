package uk.ac.ebi.pride.archive.pipeline.services.ws;

import com.fasterxml.jackson.annotation.JsonFormat;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonProperty;
import lombok.Data;
import uk.ac.ebi.pride.ws.pride.models.*;



import javax.xml.bind.annotation.XmlElement;
import java.util.*;

@JsonIgnoreProperties
@Data
public class PrideProject {

    @XmlElement
    @JsonProperty(PrideArchiveAPIField.PRIDE_PROJECT_ACCESSION)
    private String accession;

    @JsonProperty(PrideArchiveAPIField.PRIDE_PROJECT_TITLE)
    private String title;

    @JsonProperty(PrideArchiveAPIField.PRIDE_PROJECT_PUBLICATION_DATE)
    @JsonFormat(pattern="yyyy-MM-dd")
    private Date publicationDate;

    @Override
    public String toString() {
        return "PrideProject{" +
                "accession='" + accession + '\'' +
                ", title='" + title + '\'' +
                ", publicationDate=" + publicationDate +
                '}';
    }
}
