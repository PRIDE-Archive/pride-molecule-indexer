package uk.ac.ebi.pride.archive.pipeline.services.ws;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Configuration;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;
import org.springframework.web.client.RestTemplate;

import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.Optional;

@Service
public class PrideArchiveWebService {

    RestTemplate restTemplate = new RestTemplate();

    @Value("${prideWSPublic}")
    String prideWSPublic;

    public Optional<PrideProject> findByAccession(String projectAccession) throws IOException {

        Optional<PrideProject> project;
        final String baseUrl = prideWSPublic + "projects/" + projectAccession;

        try {
            URI uri = new URI(baseUrl);
            ResponseEntity<PrideProject> result = restTemplate.getForEntity(uri, PrideProject.class);
            if (result.getStatusCode() == HttpStatus.OK && result.hasBody()){
                project = Optional.of(result.getBody());
            }else
                throw new IOException("Connection to pride ws unuseful for project: " + projectAccession);
        } catch (URISyntaxException e) {
            throw new IOException("Connection to pride ws unuseful for project: " + projectAccession);
        }
        return project;
    }

    public List<PrideFile> findFilesByProjectAccession(String projectAccession) throws IOException {

        final String baseUrl = prideWSPublic + "files/byProject?accession=" + projectAccession;
        List<PrideFile> files;

        try {
            URI uri = new URI(baseUrl);
            ResponseEntity<PrideFile[]> result = restTemplate.getForEntity(uri, PrideFile[].class);
            if (result.getStatusCode() == HttpStatus.OK && result.hasBody())
                files = Arrays.asList(result.getBody());
            else
                throw new IOException("Connection to pride ws unuseful for project: " + projectAccession);
        } catch (URISyntaxException e) {
            throw new IOException("Connection to pride ws unuseful for project: " + projectAccession);
        }
        return files;
    }

    public void writeResultForProjectAccession(String projectAccession, String fileOutput){
        try {
            Optional<PrideProject> projectOption = findByAccession(projectAccession);
            if(projectOption.isPresent()){
                List<PrideFile> files = findFilesByProjectAccession(projectAccession);
                try (PrintWriter writer = new PrintWriter(
                        Files.newBufferedWriter(Paths.get(fileOutput)))) {
                    files.forEach(x -> {
                        if(Objects.equals(x.getFileCategory().getValue(), "RESULT")){
                            writer.printf("%s\t%s\t", x.getFileName(), projectOption.get().getPublicationDate());
                            writer.println();
                        }
                    });
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
