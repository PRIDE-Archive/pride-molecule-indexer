package uk.ac.ebi.pride.archive.pipeline.services.ws;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Configuration;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Service;
import org.springframework.web.client.RestTemplate;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParamProvider;

import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.Optional;

@Service
public class PrideArchiveWebService {

    RestTemplate restTemplate = new RestTemplate();

    @Value("${prideWSPublic}")
    String prideWSPublic;

    /**
     * Find the {@link PrideProject} for a specific project accession
     * @param projectAccession Project Accession
     * @return Optional {@link PrideProject}
     * @throws IOException
     */
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
                files = Arrays.asList(Objects.requireNonNull(result.getBody()));
            else
                throw new IOException("Connection to pride ws unuseful for project: " + projectAccession);
        } catch (URISyntaxException e) {
            throw new IOException("Connection to pride ws unuseful for project: " + projectAccession);
        }
        return files;
    }

    /**
     * Write the result files from a specific project into a tab-delimited file including the following properties:
     *  - Name file
     *  - date of publication
     *  - accession of the file
     *  - ftp location
     *
     * @param projectAccession project accession
     * @param fileOutput File output
     */
    public void writeResultForProjectAccession(String projectAccession, String fileOutput){
        try {
            Optional<PrideProject> projectOption = findByAccession(projectAccession);

            if(projectOption.isPresent()){
                String pattern = "yyyy-MM-dd";
                SimpleDateFormat simpleDateFormat = new SimpleDateFormat(pattern);

                String date = simpleDateFormat.format(projectOption.get().getPublicationDate());
                List<PrideFile> files = findFilesByProjectAccession(projectAccession);
                try (PrintWriter writer = new PrintWriter(
                        Files.newBufferedWriter(Paths.get(fileOutput)))) {
                    writer.printf("%s\t%s\t%s\t%s", "name", "date", "accession", "ftp");
                    writer.println();
                    files.forEach(x -> {
                        Optional<CvParamProvider> location = x.publicFileLocations.stream().filter(y -> Objects.equals(y.getAccession(), "PRIDE:0000469")).findFirst();
                        if(Objects.equals(x.getFileCategory().getValue(), "RESULT")){
                            writer.printf("%s\t%s\t%s\t%s", x.getFileName(),date, x.getAccession(), location.get().getValue());
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
