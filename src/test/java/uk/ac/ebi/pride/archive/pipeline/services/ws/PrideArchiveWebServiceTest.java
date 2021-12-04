package uk.ac.ebi.pride.archive.pipeline.services.ws;

import lombok.extern.slf4j.Slf4j;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.TestPropertySource;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideArchiveWebService;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideFile;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideProject;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(classes = {PrideArchiveWebService.class})
@TestPropertySource(value = "classpath:application.yaml")
@Slf4j
public class PrideArchiveWebServiceTest {

    String projectAccession = "PXD029360";

    @Autowired
    PrideArchiveWebService webService;

    @Test
    public void findByAccession() throws IOException {
        PrideProject project = webService.findByAccession(projectAccession).get();
        assertEquals(project.getAccession(), projectAccession);
        log.info(project.toString());
    }

    @Test
    public void findFilesByProjectAccession() throws IOException {
        List<PrideFile> files = webService.findFilesByProjectAccession(projectAccession);
        assertEquals(files.size(), 1138);
        log.info(files.toString());
    }
}