package uk.ac.ebi.pride.archive.indexer.services;

import lombok.extern.slf4j.Slf4j;
import org.springframework.stereotype.Service;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.MGFPRIDEWriter;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PrideJsonRandomAccess;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Paths;

@Slf4j
@Service
public class PSMClusteringService {

    public PSMClusteringService() {
    }

    public void convertToMgf(String prideJsonFile, String mgfOutputFile) {
        try {
            PrideJsonRandomAccess pridePSMJsonReader = new PrideJsonRandomAccess(prideJsonFile);
            pridePSMJsonReader.parseIndex();
            OutputStream outputStream       = Files.newOutputStream(Paths.get(mgfOutputFile));
            OutputStreamWriter outputStreamWriter = new OutputStreamWriter(outputStream);
            for(String usi: pridePSMJsonReader.getKeys())
                MGFPRIDEWriter.appendSpectrum(outputStreamWriter, pridePSMJsonReader.readArchiveSpectrum(usi));

        } catch (IOException e) {
            throw new RuntimeException(e);
        }


    }
}
