package uk.ac.ebi.pride.archive.indexer.services;

import de.mpc.pia.modeller.PIAModeller;
import de.mpc.pia.modeller.exporter.MzIdentMLExporter;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.stereotype.Service;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PIAModelerService;
import uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;

@Configuration
@Slf4j
@Service
public class PIAAnalysisService {

    private PIAModelerService piaModellerInference;

    @Value("${productionPath}")
    String productionPath;

    @Value("${qValueThreshold:#{0.01}}")
    private Double qValueThreshold;

    @Value("${qFilterProteinFDR:#{0.01}}")
    private Double qFilterProteinFDR;

    @Value("${peptideLength:#{7}}")
    private Double peptideLength;

    @Bean
    PIAModelerService getPiaModellerInference() {
        piaModellerInference = new PIAModelerService();
        return piaModellerInference;
    }

    public void performProteinInference(List<String> filePaths, SubmissionPipelineUtils.FileType fileType, String outputPath) throws IOException {
        PIAModeller modeller = piaModellerInference.performProteinInference(filePaths, fileType, qValueThreshold, qFilterProteinFDR);
        MzIdentMLExporter exporter = new MzIdentMLExporter(modeller);
        boolean exportOK = exporter.exportToMzIdentML(0L, new File(outputPath), false, true);
        if (exportOK)
            log.info("Analysis completed without errors.");
    }

    public void setqValueThreshold(Double qValueThreshold) {
        this.qValueThreshold = qValueThreshold;
    }

    public void setqFilterProteinFDR(Double qFilterProteinFDR) {
        this.qFilterProteinFDR = qFilterProteinFDR;
    }

    public void setPeptideLength(Double peptideLength) {
        this.peptideLength = peptideLength;
    }
}
