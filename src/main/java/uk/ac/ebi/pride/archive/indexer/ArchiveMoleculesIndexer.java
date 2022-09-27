package uk.ac.ebi.pride.archive.indexer;

import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.ConfigurableApplicationContext;
import org.springframework.context.annotation.ComponentScan;
import uk.ac.ebi.pride.archive.indexer.services.InferenceService;
import uk.ac.ebi.pride.archive.indexer.services.PSMClusteringService;
import uk.ac.ebi.pride.archive.indexer.services.PrideAnalysisAssayService;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideArchiveWebService;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

@Slf4j
@SpringBootApplication
@ComponentScan(basePackages = "uk.ac.ebi.pride.archive.indexer.services")
public class ArchiveMoleculesIndexer implements ApplicationRunner {

    private final String[] options = {"get-result-files", "get-related-files",
            "generate-index-files", "perform-inference",
            "generate-mgf-files"};

    @Autowired
    private ConfigurableApplicationContext context;

    @Autowired
    private PrideArchiveWebService prideArchiveWebService;

    @Autowired
    PrideAnalysisAssayService analysisAssayService;

    @Autowired
    InferenceService inferenceAnalysisService;

    @Autowired
    PSMClusteringService clusteringService;

    /**
     * usage outputDirectory filterFile directoryToProcess
     * @param args
     */
    public static void main(String[] args){
        SpringApplication.run(ArchiveMoleculesIndexer.class, args);
    }

    @Override
    public void run(ApplicationArguments args) throws Exception {
        if (args.getNonOptionArgs().size() != 1 || !Arrays.asList(options).contains(args.getNonOptionArgs().get(0))){
            throw new Exception("Available commands are: " + Arrays.asList(options));
        }
        // Get results files command
        String command = args.getNonOptionArgs().get(0);
        if(Objects.equals(command, "get-result-files")){
            List<String> projectAccessionOptions = args.getOptionValues("app.project-accession");
            if(projectAccessionOptions.size() != 1){
                throw new Exception("Project accession should be provided for command " +
                        "get-result-files with parameter --app.project-accession");
            }
            String projectAccession = projectAccessionOptions.get(0);

            List<String> fileOutputOptions = args.getOptionValues("app.file-output");
            if(fileOutputOptions.size() != 1){
                throw new Exception("Project accession should be provided for command " +
                        "get-result-files with parameter --app.file-output");
            }
            String fileOutput = fileOutputOptions.get(0);
            prideArchiveWebService.writeResultForProjectAccession(projectAccession, fileOutput);
        }

        // Get the related files to result files
        else if(Objects.equals(command, "get-related-files")){
            List<String> projectAccessionOptions = args.getOptionValues("app.project-accession");
            if(projectAccessionOptions.size() != 1){
                throw new Exception("Project accession should be provided for command " +
                        "get-related-files with parameter --app.project-accession");
            }
            String projectAccession = projectAccessionOptions.get(0);

            List<String> fileOutputOptions = args.getOptionValues("app.file-output");
            if(fileOutputOptions.size() != 1){
                throw new Exception("Project accession should be provided for command " +
                        "get-related-files with parameter --app.file-output");
            }
            String fileOutput = cleanFileName(fileOutputOptions.get(0));

            List<String> resultFileOptions = args.getOptionValues("app.result-file");
            if(resultFileOptions.size() == 0){
                throw new Exception("Project accession should be provided for command " +
                        "get-related-files with parameter --app.result-file");
            }
            resultFileOptions = resultFileOptions.stream().map(this::cleanFileName).collect(Collectors.toList());
            analysisAssayService.writeRelatedSpectraFiles(projectAccession, resultFileOptions, fileOutput);
        }

        // Generate indexes for experiments
        else if(Objects.equals(command, "generate-index-files")){

            List<String> minPSMsOption = args.getOptionValues("app.minPSMs");
            try{
                if(minPSMsOption != null && minPSMsOption.size() > 0) {
                    int minPSMs = Integer.parseInt(minPSMsOption.get(0));
                    analysisAssayService.setMinPSMs(minPSMs);
                }
            }catch (NumberFormatException e){
                throw new Exception("The minimum number of PSMs parameter must be an integer number, example: --app.minPMS=300");
            }

            List<String> valueOption = args.getOptionValues("app.qValueThreshold");
            try{
                if(valueOption != null && valueOption.size() > 0) {
                    double qValueThreshold = Double.parseDouble(valueOption.get(0));
                    analysisAssayService.setQValueThreshold(qValueThreshold);
                }
            }catch (NumberFormatException e){
                throw new Exception("The PSM FDR q-value Threshold (default: 0.01): --app.qValueThreshold=0.05");
            }

            valueOption = args.getOptionValues("app.qFilterProteinFDR");
            try{
                if(valueOption != null && valueOption.size() > 0) {
                    double qValueThreshold = Double.parseDouble(valueOption.get(0));
                    analysisAssayService.setqFilterProteinFDR(qValueThreshold);
                }
            }catch (NumberFormatException e){
                throw new Exception("The Protein FDR q-value Threshold (default: 0.01): --app.qFilterProteinFDR=0.05");
            }

            valueOption = args.getOptionValues("app.peptideLength");
            try{
                if(valueOption != null && valueOption.size() > 0) {
                    double qValueThreshold = Double.parseDouble(valueOption.get(0));
                    analysisAssayService.setPeptideLength(qValueThreshold);
                }
            }catch (NumberFormatException e){
                throw new Exception("The peptide length threshold (default: 7): --app.peptideLength=7");
            }

            valueOption = args.getOptionValues("app.uniquePeptides");
            try{
                if(valueOption != null && valueOption.size() > 0) {
                    double qValueThreshold = Double.parseDouble(valueOption.get(0));
                    analysisAssayService.setPeptideLength(qValueThreshold);
                }
            }catch (NumberFormatException e){
                throw new Exception("The number of unique peptides per protein: --app.uniquePeptides=0");
            }

            List<String> projectAccessionOptions = args.getOptionValues("app.project-accession");
            if(projectAccessionOptions.size() != 1){
                throw new Exception("Project accession should be provided for command " +
                        "generate-index-files with parameter --app.project-accession");
            }
            String projectAccession = projectAccessionOptions.get(0);

            List<String> fileOutputOptions = args.getOptionValues("app.folder-output");
            if(fileOutputOptions.size() != 1){
                throw new Exception("Project accession should be provided for command " +
                        "generate-index-files with parameter --app.folder-output");
            }
            String fileOutput = fileOutputOptions.get(0);

            List<String> resultFileOptions = args.getOptionValues("app.result-file");
            if(resultFileOptions.size() == 0){
                throw new Exception("Result files must be provided for " +
                        "generate-index-files with parameter --app.result-file");
            }
            resultFileOptions = resultFileOptions.stream()
                    .map(this::cleanFileName)
                    .collect(Collectors.toList());

            List<String> spectraFiles = args.getOptionValues("app.spectra-files");
            if(spectraFiles == null)
                spectraFiles = new ArrayList<>();
            else
                spectraFiles = spectraFiles
                        .stream()
                        .flatMap( x-> Arrays.stream(x.split(",")))
                        .map(this::cleanFileName)
                        .filter(x -> (new File(x)).exists())
                        .collect(Collectors.toList());

            List<String> sampleFileOptions = args.getOptionValues("app.sample-file");
            if(sampleFileOptions == null){
                sampleFileOptions = new ArrayList<>();
            }
            List<String> reanalysisAccessionOptions = args.getOptionValues("app.reanalysis-accession");
            if(reanalysisAccessionOptions != null && reanalysisAccessionOptions.size() != 1){
                throw new Exception("Project reanalysis accession should be provided for command " +
                        "generate-index-files with parameter --app.reanalysis-accession");
            }
            String reanalysisAccession = (reanalysisAccessionOptions == null)?null:reanalysisAccessionOptions.get(0);
            try{
                analysisAssayService.writeAnalysisOutputFromResultFiles(projectAccession, resultFileOptions, new HashSet<>(spectraFiles), new HashSet<>(sampleFileOptions), fileOutput, reanalysisAccession);
            }catch (IOException e){
                log.error("Project --- " + projectAccession + "can't be analyzed due the following error --- " + e.getMessage());
            }
        }

        // Perform the protein inference
        else if(Objects.equals(command, "perform-inference")){

            List<String> valueOption = args.getOptionValues("app.qValueThreshold");
            try{
                if(valueOption != null && valueOption.size() > 0) {
                    double qValueThreshold = Double.parseDouble(valueOption.get(0));
                    inferenceAnalysisService.setqValueThreshold(qValueThreshold);
                }
            }catch (NumberFormatException e){
                throw new Exception("The PSM FDR q-value Threshold (default: 0.01): --app.qValueThreshold=0.05");
            }

            valueOption = args.getOptionValues("app.qFilterProteinFDR");
            try{
                if(valueOption != null && valueOption.size() > 0) {
                    double qValueThreshold = Double.parseDouble(valueOption.get(0));
                    inferenceAnalysisService.setqFilterProteinFDR(qValueThreshold);
                }
            }catch (NumberFormatException e){
                throw new Exception("The Protein FDR q-value Threshold (default: 0.01): --app.qFilterProteinFDR=0.05");
            }

            List<String> resultFileOptions = args.getOptionValues("app.archive-spectra");
            if(resultFileOptions.size() != 1){
                throw new Exception("The archive spectra file must be provided --app.archive-spectra");
            }

            List<String> maraClusterOption = args.getOptionValues("app.cluster-file");
            if(resultFileOptions.size() != 1){
                throw new Exception("A file containing the clusters of the spectra --app.cluster-file");
            }

            List<String> outputFolderOption = args.getOptionValues("app.output-folder");
            if(outputFolderOption == null || outputFolderOption.size() > 1)
                throw new Exception("Output folder must be specified --app.output-folder");

            List<String> projectAccessionOption = args.getOptionValues("app.project-accession");
            if(projectAccessionOption.size() != 1){
                throw new Exception("Project Accession must be only one project --app.project-accession");
            }
            String projectAccession = projectAccessionOption.get(0);

            List<String> reanalysisAccessionOption = args.getOptionValues("app.reanalysis-accession");
            String reanalysisAccession = null;
            if(reanalysisAccessionOption != null && reanalysisAccessionOption.size() != 0){
                reanalysisAccession = reanalysisAccessionOption.get(0);
            }
            inferenceAnalysisService.performProteinInference(resultFileOptions.get(0), maraClusterOption.get(0), projectAccession,
                    reanalysisAccession, outputFolderOption.get(0));

        }
        // Convert pride json files to mgf
        else if(Objects.equals(command, "generate-mgf-files")){

            List<String> resultFileOptions = args.getOptionValues("app.archive-spectra");
            if(resultFileOptions.size() != 1){
                throw new Exception("The archive spectra file must be provided --app.archive-spectra");
            }

            List<String> outputFileOptions = args.getOptionValues("app.mgf-file");
            if(resultFileOptions.size() != 1){
                throw new Exception("The mgf file containing all the spectra --app.mgf-file");
            }

            clusteringService.convertToMgf(resultFileOptions.get(0), outputFileOptions.get(0));
        }

        args.getOptionNames().forEach(optionName -> System.out.println(optionName + "=" + args.getOptionValues(optionName)));

    }

    private String cleanFileName(String fileName){
        if (fileName.startsWith("\""))
            fileName = fileName.substring(1);
        if (fileName.endsWith("\""))
            fileName = fileName.substring(0, fileName.length() -1);
        fileName = fileName.replace("\\","");
        fileName = fileName.trim();
        return fileName;
    }
}
