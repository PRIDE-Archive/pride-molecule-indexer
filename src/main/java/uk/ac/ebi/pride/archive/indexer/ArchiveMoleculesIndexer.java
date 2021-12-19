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
import uk.ac.ebi.pride.archive.indexer.services.PrideAnalysisAssayService;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideArchiveWebService;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

@Slf4j
@SpringBootApplication
@ComponentScan(basePackages = "uk.ac.ebi.pride.archive.indexer.services")
public class ArchiveMoleculesIndexer implements ApplicationRunner {

    @Value("${batchCommit}")
    private int batchComit = 1;

    private final String[] options = {"get-result-files", "get-related-files", "generate-index-files"};

    @Autowired
    private ConfigurableApplicationContext context;

    @Autowired
    private PrideArchiveWebService prideArchiveWebService;

    @Autowired
    PrideAnalysisAssayService analysisAssayService;

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
        else if(Objects.equals(command, "generate-index-files")){
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

            analysisAssayService.writeAnalysisOutputFromResultFiles(projectAccession, resultFileOptions, new HashSet<>(spectraFiles), new HashSet<>(sampleFileOptions), fileOutput);
        }

        args.getOptionNames().forEach(optionName -> System.out.println(optionName + "=" + args.getOptionValues(optionName)));

    }

    private String cleanFileName(String fileName){
        if (fileName.startsWith("'"))
            fileName = fileName.substring(1);
        if (fileName.endsWith("'"))
            fileName = fileName.substring(0, fileName.length() -1);
        return fileName;
    }
}
