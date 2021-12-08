package uk.ac.ebi.pride.archive.indexer;

import com.amazonaws.auth.AWSStaticCredentialsProvider;
import com.amazonaws.auth.BasicAWSCredentials;
import com.amazonaws.client.builder.AwsClientBuilder;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.autoconfigure.data.mongo.MongoDataAutoConfiguration;
import org.springframework.boot.autoconfigure.mongo.MongoAutoConfiguration;
import org.springframework.context.ConfigurableApplicationContext;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.FilterType;
import uk.ac.ebi.pride.archive.indexer.services.PrideAnalysisAssayService;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideArchiveWebService;
import uk.ac.ebi.pride.archive.spectra.configs.AWS3Configuration;
import uk.ac.ebi.pride.archive.spectra.services.S3SpectralArchive;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

@Slf4j
@SpringBootApplication(exclude={MongoAutoConfiguration.class, MongoDataAutoConfiguration.class})
@ComponentScan(basePackages = "uk.ac.ebi.pride.archive.indexer.services", excludeFilters = @ComponentScan.Filter(type = FilterType.ASSIGNABLE_TYPE, value = {
        S3SpectralArchive.class}))
public class ArchiveMoleculesIndexer implements ApplicationRunner {

    @Value("${batchCommit}")
    private int batchComit = 1;

    private String[] options = {"get-result-files", "get-related-files", "generate-index-files"};

    @Autowired
    private ConfigurableApplicationContext context;

    @Autowired
    private PrideArchiveWebService prideArchiveWebService;

    @Autowired
    PrideAnalysisAssayService analysisAssayService;

    /**
     * usage outputDirectory filterFile directoryToProcess
     *
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
            String fileOutput = fileOutputOptions.get(0);

            List<String> resultFileOptions = args.getOptionValues("app.result-file");
            if(resultFileOptions.size() == 0){
                throw new Exception("Project accession should be provided for command " +
                        "get-related-files with parameter --app.result-file");
            }
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
                throw new Exception("Project accession should be provided for command " +
                        "generate-index-files with parameter --app.result-file");
            }

            List<String> spectraFiles = args.getOptionValues("app.spectra-files");
            if(spectraFiles == null)
                spectraFiles = new ArrayList<>();
            analysisAssayService.writeAnalysisOutputFromResultFiles(projectAccession, resultFileOptions, spectraFiles, fileOutput);

        }

        args.getOptionNames().forEach(optionName -> {
            System.out.println(optionName + "=" + args.getOptionValues(optionName));
        });

    }

    @Bean
    @Qualifier("amazonS3")
    AmazonS3 amazonS3(){
        return AmazonS3ClientBuilder.standard().build();
    }
}
