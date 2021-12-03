package uk.ac.ebi.pride.archive.pipeline;

import com.compomics.util.pride.PrideWebService;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.ApplicationRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.ConfigurableApplicationContext;
import org.springframework.context.annotation.ComponentScan;
import uk.ac.ebi.pride.archive.pipeline.services.ws.PrideArchiveWebService;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

@Slf4j
@SpringBootApplication
@ComponentScan("uk.ac.ebi.pride.archive.pipeline.services")
public class ArchiveMoleculesIndexer implements ApplicationRunner {

    @Value("${batchCommit}")
    private int batchComit = 1;

    private String[] options = {"get-result-files", ""};

    @Autowired
    private ConfigurableApplicationContext context;

    @Autowired
    private PrideArchiveWebService prideArchiveWebService;

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
                throw new Exception("Project accession should be provided for command get-result-files with parameter --app.project-accession");
            }
            String projectAccession = projectAccessionOptions.get(0);

            List<String> fileOutputOptions = args.getOptionValues("app.file-output");
            if(fileOutputOptions.size() != 1){
                throw new Exception("Project accession should be provided for command get-result-files with parameter --app.file-output");
            }
            String fileOutput = fileOutputOptions.get(0);
            prideArchiveWebService.writeResultForProjectAccession(projectAccession, fileOutput);
        }

        args.getOptionNames().forEach(optionName -> {
            System.out.println(optionName + "=" + args.getOptionValues(optionName));
        });
        SpringApplication.exit(context);

    }

}
