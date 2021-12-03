package uk.ac.ebi.pride.archive.pipeline.services;


import de.mpc.pia.intermediate.Accession;
import de.mpc.pia.intermediate.AccessionOccurrence;
import de.mpc.pia.intermediate.Modification;
import de.mpc.pia.intermediate.PeptideSpectrumMatch;
import de.mpc.pia.modeller.PIAModeller;
import de.mpc.pia.modeller.peptide.ReportPeptide;
import de.mpc.pia.modeller.protein.ReportProtein;
import de.mpc.pia.modeller.psm.ReportPSM;
import de.mpc.pia.modeller.report.filter.AbstractFilter;
import de.mpc.pia.modeller.report.filter.FilterComparator;
import de.mpc.pia.modeller.report.filter.RegisteredFilters;
import de.mpc.pia.modeller.report.filter.impl.PSMScoreFilter;
import de.mpc.pia.modeller.score.ScoreModelEnum;
import de.mpc.pia.tools.pride.PRIDETools;
import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.dao.DuplicateKeyException;
import org.springframework.stereotype.Service;
import uk.ac.ebi.jmzidml.model.mzidml.AbstractParam;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.pride.archive.dataprovider.assay.AssayType;
import uk.ac.ebi.pride.archive.dataprovider.common.Tuple;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PSMProvider;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PeptideSpectrumOverview;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModification;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModificationProvider;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParam;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParamProvider;
import uk.ac.ebi.pride.archive.pipeline.services.proteomics.JmzReaderSpectrumService;
import uk.ac.ebi.pride.archive.pipeline.services.proteomics.PIAModelerService;
import uk.ac.ebi.pride.archive.pipeline.services.ws.PrideArchiveWebService;
import uk.ac.ebi.pride.archive.pipeline.services.ws.PrideFile;
import uk.ac.ebi.pride.archive.pipeline.services.ws.PrideProject;
import uk.ac.ebi.pride.archive.pipeline.utility.BackupUtil;
import uk.ac.ebi.pride.archive.pipeline.utility.SubmissionPipelineConstants;
import uk.ac.ebi.pride.archive.spectra.model.ArchiveSpectrum;
import uk.ac.ebi.pride.mongodb.archive.model.assay.MongoAssayFile;
import uk.ac.ebi.pride.mongodb.archive.model.assay.MongoPrideAssay;
import uk.ac.ebi.pride.mongodb.molecules.model.peptide.PrideMongoPeptideEvidence;
import uk.ac.ebi.pride.mongodb.molecules.model.protein.PrideMongoProteinEvidence;
import uk.ac.ebi.pride.mongodb.molecules.model.psm.PrideMongoPsmSummaryEvidence;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.utilities.term.CvTermReference;
import uk.ac.ebi.pride.utilities.util.MoleculeUtilities;
import uk.ac.ebi.pride.utilities.util.Triple;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.AccessDeniedException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;
import java.util.stream.Collectors;

@Configuration
@Slf4j
@Service
public class PrideAnalysisAssayService {

    private static final Long MERGE_FILE_ID = 1L;
    private static final Long FILE_ID = 1L;
    private static final int PIPELINE_RETRY_LIMIT = 20;

    @Autowired
    PrideArchiveWebService prideFileMongoRepository;

    private PIAModelerService piaModellerService;

    Map<String, Long> taskTimeMap = new HashMap<>();

    @Value("${productionPath}")
    String productionPath;

    String projectPath;

    @Value("${backupPath}")
    String backupPath;

    PrideProject project;

    //@Value("#{jobParameters['project']}")
    private String projectAccession;

    private String resultFileName;

    @Value("${qValueThreshold:#{0.01}}")
    private Double qValueThreshold;

    @Value("${qFilterProteinFDR:#{1.0}}")
    private Double qFilterProteinFDR;

    List<PrideFile> mongoResultFiles;

    Map<String, Map<String, Object>> assayObjectMap;

    List<MongoPrideAssay> assayList = new ArrayList<>();

    String buildPath;

    DecimalFormat df = new DecimalFormat("###.#####");
    private List<PrideFile> mongoPrideFiles;

    @Bean
    PIAModelerService getPIAModellerService() {
        piaModellerService = new PIAModelerService();
        return piaModellerService;
    }

    public void initJobPRIDEReanalysisAnalyzeAssayJob() throws IOException {

        Optional<PrideProject> project = prideFileMongoRepository.findByAccession(projectAccession);
        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM");
        String allDate = dateFormat.format(project.get().getPublicationDate());
        String[] allDateString = allDate.split("-");
        String year = null, month = null;

        if (allDateString.length == 2) {
            year = allDateString[0];
            month = allDateString[1];
        }
        if (year != null && month != null) {
            projectPath = String.join("/",productionPath, year, month, projectAccession);
        }
        mongoPrideFiles = prideFileMongoRepository.findFilesByProjectAccession(projectAccession);
        this.mongoResultFiles = mongoPrideFiles.parallelStream().filter(mongoFile -> mongoFile.getFileCategory().getValue().equals("RESULT")).collect(Collectors.toList());
        if(resultFileName != null){
            this.mongoResultFiles = mongoResultFiles.stream().filter( x-> Objects.equals(x.getFileName(), resultFileName)).collect(Collectors.toList());
        }
        this.assayObjectMap = new HashMap<>();
        for (PrideFile mongoResultFile : mongoResultFiles) {
            log.info("Analyzing assay file  -- " + mongoResultFile.getFileName());
            Map<String, Object> assayObjects = new HashMap<>();
            assayObjectMap.put(mongoResultFile.getAccession(), assayObjects);
        }
        System.out.println(String.format("==================>>>>>>> PRIDEReanalysisAnalyzeAssayJob - Run the job for Project %s", projectAccession));
    }

    private void createBackupFiles() throws IOException {
        createBackupDir();
        for (PrideFile mongoResultFile : mongoResultFiles) {

            log.info("Analyzing assay file  -- " + mongoResultFile.getFileName());

            Map<String, Object> assayObjects = assayObjectMap.get(mongoResultFile.getAccession());
            final String peptideEvidenceFileName = BackupUtil.getPrideMongoPeptideEvidenceFile(backupPath, projectAccession, mongoResultFile.getAccession());

            assayObjects.put("peptideEvidenceBufferedWriter", new BufferedWriter(new FileWriter(peptideEvidenceFileName, false)));

            final String proteinEvidenceFileName = BackupUtil.getPrideMongoProteinEvidenceFile(backupPath, projectAccession, mongoResultFile.getAccession());
            assayObjects.put("proteinEvidenceBufferedWriter", new BufferedWriter(new FileWriter(proteinEvidenceFileName, false)));

            final String archiveSpectrumFileName = BackupUtil.getArchiveSpectrumFile(backupPath, projectAccession, mongoResultFile.getAccession());
            assayObjects.put("archiveSpectrumBufferedWriter", new BufferedWriter(new FileWriter(archiveSpectrumFileName, false)));

            final String psmSummaryEvidenceFileName = BackupUtil.getPrideMongoPsmSummaryEvidenceFile(backupPath, projectAccession, mongoResultFile.getAccession());
            assayObjects.put("psmSummaryEvidenceBufferedWriter", new BufferedWriter(new FileWriter(psmSummaryEvidenceFileName, false)));

            assayObjectMap.put(mongoResultFile.getAccession(), assayObjects);
        }
    }

    private void createBackupDir() throws AccessDeniedException {
        String path = backupPath;
        if (!path.endsWith(File.separator)) {
            path = backupPath + File.separator;
        }
        path = path + projectAccession;
        File file = new File(path);
        if (file.exists() && file.isDirectory()) {
            return;
        }
        boolean mkdirs = file.mkdirs();
        if (!mkdirs) {
            throw new AccessDeniedException("Failed to create Dir : " + backupPath);
        }
    }

    public void analyzeAssayInformationStep() throws IOException {

        long initAnalysisAssay = System.currentTimeMillis();

        log.info("Analyzing project -- " + projectAccession);
        log.info("creating backup files");
        createBackupFiles();
        log.info("creating backup files: finished");

        Optional<PrideProject> optionalProject = prideFileMongoRepository.findByAccession(projectAccession);

        if (mongoResultFiles.size() < 0 || !optionalProject.isPresent()) {
            String errorMessage = "No Project or Assay found!";
            log.error(errorMessage);
            throw new IOException(errorMessage);
        }

        project = optionalProject.get();

        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM");
        String allDate = dateFormat.format(project.getPublicationDate());
        String[] allDateString = allDate.split("-");
        String year = null, month = null;

        if (allDateString.length == 2) {
            year = allDateString[0];
            month = allDateString[1];
        }

        if (year != null && month != null) {
            //change it local path to test
            String buildPath = SubmissionPipelineConstants.buildInternalPath(productionPath,projectAccession, year, month);
            for (PrideFile mongoResultFile : mongoResultFiles) {

                log.info("Analyzing assay file  -- " + mongoResultFile.getFileName());

                long nrDecoys = 0;

                Map<String, Object> assayObjects = assayObjectMap.get(mongoResultFile.getAccession());
                SubmissionPipelineConstants.FileType fileType = SubmissionPipelineConstants.FileType.getFileTypeFromPRIDEFileName(mongoResultFile.getFileName());

                /**
                 * The first threshold for modeller is not threshold at PSM and Protein level.
                 */
                PIAModeller modeller = piaModellerService.performProteinInference(mongoResultFile.getAccession(),
                        SubmissionPipelineConstants.returnUnCompressPath(buildPath + mongoResultFile.getFileName()),
                        fileType, 1.0, 1.0);

                nrDecoys = modeller.getPSMModeller().getReportPSMSets().entrySet().stream()
                        .filter(entry -> entry.getValue().getIsDecoy())
                        .count();

                List<AbstractFilter> filters = new ArrayList<>();
                filters.add(RegisteredFilters.PSM_SOURCE_ID_FILTER
                        .newInstanceOf(FilterComparator.equal, "index=null", true));

                assayObjects.put("allPsms", modeller.getPSMModeller()
                        .getFilteredReportPSMs(FILE_ID, filters));
                assayObjects.put("allPeptides", modeller.getPeptideModeller()
                        .getFilteredReportPeptides(MERGE_FILE_ID, filters));
                assayObjects.put("allProteins", modeller.getProteinModeller()
                        .getFilteredReportProteins(filters));


                // setting filter for peptide level filtering
                modeller = piaModellerService.performFilteringInference(modeller, qValueThreshold, qFilterProteinFDR);
                filters.add(new PSMScoreFilter(FilterComparator.less_equal, false,
                        qValueThreshold, ScoreModelEnum.PSM_LEVEL_Q_VALUE.getShortName()));              // you can also use fdr score here

                // get the FDR filtered highQualityPeptides
                List<ReportPSM> highQualityPsms = modeller.getPSMModeller()
                        .getFilteredReportPSMs(MERGE_FILE_ID, filters);
                List<ReportPeptide> highQualityPeptides = modeller.getPeptideModeller()
                        .getFilteredReportPeptides(MERGE_FILE_ID, filters);
                List<ReportProtein> highQualityProteins = modeller.getProteinModeller()
                        .getFilteredReportProteins(filters);


                if (!(nrDecoys > 0 && highQualityProteins.size() > 0 && highQualityPeptides.size() > 0 && highQualityPsms.size() > 0 &&
                        highQualityPsms.size() >= highQualityPeptides.size())) {
                    highQualityPeptides = new ArrayList<>();
                    highQualityProteins = new ArrayList<>();
                    highQualityPsms = new ArrayList<>();
                }

                assayObjects.put("modeller", modeller);
                assayObjects.put("highQualityPsms", highQualityPsms);
                assayObjects.put("highQualityPeptides", highQualityPeptides);
                assayObjects.put("highQualityProteins", highQualityProteins);
            }
        } else {
            String errorMessage = "The Year and Month for Project Accession can't be found -- " + project.getAccession();
            log.error(errorMessage);
            throw new IOException(errorMessage);
        }

        taskTimeMap.put(SubmissionPipelineConstants.PrideArchiveStepNames.PRIDE_ARCHIVE_MONGODB_ASSAY_INFERENCE.getName(),
                System.currentTimeMillis() - initAnalysisAssay);

    }


    void updateAssayInformationStep() {
        for (PrideFile mongoResultFile : mongoResultFiles) {
            Map<String, Object> assayObject = assayObjectMap.get(mongoResultFile.getAccession());
            if (mongoResultFiles != null && assayObject.get("modeller") != null) {

                PIAModeller modeller = (PIAModeller) assayObject.get("modeller");

                List<ReportPeptide> highQualityPeptides = (List<ReportPeptide>) assayObject.get("highQualityPeptides");
                List<ReportProtein> highQualityProteins = (List<ReportProtein>) assayObject.get("highQualityProteins");
                List<ReportPSM> highQualityPsms = (List<ReportPSM>) assayObject.get("highQualityPsms");

                List<ReportPeptide> modifiedPeptides = highQualityPeptides.
                        stream().filter(x -> x.getModifications().size() > 0)
                        .collect(Collectors.toList());

                Set<CvParam> summaryResults = new HashSet<>();


                CvParam peptideParam = new CvParam(CvTermReference.PRIDE_NUMBER_ID_PEPTIDES.getCvLabel(),
                        CvTermReference.PRIDE_NUMBER_ID_PEPTIDES.getAccession(),
                        CvTermReference.PRIDE_NUMBER_ID_PEPTIDES.getName(),
                        String.valueOf(highQualityPeptides.size()));
                summaryResults.add(peptideParam);

                CvParam proteinParam = new CvParam(CvTermReference.PRIDE_NUMBER_ID_PROTEINS.getCvLabel(),
                        CvTermReference.PRIDE_NUMBER_ID_PROTEINS.getAccession(),
                        CvTermReference.PRIDE_NUMBER_ID_PROTEINS.getName(),
                        String.valueOf(highQualityProteins.size()));
                summaryResults.add(proteinParam);

                CvParam psmParam = new CvParam(CvTermReference.PRIDE_NUMBER_ID_PSMS.getCvLabel(),
                        CvTermReference.PRIDE_NUMBER_ID_PSMS.getAccession(),
                        CvTermReference.PRIDE_NUMBER_ID_PSMS.getName(),
                        String.valueOf(highQualityPsms.size()));
                summaryResults.add(psmParam);

                CvParam modifiedPeptideParam = new CvParam(CvTermReference.PRIDE_NUMBER_ID_PEPTIDES.getCvLabel(),
                        CvTermReference.PRIDE_NUMBER_ID_PEPTIDES.getAccession(),
                        CvTermReference.PRIDE_NUMBER_ID_PEPTIDES.getName(),
                        String.valueOf(modifiedPeptides.size()));
                summaryResults.add(modifiedPeptideParam);


                List<Tuple<CvParam, Integer>> modificationCount = modifiedPeptides.stream()
                        .flatMap(x -> x.getModifications().values().stream())
                        .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                        .entrySet().stream()
                        .map(entry -> new Tuple<>(new CvParam(entry.getKey().getCvLabel(), entry.getKey().getAccession(), entry.getKey().getDescription(), String.valueOf(entry.getKey().getMass())), entry.getValue().intValue()))
                        .collect(Collectors.toList());

                boolean isValid;


                if (highQualityPeptides.size() > 0 && highQualityProteins.size() > 0 && highQualityPsms.size() > 0)
                    isValid = true;
                else
                    isValid = false;

                Set<CvParam> validationMethods = new HashSet<>();

                if (isValid) {
                    validationMethods.add(new CvParam(CvTermReference.MS_DECOY_VALIDATION_METHOD.getCvLabel(),
                            CvTermReference.MS_DECOY_VALIDATION_METHOD.getAccession(), CvTermReference.MS_DECOY_VALIDATION_METHOD.getName(), String.valueOf(true)));
                } else
                    validationMethods.add(new CvParam(CvTermReference.MS_DECOY_VALIDATION_METHOD.getCvLabel(),
                            CvTermReference.MS_DECOY_VALIDATION_METHOD.getAccession(), CvTermReference.MS_DECOY_VALIDATION_METHOD.getName(), String.valueOf(false)));

                assayObject.put("validationMethods", validationMethods);
                assayObject.put("isValid", isValid);

                List<MongoAssayFile> mongoAssayFiles = new ArrayList<>();

                //TO-DO relatedFiles

                mongoAssayFiles.add(MongoAssayFile.builder().fileName(mongoResultFile.getFileName())
                        .fileAccession(mongoResultFile.getAccession())
                        .fileCategory((CvParam) mongoResultFile.getFileCategory())
                        // .relatedFiles(mo)
                        .build());


                MongoPrideAssay mongoPrideAssay = MongoPrideAssay.builder()
                        .assayType(AssayType.IDENTIFICATION)
                        .accession(mongoResultFile.getAccession())
                        .title(project.getTitle())
                        .projectAccessions(Collections.singleton(project.getAccession()))
                        .analysisAccessions(Collections.singleton(project.getAccession()))
                        .dataAnalysisSoftwares(convertMapOfAnalysisSoftwareToSet(modeller.getAnalysisSoftwares()))
                        .summaryResults(summaryResults)
                        .ptmsResults(getPtmsResults(modifiedPeptides))
                        .assayFiles(mongoAssayFiles)
                        .build();
                assayList.add(mongoPrideAssay);
            }


        }

//                    prideProjectMongoService.saveAssays(assayList);



    }

    private List<Tuple<CvParam, Integer>> getPtmsResults(List<ReportPeptide> modifiedPeptides) {
        List<Tuple<CvParam, Integer>> ptmsResultSet = new ArrayList<>();
        Set<CvParam> CvParamSet = new HashSet<>();
        modifiedPeptides.stream().forEach(pep -> {
            Collection<? extends IdentifiedModificationProvider> imp = convertPeptideModifications(pep.getModifications());
            imp.forEach(
                    ptm -> {
                        CvParam cvParam = (CvParam) ptm.getModificationCvTerm();
                        if (!CvParamSet.contains(cvParam)) {
                            ptmsResultSet.add(new Tuple(cvParam, 0));
                            CvParamSet.add(cvParam);
                        }
                    }
            );
        });
        return ptmsResultSet;
    }

    private Set<CvParam> convertMapOfAnalysisSoftwareToSet(Map<String, AnalysisSoftware> analysisSoftwareMap) {
        Set<CvParam> analysisSoftwareSet = new HashSet<>();
        analysisSoftwareMap.entrySet().stream().forEach(entry -> {
            AnalysisSoftware analysisSoftware = entry.getValue();
            analysisSoftwareSet.add(new CvParam(entry.getKey(), analysisSoftware.getId(), analysisSoftware.getSoftwareName().getCvParam().getValue()
                    , analysisSoftware.getName()));
        });
        return analysisSoftwareSet;
    }

    /**
     * Defines the job to Sync all the projects from OracleDB into MongoDB database.
     *
     * @return the calculatePrideArchiveDataUsage job
     */
    public void reanalysisAnalyzeAssayInformationJob() throws Exception {
        initJobPRIDEReanalysisAnalyzeAssayJob();
        analyzeAssayInformationStep();
        updateAssayInformationStep();
        indexSpectraStep();
        proteinPeptideIndexStep();
        analyzeAssayPrintTraceStep();
    }

    public void analyzeAssayPrintTraceStep() throws IOException {
        taskTimeMap.forEach((key, value) -> log.info("Task: " + key + " Time: " + value));

        Set<Map.Entry<String, Map<String, Object>>> assaySet = assayObjectMap.entrySet();
        for (Map.Entry<String, Map<String, Object>> assay : assaySet) {
            Map<String, Object> assayObject = assay.getValue();
            ((BufferedWriter) assayObject.get("proteinEvidenceBufferedWriter")).close();
            ((BufferedWriter) assayObject.get("peptideEvidenceBufferedWriter")).close();
            ((BufferedWriter) assayObject.get("archiveSpectrumBufferedWriter")).close();
            ((BufferedWriter) assayObject.get("psmSummaryEvidenceBufferedWriter")).close();
        }
    }

    public void proteinPeptideIndexStep() throws Exception {
        long initInsertPeptides = System.currentTimeMillis();

        for (PrideFile mongoResultFile : mongoResultFiles) {
            Map<String, Object> assayObjects = assayObjectMap.get(mongoResultFile.getAccession());
            if (mongoResultFiles != null && assayObjectMap.get("modeller") != null) {

                log.info("proteinPeptideIndexStep assay file  -- " + mongoResultFile.getFileName());

                PIAModeller modeller = (PIAModeller) assayObjects.get("modeller");

                List<ReportPeptide> highQualityPeptides = (List<ReportPeptide>) assayObjects.get("highQualityPeptides");
                List<ReportProtein> highQualityProteins = (List<ReportProtein>) assayObjects.get("highQualityProteins");
                List<ReportPSM> highQualityPsms = (List<ReportPSM>) assayObjects.get("highQualityPsms");
                List<ReportPSM> allPsms = (List<ReportPSM>) assayObjects.get("allPsms");
                List<ReportPeptide> allPeptides = (List<ReportPeptide>) assayObjects.get("allPeptides");
                List<ReportProtein> allProteins = (List<ReportProtein>) assayObjects.get("allProteins");

                Set<String> proteinIds = new HashSet<>();
                Set<String> peptideSequences = new HashSet<>();

                List<ReportPeptide> peptides;
                List<ReportProtein> proteins;
                List<ReportPSM> psms = new ArrayList<>();

                if (highQualityPeptides.size() > 0 && highQualityPsms.size() > 0 && highQualityProteins.size() > 0) {
                    peptides = highQualityPeptides;
                    proteins = highQualityProteins;
                    psms = highQualityPsms;
                } else {
                    peptides = allPeptides;
                    proteins = allProteins;
                    psms = allPsms;
                }

                List<ReportPeptide> finalPeptides = peptides;

                for (ReportProtein protein : proteins) {
                    String proteinSequence = protein.getRepresentative().getDbSequence();
                    String proteinAccession = protein.getRepresentative().getAccession();
                    Set<String> proteinGroups = protein.getAccessions()
                            .stream().map(Accession::getAccession)
                            .collect(Collectors.toSet());

                    List<IdentifiedModificationProvider> proteinPTMs = new ArrayList<>(convertProteinModifications(
                            proteinAccession, protein.getPeptides()));

                    log.info(String.valueOf(protein.getQValue()));

                    CvParam scoreParam = null;
                    Set<CvParam> attributes = new HashSet<>();

                    if (!Double.isFinite(protein.getQValue()) && !Double.isNaN(protein.getQValue())) {

                        String value = df.format(protein.getQValue());

                        scoreParam = new CvParam(CvTermReference.MS_PIA_PROTEIN_GROUP_QVALUE.getCvLabel(),
                                CvTermReference.MS_PIA_PROTEIN_GROUP_QVALUE.getAccession(),
                                CvTermReference.MS_PIA_PROTEIN_GROUP_QVALUE.getName(), value);
                        attributes.add(scoreParam);
                    }

                    if (protein.getScore() != null && !protein.getScore().isNaN()) {
                        String value = df.format(protein.getScore());
                        scoreParam = new CvParam(CvTermReference.MS_PIA_PROTEIN_SCORE.getCvLabel(),
                                CvTermReference.MS_PIA_PROTEIN_SCORE.getAccession(),
                                CvTermReference.MS_PIA_PROTEIN_SCORE.getName(), value);
                        attributes.add(scoreParam);
                    }

                    AtomicReference<CvParam> param = new AtomicReference<>(new CvParam(PRIDETools.PrideOntologyConstants
                            .PRIDE_SUBMITTERS_THERSHOLD.getCvLabel(),
                            PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getAccession(),
                            PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getName(),
                            Boolean.toString(false)));

                    AtomicBoolean submitterValid = new AtomicBoolean(false);
                    protein.getPeptides().stream().forEach(x -> {
                        x.getPeptide().getSpectra().stream().forEach(y -> {
                            for (AbstractParam abstractParam : y.getParams()) {
                                if (abstractParam instanceof uk.ac.ebi.jmzidml.model.mzidml.CvParam) {
                                    uk.ac.ebi.jmzidml.model.mzidml.CvParam cv = (uk.ac.ebi.jmzidml.model.mzidml.CvParam) abstractParam;
                                    if (cv.getAccession().equalsIgnoreCase("PRIDE:0000511") &&
                                            cv.getValue().equalsIgnoreCase("true")) {
                                        param.set(new CvParam(PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD
                                                .getCvLabel(),
                                                PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getAccession(),
                                                PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getName(),
                                                Boolean.toString(true)));
                                        submitterValid.set(true);
                                    }
                                }
                            }
                        });
                    });

                    attributes.add(param.get());
                    proteinIds.add(proteinAccession);
                    protein.getPeptides().forEach(x -> peptideSequences.add(x.getSequence()));

                    Set<CvParam> validationMethods = (Set<CvParam>) assayObjects.get("validationMethods");

                    boolean isValid = (boolean) assayObjects.get("isValid");

                    PrideMongoProteinEvidence proteinEvidence = PrideMongoProteinEvidence
                            .builder()
                            .reportedAccession(proteinAccession)
                            .isDecoy(protein.getIsDecoy())
                            .proteinGroupMembers(proteinGroups)
                            .ptms(proteinPTMs)
                            .projectAccession(projectAccession)
                            .proteinSequence(proteinSequence)
                            .bestSearchEngineScore(scoreParam)
                            .additionalAttributes(attributes)
                            .assayAccession(mongoResultFile.getAccession())
                            .isValid(isValid)
                            .qualityEstimationMethods(validationMethods)
                            .numberPeptides(protein.getPeptides().size())
                            .numberPSMs(protein.getNrPSMs())
                            .sequenceCoverage(protein.getCoverage(proteinAccession))
                            .build();

                    if (isValid || submitterValid.get()) {
                        try {
                            BackupUtil.write(proteinEvidence, (BufferedWriter) assayObjects.get("proteinEvidenceBufferedWriter"));
//                                            moleculesService.insertProteinEvidences(proteinEvidence);
                        } catch (DuplicateKeyException ex) {
//                                            moleculesService.saveProteinEvidences(proteinEvidence);
                            log.debug("The protein was already in the database -- " + proteinEvidence.getReportedAccession());
                        } catch (Exception e) {
                            log.error(e.getMessage(), e);
                            throw new Exception(e);
                        }
                        indexPeptideByProtein(protein, finalPeptides, (Map<Long, List<PeptideSpectrumOverview>>) assayObjects.get("peptideUsi")
                                , (Set<CvParam>) assayObjects.get("validationMethods"),
                                (BufferedWriter) assayObjects.get("peptideEvidenceBufferedWriter"),
                                mongoResultFile.getAccession(), isValid);
                    }
                }
            }
        }
    }

    /**
     * This method index all the highQualityPeptides that identified a protein into the mongoDB
     *
     * @param protein  Identified Protein
     * @param peptides Collection of identified highQualityPeptides in the experiment
     */
    private void indexPeptideByProtein(ReportProtein protein, List<ReportPeptide> peptides, Map<Long, List<PeptideSpectrumOverview>> peptideUsi
            , Set<CvParam> validationMethods, BufferedWriter peptideEvidenceBufferedWriter, String assayAccession, boolean isValid) throws Exception {

        for (ReportPeptide peptide : protein.getPeptides()) {
            Optional<ReportPeptide> firstPeptide = peptides.stream()
                    .filter(globalPeptide -> globalPeptide.getStringID().equalsIgnoreCase(peptide.getStringID()))
                    .findFirst();

            if (firstPeptide.isPresent()) {

                Set<CvParam> peptideAttributes = new HashSet<>();
                if (!Double.isInfinite(firstPeptide.get().getQValue()) && !Double.isNaN(firstPeptide.get().getQValue())) {

                    String value = df.format(firstPeptide.get().getQValue());

                    CvParam peptideScore = new CvParam(CvTermReference.MS_PIA_PEPTIDE_QVALUE
                            .getCvLabel(),
                            CvTermReference.MS_PIA_PEPTIDE_QVALUE.getAccession(),
                            CvTermReference.MS_PIA_PEPTIDE_QVALUE.getName(), value);
                    peptideAttributes.add(peptideScore);
                }


                if (!Double.isInfinite(firstPeptide.get().getScore("peptide_fdr_score"))
                        && !Double.isNaN(firstPeptide.get().getScore("peptide_fdr_score"))) {

                    String value = df.format(firstPeptide.get().getScore("peptide_fdr_score"));

                    CvParam peptideScore = new CvParam(CvTermReference.MS_PIA_PEPTIDE_FDR
                            .getCvLabel(),
                            CvTermReference.MS_PIA_PEPTIDE_FDR.getAccession(),
                            CvTermReference.MS_PIA_PEPTIDE_FDR.getName(), value);
                    peptideAttributes.add(peptideScore);
                }

                List<PeptideSpectrumOverview> usiList = peptideUsi.get(firstPeptide.get().getPeptide().getID());

                int startPosition = 0;
                int endPosition = 0;

                Optional<AccessionOccurrence> occurrence = firstPeptide.get().getPeptide().getAccessionOccurrences().stream()
                        .filter(x -> x.getAccession().getAccession().equalsIgnoreCase(protein.getRepresentative().getAccession()))
                        .findFirst();
                if (occurrence.isPresent()) {
                    startPosition = occurrence.get().getStart();
                    endPosition = occurrence.get().getEnd();
                } else {
                    log.info("Position of the corresponding peptide is not present -- " + protein.getRepresentative().getAccession());
                }

                AtomicReference<CvParam> param = new AtomicReference<>(new CvParam(PRIDETools.PrideOntologyConstants
                        .PRIDE_SUBMITTERS_THERSHOLD.getCvLabel(),
                        PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getAccession(),
                        PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getName(), Boolean.toString(false)));

                AtomicBoolean submitterValid = new AtomicBoolean(false);
                peptides.stream().forEach(x -> {
                    x.getPeptide().getSpectra().stream().forEach(y -> {
                        for (AbstractParam abstractParam : y.getParams()) {
                            if (abstractParam instanceof uk.ac.ebi.jmzidml.model.mzidml.CvParam) {
                                if (abstractParam instanceof uk.ac.ebi.jmzidml.model.mzidml.CvParam) {
                                    uk.ac.ebi.jmzidml.model.mzidml.CvParam cv = (uk.ac.ebi.jmzidml.model.mzidml.CvParam) abstractParam;
                                    if (cv.getAccession().equalsIgnoreCase("PRIDE:0000511") && cv.getValue().equalsIgnoreCase("true")) {
                                        param.set(new CvParam(PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getCvLabel(),
                                                PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getAccession(),
                                                PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getName(),
                                                Boolean.toString(true)));
                                        submitterValid.set(true);
                                    }
                                }
                            }
                        }
                    });
                });

                peptideAttributes.add(param.get());

                int misssedCleavages = firstPeptide.get().getMissedCleavages();
                if (misssedCleavages == -1){
                    misssedCleavages = uk.ac.ebi.pride.utilities.mol.MoleculeUtilities.calcMissedCleavages(peptide.getSequence());
                }

                PrideMongoPeptideEvidence peptideEvidence = PrideMongoPeptideEvidence
                        .builder()
                        .assayAccession(assayAccession)
                        .proteinAccession(protein.getRepresentative().getAccession())
                        .isDecoy(firstPeptide.get().getIsDecoy())
                        .peptideAccession(SubmissionPipelineConstants
                                .encodePeptide(peptide.getSequence(), peptide.getModifications()))
                        .peptideSequence(peptide.getSequence())
                        .additionalAttributes(peptideAttributes)
                        .projectAccession(projectAccession)
                        .psmAccessions(usiList)
                        .startPosition(startPosition)
                        .endPosition(endPosition)
                        .missedCleavages(misssedCleavages)
                        .ptmList(convertPeptideModifications(firstPeptide.get().getModifications()))
                        .isValid(isValid)
                        .qualityEstimationMethods(validationMethods)
                        .build();

                if (isValid || submitterValid.get()) {
                    try {
                        BackupUtil.write(peptideEvidence, peptideEvidenceBufferedWriter);
//                        moleculesService.insertPeptideEvidence(peptideEvidence);
                    } catch (DuplicateKeyException ex) {
//                       moleculesService.savePeptideEvidence(peptideEvidence);
                        log.debug("The peptide evidence was already in the database -- " + peptideEvidence.getPeptideAccession());
                    } catch (Exception e) {
                        log.error(e.getMessage(), e);
                        throw new Exception(e);
                    }
                }
            }
        }

    }

    /**
     * Convert Peptide Modifications from PIA modeller to PeptideEvidence modifications
     *
     * @param modifications Modifications Map
     * @return List if {@link IdentifiedModificationProvider}
     */
    private Collection<? extends IdentifiedModificationProvider> convertPeptideModifications(Map<Integer, Modification> modifications) {

        List<IdentifiedModification> ptms = new ArrayList<>();

        for (Map.Entry<Integer, Modification> ptmEntry : modifications.entrySet()) {
            Modification ptm = ptmEntry.getValue();
            Integer position = ptmEntry.getKey();
            Set<CvParam> probabilities = ptm.getProbability()
                    .stream().map(oldProbability -> new CvParam(oldProbability.getCvLabel(),
                            oldProbability.getAccession(),
                            oldProbability.getName(),
                            String.valueOf(oldProbability.getValue())))
                    .collect(Collectors.toSet());
            // ignore modifications that can't be processed correctly (can not be mapped to the protein)
            if (ptm.getAccession() == null) {
                continue;
            }

            Optional<IdentifiedModification> proteinExist = ptms.stream()
                    .filter(currentMod -> currentMod.getModificationCvTerm()
                            .getAccession().equalsIgnoreCase(ptm.getAccession()))
                    .findAny();
            if (proteinExist.isPresent()) {
                proteinExist.get().addPosition(position, probabilities);
            } else {
                CvParam ptmName = new CvParam(ptm.getCvLabel(),
                        ptm.getAccession(), ptm.getDescription(),
                        String.valueOf(ptm.getMass()));
                IdentifiedModification newPTM = new IdentifiedModification(null, null, ptmName, null);
                newPTM.addPosition(position, probabilities);
                ptms.add(newPTM);
            }
        }
        return ptms;

    }

    /**
     * Convert peptide modifications to Protein modifications. Adjust the localization using the start and end positions.
     *
     * @param proteinAccession Protein Accession
     * @param peptides         List of highQualityPeptides
     * @return List of {@link IdentifiedModificationProvider}
     */
    private Collection<? extends IdentifiedModificationProvider> convertProteinModifications(String proteinAccession, List<ReportPeptide> peptides) {

        List<IdentifiedModification> ptms = new ArrayList<>();

        for (ReportPeptide item : peptides) {

            for (Map.Entry<Integer, Modification> ptmEntry : item.getModifications().entrySet()) {

                Modification ptm = ptmEntry.getValue();
                Integer position = ptmEntry.getKey();
                Set<CvParam> probabilities = ptm.getProbability()
                        .stream().map(oldProbability -> new CvParam(oldProbability.getCvLabel(),
                                oldProbability.getAccession(),
                                oldProbability.getName(),
                                String.valueOf(oldProbability.getValue())))
                        .collect(Collectors.toSet());
                // ignore modifications that can't be processed correctly (can not be mapped to the protein)
                if (ptm.getAccession() == null) {
                    continue;
                }

                // if we can calculate the position, we add it to the modification
                // -1 to calculate properly the modification offset
                item.getPeptide().getAccessionOccurrences().forEach(peptideEvidence -> {

                    if (peptideEvidence.getAccession().getAccession() == proteinAccession) {

                        if (peptideEvidence.getStart() != null && peptideEvidence.getStart() >= 0 && position >= 0) {

                            int startPos = peptideEvidence.getStart();
                            // n-term and c-term mods are not propagated to the protein except the case that the start
                            // position is 1 (beginning of the protein)
                            int proteinPosition = startPos + position - 1;

                            Optional<IdentifiedModification> proteinExist = ptms.stream()
                                    .filter(currentMod -> currentMod.getModificationCvTerm()
                                            .getAccession().equalsIgnoreCase(ptm.getAccession()))
                                    .findAny();
                            if (proteinExist.isPresent()) {
                                proteinExist.get().addPosition(proteinPosition, probabilities);
                            } else {
                                CvParam ptmName = new CvParam(ptm.getCvLabel(),
                                        ptm.getAccession(), ptm.getDescription(),
                                        String.valueOf(ptm.getMass()));
                                IdentifiedModification newPTM = new IdentifiedModification(null, null, ptmName, null);
                                newPTM.addPosition(proteinPosition, probabilities);
                                ptms.add(newPTM);
                            }

                            if (position > 0 && position < (item.getSequence().length() + 1)) {
//                                mod.addPosition(position, null);
//                                modifications.add(mod);
//                                log.info(String.valueOf(proteinPosition));
//                                log.info(ptm.getAccession());
                            } else if (position == 0) { //n-term for protein
//                                mod.addPosition(position, null);
//                                modifications.add(mod);
//                                log.info(String.valueOf(proteinPosition));
//                                log.info(ptm.getAccession());

                            }
                        } else {
//                            modifications.add(mod);
                            //if position is not set null is reported
                        }

                    }

                });

            }
        }
        return ptms;

    }

    public void indexSpectraStep() throws Exception {

        long initSpectraStep = System.currentTimeMillis();

        for (MongoPrideAssay mongoPrideAssay : assayList) {

            log.info("indexSpectraStep assay file  -- " + mongoPrideAssay.getFileName());

            Map<String, Object> assayObjects = assayObjectMap.get(mongoPrideAssay.getAccession());
            List<ReportPeptide> highQualityPeptides = (List<ReportPeptide>) assayObjects.get("highQualityPeptides");

            List<ReportPeptide> peptides;
            if (highQualityPeptides.size() > 0)
                peptides = highQualityPeptides;
            else
                peptides = (List<ReportPeptide>) assayObjects.get("allPeptides");

            PIAModeller modeller = (PIAModeller) assayObjects.get("modeller");

            if (modeller != null && peptides.size() > 0) {

                List<SpectraData> spectrumFiles = modeller.getSpectraData()
                        .entrySet().stream().map(Map.Entry::getValue)
                        .collect(Collectors.toList());

                AtomicInteger totalPSM = new AtomicInteger();
                AtomicInteger errorDeltaPSM = new AtomicInteger();

                Optional<MongoAssayFile> assayResultFile = mongoPrideAssay.getAssayFiles()
                        .stream().filter(x -> x.getFileCategory()
                                .getValue().equalsIgnoreCase("RESULT")).findFirst();

                if (assayResultFile.isPresent() && (spectrumFiles.size() > 0
                        || assayResultFile.get().getFileCategory().getAccession().equalsIgnoreCase("PRIDE:1002848"))) {

                    JmzReaderSpectrumService service = null;
                    List<Triple<String, SpectraData, SubmissionPipelineConstants.FileType>> mongoRelatedFiles = new ArrayList<>(spectrumFiles.size());


                    if (spectrumFiles.size() == 0) {
                        throw new Exception("No spectra file found");
                    }
                    List<Triple<String, SpectraData, SubmissionPipelineConstants.FileType>> finalMongoRelatedFilesFiltered = mongoRelatedFiles;
//                                    List<MongoPrideFile> fileNamesFiltered = mongoPrideFiles.stream().filter(x -> {
//                                        boolean isNotFile = false;
//                                        for( SpectraData fileInResult: spectrumFiles){
//                                            String location = fileInResult.getLocation();
//                                            String name = fileInResult.getName();
//                                            if (location != null){
//                                                Path p = Paths.get(location);
//                                                String file = p.getFileName().toString();
//                                                if (x.getFileName().contains(file)){
//                                                    isNotFile = true;
//                                                    String filePath  = String.join("/",projectPath, "internal", file);
//                                                    finalMongoRelatedFilesFiltered.add(new Triple<>(filePath, fileInResult, SubmissionPipelineConstants.FileType.MZML));
//                                                    SubmissionPipelineConstants.copyUncompressSpectraFile(x.getFileName(), file, projectPath);
//                                                }
//                                            }
//                                        }
//                                        return isNotFile;
//                                    }).collect(Collectors.toList());

                    service = JmzReaderSpectrumService.getInstance(finalMongoRelatedFilesFiltered);




                    JmzReaderSpectrumService finalService = service;
                    List<Triple<String, SpectraData, SubmissionPipelineConstants.FileType>> finalMongoRelatedFiles = finalMongoRelatedFilesFiltered;

                    peptides.forEach(peptide -> peptide.getPSMs().forEach(psm -> {
                        try {
                            PeptideSpectrumMatch spectrum = null;
                            if (psm instanceof ReportPSM)
                                spectrum = ((ReportPSM) psm).getSpectrum();

                            totalPSM.set(totalPSM.get() + 1);

                            PeptideSpectrumMatch finalSpectrum = spectrum;
                            System.out.println(finalSpectrum.getSourceID());

                            Spectrum fileSpectrum = null;
                            String spectrumFile = null;
                            String fileName = null;
                            Optional<Triple<String, SpectraData, SubmissionPipelineConstants.FileType>> refeFile = null;
                            String usi = null;
                            String spectraUsi = null;

                            if (spectrumFiles.size() > 0) {
                                refeFile = finalMongoRelatedFiles.stream()
                                        .filter(x -> x.getSecond().getId()
                                                .equalsIgnoreCase(finalSpectrum.getSpectrumIdentification()
                                                        .getInputSpectra().get(0).getSpectraDataRef()))
                                        .findFirst();
                                spectrumFile = refeFile.get().getFirst();
                                String spectrumId = SubmissionPipelineConstants.getSpectrumId(refeFile.get().getSecond(), (ReportPSM) psm);
                                fileSpectrum = finalService.getSpectrumById(spectrumFile, spectrumId);

                                usi = SubmissionPipelineConstants.buildUsi(projectAccession, refeFile.get(),
                                        (ReportPSM) psm);
                                Path p = Paths.get(refeFile.get().getFirst());
                                fileName = p.getFileName().toString();
                            }

                            spectraUsi = SubmissionPipelineConstants.getSpectraUsiFromUsi(usi);

                            log.info(fileSpectrum.getId() + " " + (psm.getMassToCharge() - fileSpectrum.getPrecursorMZ()));
                            Double[] masses = new Double[fileSpectrum.getPeakList().size()];
                            Double[] intensities = new Double[fileSpectrum.getPeakList().size()];
                            int count = 0;
                            for (Map.Entry entry : fileSpectrum.getPeakList().entrySet()) {
                                masses[count] = (Double) entry.getKey();
                                intensities[count] = (Double) entry.getValue();
                                count++;
                            }

                            Set<CvParam> properties = new HashSet<>();
                            Set<CvParam> psmAttributes = new HashSet<>();

                            for (ScoreModelEnum scoreModel : ScoreModelEnum.values()) {
                                Double scoreValue = psm.getScore(scoreModel.getShortName());
                                if (scoreValue != null && !scoreValue.isNaN()) {
                                    for (CvTermReference ref : CvTermReference.values()) {
                                        if (ref.getAccession().equalsIgnoreCase(scoreModel.getCvAccession())) {
                                            CvParam cv = new CvParam(ref.getCvLabel(), ref.getAccession(),
                                                    ref.getName(), String.valueOf(scoreValue));
                                            properties.add(cv);
                                            if (ref.getAccession().equalsIgnoreCase("MS:1002355")) {
                                                CvParam bestSearchEngine = new CvParam(cv.getCvLabel(), cv.getAccession(), cv.getName(), cv.getValue());
                                                psmAttributes.add(bestSearchEngine);
                                            }

                                        }

                                    }
                                }
                            }

                            // Capturing additional parameters provided by the user.
                            boolean submitterValid = false;
                            for (AbstractParam abstractParam : spectrum.getParams()) {
                                if (abstractParam != null) {
                                    if (abstractParam instanceof uk.ac.ebi.jmzidml.model.mzidml.CvParam) {
                                        uk.ac.ebi.jmzidml.model.mzidml.CvParam cvParam = (uk.ac.ebi.jmzidml.model.mzidml.CvParam) abstractParam;
                                        if (cvParam.getAccession() != null) {
                                            CvParam cv = new CvParam(cvParam.getCvRef(),
                                                    cvParam.getAccession(), cvParam.getName(), cvParam.getValue());
                                            if (cv.getAccession().equalsIgnoreCase("PRIDE:0000511")) {
                                                psmAttributes.add(cv);
                                                if (cv.getValue().equalsIgnoreCase("true"))
                                                    submitterValid = true;
                                            }
                                            properties.add(cv);
                                        }
                                    }
                                }
                            }

                            properties.add(new CvParam(CvTermReference.MS_PIA_PEPTIDE_QVALUE.getCvLabel(),
                                    CvTermReference.MS_PIA_PEPTIDE_QVALUE.getAccession(),
                                    CvTermReference.MS_PIA_PEPTIDE_QVALUE.getName(),
                                    String.valueOf(peptide.getQValue())));
                            System.out.print(peptide.getQValue());

                            double retentionTime = Double.NaN;
                            if (psm.getRetentionTime() != null)
                                retentionTime = psm.getRetentionTime();

                            List<Double> ptmMasses = peptide.getModifications().entrySet()
                                    .stream().map(x -> x.getValue().getMass()).collect(Collectors.toList());
                            double deltaMass = MoleculeUtilities
                                    .calculateDeltaMz(peptide.getSequence(),
                                            spectrum.getMassToCharge(),
                                            spectrum.getCharge(),
                                            ptmMasses);

                            log.info("Delta Mass -- " + deltaMass);

                            if (deltaMass > 0.9) {
                                errorDeltaPSM.set(errorDeltaPSM.get() + 1);
                            }
                            properties.add(new CvParam(CvTermReference.MS_DELTA_MASS.getCvLabel(),
                                    CvTermReference.MS_DELTA_MASS.getAccession(),
                                    CvTermReference.MS_DELTA_MASS.getName(),
                                    String.valueOf(deltaMass))
                            );

                            List<IdentifiedModification> mods = new ArrayList<>();
                            if (psm.getModifications() != null && psm.getModifications().size() > 0)
                                mods = PrideAnalysisAssayService.this.convertPeptideModifications(psm.getModifications()).stream().map(x -> {

                                    CvParam neutralLoss = null;
                                    if (x.getNeutralLoss() != null)
                                        neutralLoss = new CvParam(x.getNeutralLoss().getCvLabel(),
                                                x.getNeutralLoss().getAccession(),
                                                x.getNeutralLoss().getName(), x.getNeutralLoss().getValue());

                                    List<Tuple<Integer, Set<? extends CvParamProvider>>> positionMap = new ArrayList<>();
                                    if (x.getPositionMap() != null && x.getPositionMap().size() > 0)
                                        positionMap = x.getPositionMap().stream()
                                                .map(y -> new Tuple<Integer, Set<? extends CvParamProvider>>(y.getKey(),
                                                        y.getValue().stream()
                                                                .map(z -> new CvParam(z.getCvLabel(),
                                                                        z.getAccession(), z.getName(), z.getValue()))
                                                                .collect(Collectors.toSet())))
                                                .collect(Collectors.toList());

                                    CvParam modCv = null;
                                    if (x.getModificationCvTerm() != null)
                                        modCv = new CvParam(x.getModificationCvTerm().getCvLabel(),
                                                x.getModificationCvTerm().getAccession(),
                                                x.getModificationCvTerm().getName(),
                                                x.getModificationCvTerm().getValue());

                                    Set<CvParamProvider> modProperties = new HashSet<>();

                                    return new IdentifiedModification(neutralLoss, positionMap, modCv, modProperties);
                                }).collect(Collectors.toList());

                            Set<CvParam> validationMethods = (Set<CvParam>) assayObjects.get("validationMethods");

                            boolean isValid = (boolean) assayObjects.get("isValid");

                            int misssedCleavages = ((ReportPSM) psm).getMissedCleavages();
                            if (misssedCleavages == -1){
                                misssedCleavages = uk.ac.ebi.pride.utilities.mol.MoleculeUtilities.calcMissedCleavages(psm.getSequence());
                            }

                            PSMProvider archivePSM = ArchiveSpectrum
                                    .builder()
                                    .projectAccession(projectAccession)
                                    .assayAccession(mongoPrideAssay.getAccession())
                                    .peptideSequence(psm.getSequence())
                                    .isDecoy(psm.getIsDecoy())
                                    .retentionTime(retentionTime)
                                    .msLevel(fileSpectrum.getMsLevel())
                                    .precursorCharge(fileSpectrum.getPrecursorCharge())
                                    .masses(masses)
                                    .numPeaks(intensities.length)
                                    .intensities(intensities)
                                    .properties(properties)
                                    .spectrumFile(spectrumFile)
                                    .modifications(mods)
                                    .precursorMz(fileSpectrum.getPrecursorMZ())
                                    .usi(usi)
                                    .spectrumFile(spectrumFile)
                                    .isValid(isValid)
                                    .missedCleavages(misssedCleavages)
                                    .qualityEstimationMethods(validationMethods.stream()
                                            .map(x -> new CvParam(x.getCvLabel(),
                                                    x.getAccession(), x.getName(), x.getValue()))
                                            .collect(Collectors.toSet()))
                                    .build();

                            PrideMongoPsmSummaryEvidence psmMongo = PrideMongoPsmSummaryEvidence
                                    .builder()
                                    .usi(usi)
                                    .spectraUsi(spectraUsi)
                                    .peptideSequence(psm.getSequence())
                                    .assayAccession(mongoPrideAssay.getAccession())
                                    .isDecoy(psm.getIsDecoy())
                                    .charge(psm.getCharge())
                                    .isValid(isValid)
                                    .projectAccession(projectAccession)
                                    .fileName(fileName)
                                    .additionalAttributes(psmAttributes)
                                    .precursorMass(psm.getMassToCharge())
                                    .modifiedPeptideSequence(SubmissionPipelineConstants
                                            .encodePeptide(psm.getSequence(), psm.getModifications()))
                                    .build();

                            if (isValid || submitterValid) {
                                try {
                                    BackupUtil.write(archivePSM, (BufferedWriter) assayObjects.get("archiveSpectrumBufferedWriter"));
                                    BackupUtil.write(psmMongo, (BufferedWriter) assayObjects.get("psmSummaryEvidenceBufferedWriter"));
//                                                    moleculesService.insertPsmSummaryEvidence(psmMongo);
                                } catch (DuplicateKeyException ex) {
//                                                    moleculesService.savePsmSummaryEvidence(psmMongo);
                                    log.debug("The psm evidence was already in the database -- " + psmMongo.getUsi());
                                }

//                                                pushToS3(archivePSM.getUsi(), archivePSM, 0);

                                Map<Long, List<PeptideSpectrumOverview>> peptideUsi = new HashMap<>();

                                List<PeptideSpectrumOverview> usis = new ArrayList<>();
                                if (peptideUsi.containsKey(peptide.getPeptide().getID())) {
                                    usis = peptideUsi.get(peptide.getPeptide().getID());
                                }
                                usis.add(new PeptideSpectrumOverview(psm.getCharge(), psm.getMassToCharge(), usi));
                                peptideUsi.put(peptide.getPeptide().getID(), usis);
                                assayObjects.put("peptideUsi", peptideUsi);
                            }

                        } catch (Exception e) {
                            log.error(e.getMessage(), e);
                            if (!(e instanceof JMzReaderException))
                                throw new RuntimeException(e);
                        }
                    }));
                }
                log.info("Delta Mass Rate -- " + (errorDeltaPSM.get() / totalPSM.get()));
            }

        }

        taskTimeMap.put(SubmissionPipelineConstants.PrideArchiveStepNames.PRIDE_ARCHIVE_MONGODB_SPECTRUM_UPDATE.getName(),
                System.currentTimeMillis() - initSpectraStep);


    }


//    private void pushToS3(String usi, PSMProvider psm, int retry) throws Exception {
//        try {
//            spectralArchive.writePSM(usi, psm);
//        } catch (com.amazonaws.SdkClientException e) {
//            Thread.sleep(10000);
//            if (retry < PIPELINE_RETRY_LIMIT)
//                pushToS3(usi, psm, retry + 1);
//            else
//                throw new Exception("The S3 is not working properly, multiple retries failed --");
//        }
//    }

    private String getSpectraLocation(SpectraData spectraData) {
        return null;
    }

    private CvParam updateValueOfMongoParamter(CvParam param, CvTermReference cvTerm, Integer value) {
        if (param.getAccession().equalsIgnoreCase(cvTerm.getAccession())) {
            param.setValue(String.valueOf(value));
        }
        return param;
    }


}
