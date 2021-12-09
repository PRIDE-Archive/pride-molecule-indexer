package uk.ac.ebi.pride.archive.indexer.services;


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
import org.apache.commons.io.FilenameUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.stereotype.Service;
import uk.ac.ebi.jmzidml.model.mzidml.AbstractParam;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.pride.archive.dataprovider.common.Tuple;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PSMProvider;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PeptideSpectrumOverview;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModification;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModificationProvider;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParam;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParamProvider;
import uk.ac.ebi.pride.archive.dataprovider.param.Param;
import uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.JmzReaderSpectrumService;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PIAModelerService;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideArchiveWebService;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideFile;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideProject;
import uk.ac.ebi.pride.archive.indexer.utility.BackupUtil;
import uk.ac.ebi.pride.archive.spectra.model.ArchiveSpectrum;
import uk.ac.ebi.pride.mongodb.molecules.model.peptide.PrideMongoPeptideEvidence;
import uk.ac.ebi.pride.mongodb.molecules.model.protein.PrideMongoProteinEvidence;
import uk.ac.ebi.pride.mongodb.molecules.model.psm.PrideMongoPsmSummaryEvidence;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.utilities.term.CvTermReference;
import uk.ac.ebi.pride.utilities.util.MoleculeUtilities;
import uk.ac.ebi.pride.utilities.util.Triple;

import java.io.*;
import java.nio.file.AccessDeniedException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;

@Configuration
@Slf4j
@Service
public class PrideAnalysisAssayService {

    private static final Long MERGE_FILE_ID = 1L;
    private static final Long FILE_ID = 1L;

    @Autowired
    PrideArchiveWebService prideArchiveWebService;

    private PIAModelerService piaModellerService;

    @Value("${productionPath}")
    String productionPath;

    @Value("${qValueThreshold:#{0.01}}")
    private Double qValueThreshold;

    @Value("${qFilterProteinFDR:#{1.0}}")
    private Double qFilterProteinFDR;

    DecimalFormat df = new DecimalFormat("###.#####");

    @Bean
    PIAModelerService getPIAModellerService() {
        piaModellerService = new PIAModelerService();
        return piaModellerService;
    }

    Map<String, List<Param>> sampleProperties;
    Map<String, Set<Param>> globalSampleProperties;

    /**
     * Write the Related files to one Result file
     * @param projectAccession Project Accession
     * @param resultFiles Result files to be analyzed
     * @param outputFile Output file to write the relation
     * @throws IOException
     */
    public void writeRelatedSpectraFiles(String projectAccession, List<String> resultFiles, String outputFile)
            throws IOException {

        Optional<PrideProject> projectOption = prideArchiveWebService.findByAccession(projectAccession);
        if(!projectOption.isPresent())
            throw new IOException("Project not present in the PRIDE WS for accession: " + projectAccession);
        PrideProject project = projectOption.get();

        List<PrideFile> projectFiles = prideArchiveWebService.findFilesByProjectAccession(projectAccession);
        if(projectFiles.size() == 0){
            throw new IOException("Not files found in the PRIDE WS for accession: " + projectAccession);
        }

        List<Triple<Tuple<String, String>, PrideFile, SubmissionPipelineUtils.FileType>> filesRelated = new ArrayList<>();

        resultFiles.stream().forEach( resultFile ->{
            SubmissionPipelineUtils.FileType fileType = SubmissionPipelineUtils.FileType.getFileTypeFromFileName(resultFile);
            boolean isCompressFile = SubmissionPipelineUtils.isCompressedByExtension(resultFile);
            if(fileType != null && !isCompressFile){
                try {
                    PIAModeller modeller = piaModellerService.performProteinInference(resultFile,resultFile,
                            fileType, 1.0, 1.0);
                    Map<String, SpectraData> spectraDataList = modeller.getSpectraData();
                    if(fileType == SubmissionPipelineUtils.FileType.PRIDE)
                        filesRelated.add(new Triple<>(new Tuple<>(resultFile, null), null, null));
                    else{
                        filesRelated.addAll(getFilesRelatedToResultFile(resultFile, spectraDataList, projectFiles));
                    }

                } catch (IOException e) {
                    log.info(String.format("Error reading the file %s with error %s",resultFile,e.getMessage()));
                }
            }else
                log.info("Provided Result File is not a recognized extension or is compressed: " + resultFile);
        });

        try {
            String pattern = "yyyy-MM-dd";
            SimpleDateFormat simpleDateFormat = new SimpleDateFormat(pattern);

            String date = simpleDateFormat.format(projectOption.get().getPublicationDate());
            try (PrintWriter writer = new PrintWriter(
                    Files.newBufferedWriter(Paths.get(outputFile)))) {
                writer.printf("%s\t%s\t%s\t%s\t%s\t%s", "resultFile", "date", "referenceFile", "fileType", "ftpName", "ftp");
                writer.println();
                filesRelated.forEach(x -> {
                    if(x.getSecond() != null){
                        Optional<CvParamProvider> location = x.getSecond().getPublicFileLocations().stream().filter(y -> Objects.equals(y.getAccession(), "PRIDE:0000469")).findFirst();
                        writer.printf("%s\t%s\t%s\t%s\t%s\t%s", x.getFirst().getKey(), date, x.getFirst().getValue(), x.getThird().name(), x.getSecond().getFileName(), location.get().getValue());
                    }else{
                        writer.printf("%s\t%s\t%s\t%s\t%s\t%s", x.getFirst().getKey(), date, null, null, null, null);
                    }
                    writer.println();
                });
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * Create backup directory for a specific project. General folder for all the backup projects.
     * The Backup root folder will contain a subfolder for each project.
     * @param folderOutput Root folder
     * @param projectAccession project accession. This accession will be used to construct the backup path
     * @return Final path folderOutput + projectAccession
     * @throws AccessDeniedException
     */
    private void createBackupDir(String folderOutput, String projectAccession) throws IOException {
        String path = folderOutput;

        if(!new File(folderOutput).isDirectory() || !new File(folderOutput).exists())
            throw new IOException("The provided path do not exists : " + folderOutput);

        if (!path.endsWith(File.separator)) {
            path = folderOutput + File.separator;
        }
        path = path + projectAccession;
        File file = new File(path);
        if (file.exists() && file.isDirectory()) {
            return;
        }
        boolean mkdirs = file.mkdirs();
        if (!mkdirs) {
            throw new AccessDeniedException("Failed to create Dir : " + folderOutput);
        }
    }

    /**
     * Create Backup folders for an specific resultFile and assay.
     * @param prideFile {@link PrideFile} from PRIDE API WS
     * @param assayObjects AssayObjects to store the protein/peptide/psms information
     * @param folderOutput Root folder containing all the backup files
     * @param projectAccession Project assay
     * @return Object Map updated with the {@link BufferedWriter}s for each object
     * @throws IOException
     */
    private Map<String, Object> createBackupFiles(PrideFile prideFile, Map<String, Object> assayObjects, String folderOutput, String projectAccession) throws IOException {

        // Create first the root folder for the project
        createBackupDir(folderOutput, projectAccession);

        log.info("Creating assay file  -- " + prideFile.getFileName());

        final String peptideEvidenceFileName = BackupUtil.getPrideMongoPeptideEvidenceFile(folderOutput, projectAccession, prideFile.getAccession());
        assayObjects.put("peptideEvidenceBufferedWriter", new BufferedWriter(new FileWriter(peptideEvidenceFileName, false)));

        final String proteinEvidenceFileName = BackupUtil.getPrideMongoProteinEvidenceFile(folderOutput, projectAccession, prideFile.getAccession());
        assayObjects.put("proteinEvidenceBufferedWriter", new BufferedWriter(new FileWriter(proteinEvidenceFileName, false)));

        final String archiveSpectrumFileName = BackupUtil.getArchiveSpectrumFile(folderOutput, projectAccession, prideFile.getAccession());
        assayObjects.put("archiveSpectrumBufferedWriter", new BufferedWriter(new FileWriter(archiveSpectrumFileName, false)));

        final String psmSummaryEvidenceFileName = BackupUtil.getPrideMongoPsmSummaryEvidenceFile(folderOutput, projectAccession, prideFile.getAccession());
        assayObjects.put("psmSummaryEvidenceBufferedWriter", new BufferedWriter(new FileWriter(psmSummaryEvidenceFileName, false)));

        return assayObjects;

    }


    public void writeAnalysisOutputFromResultFiles(String projectAccession, List<String> resultFiles, List<String> spectraFiles, String folderOutput) throws IOException {

        Optional<PrideProject> projectOption = prideArchiveWebService.findByAccession(projectAccession);
        if(!projectOption.isPresent())
            throw new IOException("Project not present in the PRIDE WS for accession: " + projectAccession);
        PrideProject project = projectOption.get();

        List<PrideFile> projectFiles = prideArchiveWebService.findFilesByProjectAccession(projectAccession);
        if(projectFiles.size() == 0){
            throw new IOException("Not files found in the PRIDE WS for accession: " + projectAccession);
        }

        initGlobalSampleMetadata(projectOption.get(), spectraFiles);

        resultFiles.stream().forEach( resultFile -> {
            SubmissionPipelineUtils.FileType fileType = SubmissionPipelineUtils.FileType.getFileTypeFromFileName(resultFile);
            boolean isCompressFile = SubmissionPipelineUtils.isCompressedByExtension(resultFile);
            if(fileType == SubmissionPipelineUtils.FileType.PRIDE && !isCompressFile){
                try {
                    Optional<PrideFile> prideFile = findPrideFileInProjectFiles(resultFile, projectFiles);
                    if(prideFile.isPresent()){
                        Map<String, Object> assayObjectMap = analyzeAssayInformationStep(resultFile, prideFile.get(), fileType);
                        List<ReportPSM> allPSMs = (List<ReportPSM>) assayObjectMap.get("allPsms");
                        if(allPSMs.size() == 0)
                            throw new IOException(String.format("The result file analyzed %s don't have peptides", resultFile));
                        assayObjectMap = createBackupFiles(prideFile.get(), assayObjectMap, folderOutput, projectAccession);

                        indexSpectraStep(projectAccession, resultFile, prideFile.get(), assayObjectMap, fileType, spectraFiles);
                        proteinPeptideIndexStep(prideFile.get(), assayObjectMap, projectAccession);

                        closeBackupFiles(assayObjectMap);

                    }else{
                        log.info(String.format("The file %s can be found in the result file list %s",resultFile, projectFiles.toString()));
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

        });
    }

    private void initGlobalSampleMetadata(PrideProject prideProject, List<String> spectraFiles) {
        globalSampleProperties = new HashMap<>();

        spectraFiles.stream().forEach(x-> {
            String fileName = SubmissionPipelineUtils.getFileNameNoExtension(x);
            Set<Param> properties = new HashSet<>(prideProject.getOrganisms()
                    .stream().map(y -> new Param("organism", y.getName()))
                    .collect(Collectors.toSet()));
            globalSampleProperties.put(fileName, properties);
        });
        log.info(globalSampleProperties.toString());
    }

    /**
     * Close buffer writers which contains all the peptides/proteins and psms
     * @param assayObjects Map of array objects
     * @throws IOException
     */
    private void closeBackupFiles(Map<String, Object> assayObjects) throws IOException {
        for(String bufferName: Arrays.asList("peptideEvidenceBufferedWriter","proteinEvidenceBufferedWriter", "archiveSpectrumBufferedWriter","psmSummaryEvidenceBufferedWriter")){
            BufferedWriter bufferedWriter = (BufferedWriter) assayObjects.get(bufferName);
            bufferedWriter.flush();
            bufferedWriter.close();
        }
    }

    /**
     * Find a result file name in the list of PRIDE Web services files for an specific project
     * @param resultFile Result file to be analyzed
     * @param projectFiles List of files for one project from PRIDE API
     * @return uk.ac.ebi.pride.archive.indexer.services.ws.PrideFile
     */
    private Optional<PrideFile> findPrideFileInProjectFiles(String resultFile, List<PrideFile> projectFiles) {
        return projectFiles.stream().filter( x-> {
            String resultFileName = FilenameUtils.getName(resultFile);
            if(x.getFileName().toLowerCase().contains(resultFileName.toLowerCase()))
                return true;
            return false;
        }).findFirst();
    }

    /**
     * Analyze result file to get the list of peptides, psms and proteins identified.
     * @param resultFile Result file, supported formats [PRIDE, MZTAB, MZIDENTML]
     * @param prideFile {@link PrideFile} project File retrieved from the PRIDE Web service
     * @param fileType {@link uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils.FileType}
     * @return Object Map with all psms, peptides and proteins lists
     * @throws IOException
     */
    public Map<String, Object> analyzeAssayInformationStep(String resultFile, PrideFile prideFile, SubmissionPipelineUtils.FileType fileType) throws IOException {

        long initAnalysisAssay = System.currentTimeMillis();

        log.info("Analyzing assay file  -- " + resultFile);

        Map<String, Object> assayObjectMap = new HashMap<>();

        // The first threshold for modeller is not threshold at PSM and Protein level.
        PIAModeller modeller = piaModellerService.performProteinInference(prideFile.getAccession(),resultFile,
                fileType, 1.0, 1.0);

        long nrDecoys = modeller.getPSMModeller().getReportPSMSets().entrySet().stream()
                .filter(entry -> entry.getValue().getIsDecoy())
                .count();

        assayObjectMap.put("isValid", nrDecoys > 0);

        Set<CvParam> validationMethods = new HashSet<>();
        validationMethods.add(new CvParam(CvTermReference.MS_DECOY_VALIDATION_METHOD.getCvLabel(),
                CvTermReference.MS_DECOY_VALIDATION_METHOD.getAccession(), CvTermReference.MS_DECOY_VALIDATION_METHOD.getName(), String.valueOf(nrDecoys > 0)));

        assayObjectMap.put("validationMethods", validationMethods);

        List<AbstractFilter> filters = new ArrayList<>();
        filters.add(RegisteredFilters.PSM_SOURCE_ID_FILTER
                .newInstanceOf(FilterComparator.equal, "index=null", true));

        assayObjectMap.put("allPsms", modeller.getPSMModeller()
                .getFilteredReportPSMs(FILE_ID, filters));
        assayObjectMap.put("allPeptides", modeller.getPeptideModeller()
                .getFilteredReportPeptides(MERGE_FILE_ID, filters));
        assayObjectMap.put("allProteins", modeller.getProteinModeller()
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

        assayObjectMap.put("modeller", modeller);
        assayObjectMap.put("highQualityPsms", highQualityPsms);
        assayObjectMap.put("highQualityPeptides", highQualityPeptides);
        assayObjectMap.put("highQualityProteins", highQualityProteins);
        log.info(String.valueOf(System.currentTimeMillis() - initAnalysisAssay));
        return assayObjectMap;
    }

    public void indexSpectraStep(String projectAccession, String resultFile, PrideFile prideFile,
                                 Map<String, Object> assayObjects,
                                 SubmissionPipelineUtils.FileType fileType,
                                 List<String> spectraFiles) throws Exception {

        long initSpectraStep = System.currentTimeMillis();
        log.info("indexSpectraStep assay file  -- " + assayObjects.get("modeller").toString());

        List<ReportPeptide> highQualityPeptides = (List<ReportPeptide>) assayObjects.get("highQualityPeptides");

        List<ReportPeptide> peptides;
        if (highQualityPeptides.size() > 0)
            peptides = highQualityPeptides;
        else
            peptides = (List<ReportPeptide>) assayObjects.get("allPeptides");

        PIAModeller modeller = (PIAModeller) assayObjects.get("modeller");
        JmzReaderSpectrumService service = null;

        if (modeller != null && peptides.size() > 0) {

            List<SpectraData> spectrumFiles = new ArrayList<>(modeller.getSpectraData().values());

            /** Qvalues and FDR values will be used as the main bestSearchEngine Score **/
            Set<Double> qvalues = peptides.stream().map(ReportPeptide::getQValue).collect(Collectors.toSet());
            Set<Double> fdrValues = peptides.stream().map(x -> x.getFDRScore().getValue()).collect(Collectors.toSet());

            AtomicInteger totalPSM = new AtomicInteger();
            AtomicInteger errorDeltaPSM = new AtomicInteger();

            if(fileType == SubmissionPipelineUtils.FileType.PRIDE)
                service = JmzReaderSpectrumService.getInstance(Collections.singletonList(new uk.ac.ebi.pride.utilities.util.Tuple<>(resultFile, fileType)));
            else{
                List<uk.ac.ebi.pride.utilities.util.Tuple<String, SubmissionPipelineUtils.FileType>> mongoRelatedFiles = new ArrayList<>(spectrumFiles.size());

                if (spectrumFiles.size() == 0) {
                    throw new Exception("No spectra file found");
                }
                List<uk.ac.ebi.pride.utilities.util.Tuple<String, SubmissionPipelineUtils.FileType>> finalMongoRelatedFilesFiltered = mongoRelatedFiles;
                service = JmzReaderSpectrumService.getInstance(finalMongoRelatedFilesFiltered);
                List<uk.ac.ebi.pride.utilities.util.Tuple<String, SubmissionPipelineUtils.FileType>> finalMongoRelatedFiles = finalMongoRelatedFilesFiltered;
            }

            JmzReaderSpectrumService finalService = service;
            Map<String, List<PeptideSpectrumOverview>> peptideUsi = new HashMap<>();
            peptides.forEach(peptide -> peptide.getPSMs().forEach(psm -> {

                try {
                    PeptideSpectrumMatch spectrum = null;
                    if (psm instanceof ReportPSM)
                        spectrum = ((ReportPSM) psm).getSpectrum();

                    totalPSM.set(totalPSM.get() + 1);

                    PeptideSpectrumMatch finalSpectrum = spectrum;

                    Spectrum fileSpectrum = null;
                    String usi = null;
                    String spectraUsi = null;
                    String fileName = null;
                    if(fileType == SubmissionPipelineUtils.FileType.PRIDE){
                        fileSpectrum = finalService.getSpectrumById(resultFile, finalSpectrum.getSourceID());
                        fileName = FilenameUtils.getName(resultFile);
                        usi = SubmissionPipelineUtils.buildUsi(projectAccession, fileName, (ReportPSM) psm);
                        spectraUsi = SubmissionPipelineUtils.getSpectraUsiFromUsi(usi);
                    }

                    if(fileSpectrum != null){

                        log.info(fileSpectrum.getId() + " " + (psm.getMassToCharge() - fileSpectrum.getPrecursorMZ()));
                        Double[] masses = new Double[fileSpectrum.getPeakList().size()];
                        Double[] intensities = new Double[fileSpectrum.getPeakList().size()];
                        int count = 0;
                        for (Map.Entry entry : fileSpectrum.getPeakList().entrySet()) {
                            masses[count] = (Double) entry.getKey();
                            intensities[count] = (Double) entry.getValue();
                            count++;
                        }

                        /** Add all scores for the PTMs **/
                        Set<Param> scores = new HashSet<>();
                        for (ScoreModelEnum scoreModel : ScoreModelEnum.values()) {
                            Double scoreValue = psm.getScore(scoreModel.getShortName());
                            if (scoreValue != null && !scoreValue.isNaN() && scoreValue != 0.0
                                    && !Objects.equals(scoreModel.getCvAccession(), "MS:1002355")
                                    && !Objects.equals(scoreModel.getCvAccession(), "MS:1002354")) {
                                for (CvTermReference ref : CvTermReference.values()) {
                                    if (ref.getAccession().equalsIgnoreCase(scoreModel.getCvAccession()))
                                        scores.add(new Param(ref.getName(),String.valueOf(scoreValue)));
                                }
                            }
                        }

                        /** Capture best search engine score **/
                        double piaQvalue = SubmissionPipelineUtils.getQValueLower(psm.getQValue(), qvalues);
                        Param bestSearchEngineScore = new Param(CvTermReference.MS_PIA_PSM_LEVEL_QVALUE.getName(), String.valueOf(piaQvalue));
                        scores.add(bestSearchEngineScore);


                        // Capturing additional parameters provided by the user.
                        Set<Param> properties = new HashSet<>();
                        for (AbstractParam abstractParam : spectrum.getParams()) {
                            if (abstractParam != null) {
                                if (abstractParam instanceof uk.ac.ebi.jmzidml.model.mzidml.CvParam) {
                                    uk.ac.ebi.jmzidml.model.mzidml.CvParam cvParam = (uk.ac.ebi.jmzidml.model.mzidml.CvParam) abstractParam;
                                    if (cvParam.getAccession() != null)
                                        properties.add(new Param(cvParam.getName(), cvParam.getValue()));
                                }
                            }
                        }

                        double piaFDR = SubmissionPipelineUtils.getQValueLower(((ReportPSM) psm).getFDRScore().getValue(), fdrValues);
                        scores.add(new Param(CvTermReference.MS_PIA_PSM_LEVEL_FDRSCORE.getName(), String.valueOf(piaFDR)));
                        log.info(String.valueOf(piaQvalue));

                        double retentionTime = Double.NaN;
                        if (psm.getRetentionTime() != null)
                            retentionTime = psm.getRetentionTime();

                        List<Double> ptmMasses = peptide.getModifications().values()
                                .stream().map(Modification::getMass).collect(Collectors.toList());
                        double deltaMass = MoleculeUtilities
                                .calculateDeltaMz(peptide.getSequence(),
                                        spectrum.getMassToCharge(),
                                        spectrum.getCharge(),
                                        ptmMasses);

                        log.info("Delta Mass -- " + deltaMass);

                        if (deltaMass > 0.9) {
                            errorDeltaPSM.set(errorDeltaPSM.get() + 1);
                        }else if (deltaMass > 10){
                            throw new Exception(String.format("The delta mass for the following PSM --- %s is over 10", usi));
                        }

                        /** PTMs parsing **/
                        List<IdentifiedModification> mods = new ArrayList<>();
                        if (psm.getModifications() != null && psm.getModifications().size() > 0)
                            mods = convertPeptideModifications(psm.getModifications()).stream().map(x -> {

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
                                .assayAccession(prideFile.getAccession())
                                .peptideSequence(psm.getSequence())
                                .isDecoy(psm.getIsDecoy())
                                .retentionTime(retentionTime)
                                .msLevel(fileSpectrum.getMsLevel())
                                .precursorCharge(fileSpectrum.getPrecursorCharge())
                                .masses(masses)
                                .numPeaks(intensities.length)
                                .intensities(intensities)
                                .spectrumFile(fileName)
                                .modifications(mods)
                                .precursorMz(fileSpectrum.getPrecursorMZ())
                                .usi(usi)
                                .isValid(isValid)
                                .missedCleavages(misssedCleavages)
                                .qualityEstimationMethods(validationMethods.stream()
                                        .map(x -> new Param(x.getName(), x.getValue()))
                                        .collect(Collectors.toSet()))
                                .properties(properties)
                                .bestSearchEngineScore(bestSearchEngineScore)
                                .scores(scores)
                                .build();

                        PrideMongoPsmSummaryEvidence psmMongo = PrideMongoPsmSummaryEvidence
                                .builder()
                                .usi(usi)
                                .spectraUsi(spectraUsi)
                                .peptideSequence(psm.getSequence())
                                .assayAccession(prideFile.getAccession())
                                .isDecoy(psm.getIsDecoy())
                                .charge(psm.getCharge())
                                .isValid(isValid)
                                .projectAccession(projectAccession)
                                .fileName(fileName)
                                .scores(scores)
                                .bestSearchEngineScore(bestSearchEngineScore)
                                .precursorMass(psm.getMassToCharge())
                                .modifiedPeptideSequence(SubmissionPipelineUtils
                                        .encodePeptide(psm.getSequence(), psm.getModifications()))
                                .build();

                        log.info(psmMongo.toString());
                        log.info(archivePSM.toString());

                        try {
                            BackupUtil.write(archivePSM, (BufferedWriter) assayObjects.get("archiveSpectrumBufferedWriter"));
                            BackupUtil.write(psmMongo, (BufferedWriter) assayObjects.get("psmSummaryEvidenceBufferedWriter"));
                        } catch (Exception ex) {
                            log.debug("Error writing the PSMs in the files -- " + psmMongo.getUsi());
                        }

                        List<PeptideSpectrumOverview> usis = new ArrayList<>();
                        if (peptideUsi.containsKey(peptide.getStringID())) {
                            usis = peptideUsi.get(peptide.getStringID());
                        }
                        usis.add(new PeptideSpectrumOverview(psm.getCharge(), psm.getMassToCharge(), usi));
                        peptideUsi.put(peptide.getStringID(), usis);
                    }else{
                        log.info(String.format("The following spectrum ID is not found in the PRIDE XML -- %s", finalSpectrum.getSourceID()));
                    }
                } catch (Exception e) {
                    log.error(e.getMessage(), e);
                    if (!(e instanceof JMzReaderException))
                        throw new RuntimeException(e);
                }
            }));
            assayObjects.put("peptideUsi", peptideUsi);

            log.info("Delta Mass Rate -- " + (errorDeltaPSM.get() / totalPSM.get()));
            log.info(String.valueOf(System.currentTimeMillis() - initSpectraStep));
        }
    }


    /**
     * This method use the SpectraDataFile to find the list of files in the PRIDE WS that correspond to those referenced
     * files.
     * @param resultFile Result file under analysis
     * @param spectraDataList List of {@link SpectraData} objects
     * @param projectFiles Project Files from the PRIDE WS
     * @return Triple
     */
    private List<Triple<Tuple<String, String>, PrideFile, SubmissionPipelineUtils.FileType>> getFilesRelatedToResultFile(String resultFile,
                                                                                                                         Map<String, SpectraData> spectraDataList,
                                                                                                                         List<PrideFile> projectFiles) {
        return projectFiles.stream().map(x -> {
                    for( SpectraData fileInResult: spectraDataList.values()){
                        String location = fileInResult.getLocation();
                        if (location != null){
                            String file = FilenameUtils.getName(location);
                            if (x.getFileName().contains(file)){
                                SubmissionPipelineUtils.FileType fileType = SubmissionPipelineUtils
                                        .FileType.getFileTypeFromFileName(file);
                                return new Triple<>(new Tuple<>(resultFile, file), x, fileType);
                            }
                        }
                    }
                    return null;
                }).filter(Objects::nonNull)
                .collect(Collectors.toList());
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

    public void proteinPeptideIndexStep(PrideFile prideFile, Map<String, Object> assayObjects, String projectAccession) throws Exception {
        long initInsertPeptides = System.currentTimeMillis();

        if (assayObjects.get("modeller") != null) {
            log.info("proteinPeptideIndexStep assay file  -- " + assayObjects.get("modeller").toString());

            List<ReportPeptide> highQualityPeptides = (List<ReportPeptide>) assayObjects.get("highQualityPeptides");
            List<ReportProtein> highQualityProteins = (List<ReportProtein>) assayObjects.get("highQualityProteins");
            List<ReportPSM> highQualityPsms = (List<ReportPSM>) assayObjects.get("highQualityPsms");
            List<ReportPeptide> allPeptides = (List<ReportPeptide>) assayObjects.get("allPeptides");
            List<ReportProtein> allProteins = (List<ReportProtein>) assayObjects.get("allProteins");

            Set<String> proteinIds = new HashSet<>();
            Set<String> peptideSequences = new HashSet<>();

            List<ReportPeptide> peptides;
            List<ReportProtein> proteins;


            if (highQualityPeptides.size() > 0 && highQualityPsms.size() > 0 && highQualityProteins.size() > 0) {
                peptides = highQualityPeptides;
                proteins = highQualityProteins;
            } else {
                peptides = allPeptides;
                proteins = allProteins;
            }

            Set<Double> qValues = proteins.stream().map(ReportProtein::getQValue).collect(Collectors.toSet());
            Set<Double> fdrValues = proteins.stream().map(ReportProtein::getFDR).collect(Collectors.toSet());

            Set<Double> peptideQValues = peptides.stream().map(ReportPeptide::getQValue).collect(Collectors.toSet());
            Set<Double> peptideFDRs = peptides.stream().map(x -> x.getScore("peptide_fdr_score")).collect(Collectors.toSet());

            List<ReportPeptide> finalPeptides = peptides;
            for (ReportProtein protein : proteins) {
                String proteinSequence = protein.getRepresentative().getDbSequence();
                String proteinAccession = protein.getRepresentative().getAccession();
                Set<String> proteinGroups = protein.getAccessions()
                        .stream().map(Accession::getAccession)
                        .collect(Collectors.toSet());

                List<IdentifiedModificationProvider> proteinPTMs = new ArrayList<>(convertProteinModifications(
                        proteinAccession, protein.getPeptides()));

                Set<Param> scores = new HashSet<>();
                Param bestSearchEngineScore = null;


                if (Double.isFinite(protein.getQValue()) && !Double.isNaN(protein.getQValue())) {

                    Double qValue = SubmissionPipelineUtils.getQValueLower(protein.getQValue(), qValues);
                    Param score = new Param(CvTermReference.MS_PIA_PROTEIN_GROUP_QVALUE.getName(), String.valueOf(qValue));
                    scores.add(score);
                }

                if (protein.getScore() != null && !protein.getScore().isNaN()) {
                    String value = df.format(protein.getScore());
                    bestSearchEngineScore = new Param(CvTermReference.MS_PIA_PROTEIN_SCORE.getName(), value);
                    scores.add(bestSearchEngineScore);
                }

                if (Double.isFinite(protein.getFDR()) && !Double.isNaN(protein.getFDR())) {

                    Double fdr = SubmissionPipelineUtils.getQValueLower(protein.getFDR(), fdrValues);
                    Param score = new Param(CvTermReference.MS_FDR_PROTEIN.getName(), String.valueOf(fdr));
                    scores.add(score);
                }


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
                        .assayAccession(prideFile.getAccession())
                        .isValid(isValid)
                        .numberPeptides(protein.getPeptides().size())
                        .numberPSMs(protein.getNrPSMs())
                        .scores(scores)
                        .bestSearchEngineScore(bestSearchEngineScore)
                        .sequenceCoverage(protein.getCoverage(proteinAccession))
                        .qualityEstimationMethods(validationMethods.stream().map(x -> new Param(x.getName(), x.getValue())).collect(Collectors.toSet()))
                        .build();

                try {
                    BackupUtil.write(proteinEvidence, (BufferedWriter) assayObjects.get("proteinEvidenceBufferedWriter"));
                }catch (Exception e) {
                    log.error(e.getMessage(), e);
                    throw new Exception(e);
                }
                indexPeptideByProtein(protein, finalPeptides, (Map<Long, List<PeptideSpectrumOverview>>) assayObjects.get("peptideUsi")
                        , (Set<CvParam>) assayObjects.get("validationMethods"),
                        (BufferedWriter) assayObjects.get("peptideEvidenceBufferedWriter"), prideFile.getAccession(), isValid,
                        projectAccession, peptideQValues, peptideFDRs);
            }
        }
    }

    /**
     * This method index all the highQualityPeptides that identified a protein into the mongoDB
     *
     * @param protein  Identified Protein
     * @param peptides Collection of identified highQualityPeptides in the experiment
     */
    private void indexPeptideByProtein(ReportProtein protein, List<ReportPeptide> peptides,
                                       Map<Long, List<PeptideSpectrumOverview>> peptideUsi,
                                       Set<CvParam> validationMethods,
                                       BufferedWriter peptideEvidenceBufferedWriter,
                                       String assayAccession, boolean isValid,
                                       String projectAccession, Set<Double> peptidesQValues,
                                       Set<Double> peptideFDRs) throws Exception {


        for (ReportPeptide peptide : protein.getPeptides()) {

            Optional<ReportPeptide> firstPeptide = peptides.stream()
                    .filter(globalPeptide -> globalPeptide
                            .getStringID()
                            .equalsIgnoreCase(peptide.getStringID()))
                    .findFirst();

            if (firstPeptide.isPresent()) {

                Param bestSearchEngine = null;
                Set<Param> scores = new HashSet<>();
                if (!Double.isInfinite(firstPeptide.get().getQValue()) && !Double.isNaN(firstPeptide.get().getQValue())) {
                    Double value = SubmissionPipelineUtils.getQValueLower(firstPeptide.get().getQValue(), peptidesQValues);
                    Param peptideScore = new Param(CvTermReference.MS_PIA_PEPTIDE_QVALUE.getName(), String.valueOf(value));
                    scores.add(peptideScore);
                    bestSearchEngine = peptideScore;
                }


                if (!Double.isInfinite(firstPeptide.get().getScore("peptide_fdr_score"))
                        && !Double.isNaN(firstPeptide.get().getScore("peptide_fdr_score"))) {

                    Double value = firstPeptide.get().getScore("peptide_fdr_score");
                    value = SubmissionPipelineUtils.getQValueLower(value, peptideFDRs);

                    Param peptideScore = new Param(CvTermReference.MS_PIA_PEPTIDE_FDR.getName(), String.valueOf(value));
                    scores.add(peptideScore);
                }

                Set<PeptideSpectrumOverview> usiList = new HashSet<>(peptideUsi.get(firstPeptide.get().getStringID()));

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

                int misssedCleavages = firstPeptide.get().getMissedCleavages();
                if (misssedCleavages == -1){
                    misssedCleavages = uk.ac.ebi.pride.utilities.mol.MoleculeUtilities.calcMissedCleavages(peptide.getSequence());
                }

                PrideMongoPeptideEvidence peptideEvidence = PrideMongoPeptideEvidence
                        .builder()
                        .assayAccession(assayAccession)
                        .proteinAccession(protein.getRepresentative().getAccession())
                        .isDecoy(firstPeptide.get().getIsDecoy())
                        .peptideAccession(SubmissionPipelineUtils
                                .encodePeptide(peptide.getSequence(), peptide.getModifications()))
                        .peptideSequence(peptide.getSequence())
                        .projectAccession(projectAccession)
                        .psmAccessions(usiList)
                        .startPosition(startPosition)
                        .endPosition(endPosition)
                        .missedCleavages(misssedCleavages)
                        .scores(scores)
                        .bestSearchEngineScore(bestSearchEngine)
                        .ptmList(convertPeptideModifications(firstPeptide.get().getModifications()))
                        .isValid(isValid)
                        .psmAccessions(usiList)
                        .qualityEstimationMethods(validationMethods.stream().map(x -> new Param(x.getName(), x.getValue())).collect(Collectors.toSet()))
                        .build();

                try {
                    BackupUtil.write(peptideEvidence, peptideEvidenceBufferedWriter);
                } catch (Exception e) {
                    log.error(e.getMessage(), e);
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
