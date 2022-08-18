package uk.ac.ebi.pride.archive.indexer.services;


import de.mpc.pia.intermediate.Accession;
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
import lombok.extern.slf4j.Slf4j;
import org.apache.commons.io.FilenameUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.stereotype.Service;
import uk.ac.ebi.jmzidml.model.mzidml.AbstractParam;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.pride.archive.dataprovider.common.Tuple;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PSMProvider;
import uk.ac.ebi.pride.archive.dataprovider.data.peptide.PeptideSpectrumOverview;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModification;
import uk.ac.ebi.pride.archive.dataprovider.data.ptm.IdentifiedModificationProvider;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParam;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParamProvider;
import uk.ac.ebi.pride.archive.dataprovider.param.Param;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.model.ArchiveSpectrum;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.model.PrideProteinEvidence;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.model.PridePsmSummaryEvidence;
import uk.ac.ebi.pride.archive.indexer.utility.HashUtils;
import uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.JmzReaderSpectrumService;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PIAModelerService;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideArchiveWebService;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideFile;
import uk.ac.ebi.pride.archive.indexer.services.ws.PrideProject;
import uk.ac.ebi.pride.archive.indexer.utility.BackupUtil;
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
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

@Configuration
@Slf4j
@Service
public class PrideAnalysisAssayService {

    private static final Long MERGE_FILE_ID = 0L;
    private static final Long FILE_ID = 0L;

    @Autowired
    PrideArchiveWebService prideArchiveWebService;

    private PIAModelerService piaModellerService;

    @Value("${productionPath}")
    String productionPath;

    @Value("${qValueThreshold:#{0.01}}")
    private Double qValueThreshold;

    @Value("${qFilterProteinFDR:#{1.0}}")
    private Double qFilterProteinFDR;

    @Value("${minPSMs:#{1000}}")
    private int minPSMs;

    final DecimalFormat df = new DecimalFormat("###.#####");

    @Bean
    PIAModelerService getPIAModellerService() {
        piaModellerService = new PIAModelerService();
        return piaModellerService;
    }

    final Map<String, Set<Param>> sampleProperties = new HashMap<>();
    final Map<String, Set<Param>> globalSampleProperties = new HashMap<>();

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

        resultFiles.forEach(resultFile ->{
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
                        location.ifPresent(cvParamProvider -> writer.printf("%s\t%s\t%s\t%s\t%s\t%s", x.getFirst().getKey(), date, x.getFirst().getValue(), x.getThird().name(), x.getSecond().getFileName(), cvParamProvider.getValue()));
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
     * @param fileAccession file accession, hash
     * @param assayObjects AssayObjects to store the protein/peptide/psms information
     * @param folderOutput Root folder containing all the backup files
     * @param projectAccession Project assay
     * @return Object Map updated with the {@link BufferedWriter}s for each object
     * @throws IOException
     */
    private Map<String, Object> createBackupFiles(String fileAccession, Map<String, Object> assayObjects, String folderOutput, String projectAccession) throws IOException {

        // Create first the root folder for the project
        createBackupDir(folderOutput, projectAccession);

        log.info("Creating assay file  -- " + fileAccession);

        final String proteinEvidenceFileName = BackupUtil.getProteinEvidenceFile(folderOutput, projectAccession, fileAccession);
        assayObjects.put("proteinEvidenceBufferedWriter", new BufferedWriter(new FileWriter(proteinEvidenceFileName, false)));

        final String archiveSpectrumFileName = BackupUtil.getArchiveSpectrumFile(folderOutput, projectAccession, fileAccession);
        assayObjects.put("archiveSpectrumBufferedWriter", new BufferedWriter(new FileWriter(archiveSpectrumFileName, false)));

        final String psmSummaryEvidenceFileName = BackupUtil.getPsmSummaryEvidenceFile(folderOutput, projectAccession, fileAccession);
        assayObjects.put("psmSummaryEvidenceBufferedWriter", new BufferedWriter(new FileWriter(psmSummaryEvidenceFileName, false)));

        return assayObjects;

    }

    public void writeAnalysisOutputFromResultFiles(String projectAccession, List<String> resultFiles, HashSet<String> spectraFiles, Set<String> sampleFiles, String folderOutput, boolean localReanalysis) throws IOException {

        Optional<PrideProject> projectOption = prideArchiveWebService.findByAccession(projectAccession);
        List<PrideFile> projectFiles = new ArrayList<>();
        if (!localReanalysis) {
            projectFiles = prideArchiveWebService.findFilesByProjectAccession(projectAccession);
            if(projectFiles.size() == 0){
                throw new IOException("Not files found in the PRIDE WS for accession: " + projectAccession);
            }
        }

        if(!projectOption.isPresent())
            throw new IOException("Project not present in the PRIDE WS for accession: " + projectAccession);

        initGlobalSampleMetadata(projectOption.get(), spectraFiles, sampleFiles,
                resultFiles.stream()
                        .filter( x-> SubmissionPipelineUtils.FileType.getFileTypeFromFileName(x) == SubmissionPipelineUtils.FileType.PRIDE)
                        .collect(Collectors.toSet()));

        List<PrideFile> finalProjectFiles = projectFiles;

        resultFiles.forEach(resultFile -> {
            SubmissionPipelineUtils.FileType fileType = SubmissionPipelineUtils.FileType.getFileTypeFromFileName(resultFile);
            boolean isCompressFile = SubmissionPipelineUtils.isCompressedByExtension(resultFile);
            if((fileType == SubmissionPipelineUtils.FileType.MZTAB || fileType == SubmissionPipelineUtils.FileType.MZID
            ) && !isCompressFile){
                try {
                    String fileAccession = null;
                    if (!localReanalysis){
                        Optional<PrideFile> prideFile = findPrideFileInProjectFiles(resultFile, finalProjectFiles);
                        fileAccession = prideFile.map(PrideFile::getAccession).orElse(null);
                    }else{
                        fileAccession = HashUtils.calculateSha1Checksum(resultFile);
                    }

                    if(fileAccession != null){
                        Map<String, Object> assayObjectMap = analyzeAssayInformationStep(resultFile, fileAccession, fileType);
                        assayObjectMap = createBackupFiles(fileAccession, assayObjectMap, folderOutput, projectAccession);

                        indexSpectraStep(projectAccession, fileAccession, assayObjectMap, spectraFiles);
                        proteinIndexStep(fileAccession, assayObjectMap, projectAccession);

                        closeBackupFiles(assayObjectMap);

                    }else{
                        log.info(String.format("The file %s can be found in the result file list %s",resultFile, finalProjectFiles));
                    }
                } catch (Exception e) {
                    log.error("Assay -- " + resultFile + " can't be process because of the following error -- " + e.getMessage());
                }
            }
        });
    }

    /**
     * Parse sample properties parse from SDRF files or from the project properties. The sample metadata will be used from
     * the sample files if provided.
     * @param prideProject Pride Project
     * @param spectraFiles SpectraFiles
     * @param sampleProperties Sample files.
     */
    private void initGlobalSampleMetadata(PrideProject prideProject, Set<String> spectraFiles, Set<String> sampleProperties, Set<String> resultFiles) {

        if(sampleProperties != null){
            sampleProperties.forEach( sampleFile -> {
                try {
                    BufferedReader br = new BufferedReader(new FileReader(sampleFile));
                    String line = br.readLine();
                    List<String> header = Arrays.asList(line.split("\t"));
                    List<List<String>> samples = new ArrayList<>();
                    while (line != null) {
                        line = br.readLine();
                        if(line != null)
                            samples.add(Arrays.asList(line.split("\t")));
                    }
                    int fileIndex = header.indexOf("comment[data file]");
                    List<Integer> sampleIndexes = new ArrayList<>();
                    for(int i=0; i < header.size(); i++)
                        if(header.get(i).contains("characteristics"))
                            sampleIndexes.add(i);
                    samples.forEach(sample -> {
                        String sampleFileName = SubmissionPipelineUtils
                                .getFileNameNoExtension(sample.get(fileIndex));
                        Set<Param> fileSamples = new HashSet<>();
                        sampleIndexes.forEach(index -> {
                            String key = header.get(index);
                            key = key.substring(key.indexOf("[")+1,key.indexOf("]"));
                            String value = sample.get(index);
                            Param param = new Param(key, value);
                            fileSamples.add(param);
                        });
                        if(!this.sampleProperties.containsKey(sampleFileName))
                            this.sampleProperties.put(sampleFileName, fileSamples);
                    });
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
        }

        if(resultFiles != null && resultFiles.size() > 0){
            if(spectraFiles == null)
                spectraFiles = new HashSet<>();
            spectraFiles.addAll(resultFiles);
        }

        spectraFiles.forEach(x-> {
            String fileName = SubmissionPipelineUtils.getFileNameNoExtension(x);
            Set<Param> properties = new HashSet<>();
            if(prideProject.getOrganisms() != null)
                properties = prideProject.getOrganisms()
                    .stream()
                    .map(y -> new Param("organism", y.getName()))
                    .collect(Collectors.toSet());
            if(prideProject.getOrganismParts() != null)
                properties.addAll(prideProject.getOrganismParts()
                    .stream()
                    .map(y -> new Param("organism part", y.getName()))
                    .collect(Collectors.toSet()));
            if(prideProject.getDiseases() != null)
                properties.addAll(prideProject.getDiseases()
                    .stream()
                    .map(y -> new Param("disease", y.getName()))
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
        for(String bufferName: Arrays.asList("proteinEvidenceBufferedWriter", "archiveSpectrumBufferedWriter","psmSummaryEvidenceBufferedWriter")){
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
            return x.getFileName().toLowerCase().contains(resultFileName.toLowerCase());
        }).findFirst();
    }

    /**
     * Analyze result file to get the list of peptides, psms and proteins identified.
     * @param resultFile Result file, supported formats [PRIDE, MZTAB, MZIDENTML]
     * @param fileAccession {@link PrideFile} hash file accession for the result file
     * @param fileType {@link uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils.FileType}
     * @return Object Map with all psms, peptides and proteins lists
     * @throws IOException
     */
    public Map<String, Object> analyzeAssayInformationStep(String resultFile, String fileAccession, SubmissionPipelineUtils.FileType fileType) throws IOException {

        long initAnalysisAssay = System.currentTimeMillis();

        log.info("Analyzing assay file  -- " + resultFile);

        Map<String, Object> assayObjectMap = new HashMap<>();

        // The first threshold for modeller is not threshold at PSM and Protein level.
        PIAModeller modeller = piaModellerService.performProteinInference(fileAccession,resultFile,
                fileType, qValueThreshold, qFilterProteinFDR);

        long nrDecoys = modeller.getPSMModeller().getAllFilteredReportPSMs(new ArrayList<>()).stream()
                .filter(ReportPSM::getIsDecoy)
                .count();
        long nrTargets = modeller.getPSMModeller().getAllFilteredReportPSMs(new ArrayList<>()).stream()
                .filter(x -> !x.getIsDecoy())
                .count();

        assayObjectMap.put("isValid", nrDecoys > 0);

        Set<CvParam> validationMethods = new HashSet<>();
        validationMethods.add(new CvParam(CvTermReference.MS_DECOY_VALIDATION_METHOD.getCvLabel(),
                CvTermReference.MS_DECOY_VALIDATION_METHOD.getAccession(), CvTermReference.MS_DECOY_VALIDATION_METHOD.getName(), String.valueOf(nrDecoys > 0)));

        assayObjectMap.put("validationMethods", validationMethods);

        List<AbstractFilter> filters = new ArrayList<>();
        // Remove PSMs with no spectrum reference
        filters.add(RegisteredFilters.PSM_SOURCE_ID_FILTER
                .newInstanceOf(FilterComparator.equal, "index=null", true));

        // Remove peptides less than 7 AAs
        filters.add(RegisteredFilters.PEPTIDE_SEQUENCE_LENGTH.newInstanceOf(FilterComparator.greater_equal, 7, false));

        // Remove evidences at 1% FDR
        filters.add(new PSMScoreFilter(FilterComparator.less_equal, false,
                qValueThreshold, ScoreModelEnum.PSM_LEVEL_Q_VALUE.getShortName()));              // you can also use fdr score here

        // get the FDR filtered highQualityPeptides
        List<ReportPSM> psms = modeller.getPSMModeller().getAllFilteredReportPSMs(filters);

        List<ReportProtein> proteins = modeller.getProteinModeller()
                .getFilteredReportProteins(filters);

        // The Assay to be considered should have decoy psms a minimun number of PSMS of 1000 (default value)
        if (!(nrDecoys > 0 && proteins.size() > 0 && psms.size() > minPSMs)) {
            throw new NumberFormatException("FDR calculation not possible, no decoys present or number of PSms not bigger than -- !!! " + minPSMs);
        }

        assayObjectMap.put("modeller", modeller);
        assayObjectMap.put("psms", psms);
        assayObjectMap.put("proteins", proteins);
        log.info(String.valueOf(System.currentTimeMillis() - initAnalysisAssay));
        return assayObjectMap;
    }

    public void indexSpectraStep(String projectAccession, String fileAccession,
                                 Map<String, Object> assayObjects,
                                 Set<String> spectraFiles) throws Exception {

        long initSpectraStep = System.currentTimeMillis();
        log.info("indexSpectraStep assay file  -- " + assayObjects.get("modeller").toString());

        List<ReportPSM> psms = (List<ReportPSM>) assayObjects.get("psms");
        List<ReportProtein> proteins = (List<ReportProtein>) assayObjects.get("proteins");

        PIAModeller modeller = (PIAModeller) assayObjects.get("modeller");
        JmzReaderSpectrumService service;

        if (modeller != null && psms.size() > 0) {

            List<SpectraData> spectrumFiles = new ArrayList<>(modeller.getSpectraData().values());

            /** Qvalues and FDR values will be used as the main bestSearchEngine Score **/
            Set<Double> qvalues = psms.stream().map(ReportPSM::getQValue).collect(Collectors.toSet());
            Set<Double> fdrValues = psms.stream().map(x -> x.getFDRScore().getValue()).collect(Collectors.toSet());

            AtomicInteger totalPSM = new AtomicInteger();
            AtomicInteger errorDeltaPSM = new AtomicInteger();

            List<uk.ac.ebi.pride.utilities.util.Triple<String, SpectraData, SubmissionPipelineUtils.FileType>> relatedFiles = null;
            relatedFiles = getRelatedFiles(spectrumFiles, spectraFiles);
            if (spectrumFiles.size() == 0) {
                    throw new Exception("No spectra file found");
            }
            service = JmzReaderSpectrumService.getInstance(relatedFiles.stream()
                    .map(x-> new uk.ac.ebi.pride.utilities.util.Tuple<>(x.getFirst(), x.getThird()))
                    .collect(Collectors.toList())
            );

            JmzReaderSpectrumService finalService = service;
            Map<String, List<PeptideSpectrumOverview>> proteinToPsms = new HashMap<>();
            List<Triple<String, SpectraData, SubmissionPipelineUtils.FileType>> finalRelatedFiles = relatedFiles;
            psms.forEach(psm -> {
                try {
                    PeptideSpectrumMatch spectrum = null;
                    if (psm != null) spectrum = psm.getSpectrum();

                    totalPSM.set(totalPSM.get() + 1);

                    PeptideSpectrumMatch finalSpectrum = spectrum;

                    Spectrum fileSpectrum = null;
                    String spectraUsi;
                    String fileName;
                    Triple<String, String, SubmissionPipelineUtils.FileType> spectrumID = SubmissionPipelineUtils.getSpectrumId(finalRelatedFiles, finalSpectrum);
                    if(spectrumID.getThird() == SubmissionPipelineUtils.FileType.MGF){
                        fileSpectrum = finalService.getSpectrumByIndex(spectrumID.getFirst(), spectrumID.getSecond());
                    }else if(spectrumID.getThird() == SubmissionPipelineUtils.FileType.MZML){
                        fileSpectrum = finalService.getSpectrumById(spectrumID.getFirst(), spectrumID.getSecond());
                    }
                    fileName = FilenameUtils.getName(spectrumID.getFirst());
                    String usi = SubmissionPipelineUtils.buildUsi(projectAccession, fileName, psm, spectrumID.getSecond(), spectrumID.getThird());

                    spectraUsi = SubmissionPipelineUtils.getSpectraUsiFromUsi(usi);

                    Set<Param> localSampleProperties = new HashSet<>();
                    String fileNameNoExtension = SubmissionPipelineUtils.getFileNameNoExtension(fileName);
                    if(sampleProperties.containsKey(fileNameNoExtension))
                        localSampleProperties = sampleProperties.get(fileNameNoExtension);
                    else if(globalSampleProperties.containsKey(fileNameNoExtension))
                        localSampleProperties = globalSampleProperties.get(fileNameNoExtension);

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

                        double piaFDR = SubmissionPipelineUtils.getQValueLower(psm.getFDRScore().getValue(), fdrValues);
                        scores.add(new Param(CvTermReference.MS_PIA_PSM_LEVEL_FDRSCORE.getName(), String.valueOf(piaFDR)));
                        log.info(String.valueOf(piaQvalue));

                        double retentionTime = Double.NaN;
                        if (psm.getRetentionTime() != null)
                            retentionTime = psm.getRetentionTime();

                        List<Double> ptmMasses = psm.getModifications().values()
                                .stream().map(Modification::getMass).collect(Collectors.toList());
                        double deltaMass = MoleculeUtilities
                                .calculateDeltaMz(psm.getSequence(),
                                        spectrum.getMassToCharge(),
                                        spectrum.getCharge(),
                                        ptmMasses);

//                        log.info("Delta Mass -- " + deltaMass);

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

                        int misssedCleavages = psm.getMissedCleavages();
                        if (misssedCleavages == -1){
                            misssedCleavages = uk.ac.ebi.pride.utilities.mol.MoleculeUtilities.calcMissedCleavages(psm.getSequence());
                        }

                        PSMProvider archivePSM = ArchiveSpectrum
                                .builder()
                                .projectAccession(projectAccession)
                                .assayAccession(fileAccession)
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
                                .proteinAccessions(psm.getAccessions().stream().map(x -> x.getAccession()).collect(Collectors.toList()))
                                .qualityEstimationMethods(validationMethods.stream()
                                        .map(x -> new Param(x.getName(), x.getValue()))
                                        .collect(Collectors.toSet()))
                                .properties(properties)
                                .bestSearchEngineScore(bestSearchEngineScore)
                                .scores(scores)
                                .sampleProperties(localSampleProperties)
                                .build();

                        PridePsmSummaryEvidence psmMongo = PridePsmSummaryEvidence
                                .builder()
                                .usi(usi)
                                .spectraUsi(spectraUsi)
                                .peptideSequence(psm.getSequence())
                                .assayAccession(fileAccession)
                                .isDecoy(psm.getIsDecoy())
                                .charge(psm.getCharge())
                                .isValid(isValid)
                                .projectAccession(projectAccession)
                                .fileName(fileName)
                                .scores(scores)
                                .bestSearchEngineScore(bestSearchEngineScore)
                                .precursorMass(psm.getMassToCharge())
                                .proteinAccessions(psm.getAccessions().stream().map(x -> x.getAccession()).collect(Collectors.toList()))
                                .modifiedPeptideSequence(SubmissionPipelineUtils
                                        .encodePeptide(psm.getSequence(), psm.getModifications()))
                                .sampleProperties(localSampleProperties)
                                .build();

//                        log.info(psmMongo.toString());
//                        log.info(archivePSM.toString());

                        try {
                            BackupUtil.write(archivePSM, (BufferedWriter) assayObjects.get("archiveSpectrumBufferedWriter"));
                            BackupUtil.write(psmMongo, (BufferedWriter) assayObjects.get("psmSummaryEvidenceBufferedWriter"));
                        } catch (Exception ex) {
                            log.debug("Error writing the PSMs in the files -- " + psmMongo.getUsi());
                        }

                        PeptideSpectrumOverview psmOverview = new PeptideSpectrumOverview(psm.getCharge(), psm.getMassToCharge(), usi,psm.getSequence(),SubmissionPipelineUtils.encodePeptide(psm.getSequence(), psm.getModifications()));

                        psm.getAccessions().forEach( x -> {
                            // For some reason for protein accessions for PSMs are not in any of the protein reported proteins.
                            if (findProteinInReports(proteins, x.getAccession()).isPresent()){
                                List<PeptideSpectrumOverview> usis = new ArrayList<>();
                                if(proteinToPsms.containsKey(x.getAccession())){
                                    usis = proteinToPsms.get(x.getAccession());
                                }
                                usis.add(psmOverview);
                                proteinToPsms.put(x.getAccession(), usis);
                            }
                        });


                    }else{
                        log.info(String.format("The following spectrum ID is not found in the PRIDE XML -- %s", finalSpectrum.getSourceID()));
                    }
                } catch (Exception e) {
                    log.error(e.getMessage(), e);
                    if (!(e instanceof JMzReaderException))
                        throw new RuntimeException(e);
                }
            });
            assayObjects.put("proteinToPsms", proteinToPsms);

            log.info("Delta Mass Rate -- " + (errorDeltaPSM.get() / totalPSM.get()));
            log.info(String.valueOf(System.currentTimeMillis() - initSpectraStep));
        }
    }

    /**
     * This function match the spectrum files referenced by the result file in the {@link SpectraData} with the
     * file names provided by the user.
     * @param spectraDataFiles {@link SpectraData} files related by the result file
     * @param spectraFiles List of files provided by the users
     * @return List of Tuple files.
     */
    private List<uk.ac.ebi.pride.utilities.util.Triple<String, SpectraData, SubmissionPipelineUtils.FileType>> getRelatedFiles(List<SpectraData> spectraDataFiles, Set<String> spectraFiles) throws IOException {

        if(Objects.isNull(spectraDataFiles) || Objects.isNull(spectraFiles) || spectraDataFiles.size() == 0 || spectraFiles.size() == 0){
            throw new IOException("The number of files provided and SpectraData referenced in the result files are different");
        }

        List<uk.ac.ebi.pride.utilities.util.Triple<String, SpectraData, SubmissionPipelineUtils.FileType>> files = spectraDataFiles.stream().map( x-> {
            for(String filePath: spectraFiles){
                String xFileName = FilenameUtils.getName(filePath);
                String spectraDataFileName = FilenameUtils.getName(x.getLocation());
                if(xFileName.equalsIgnoreCase(spectraDataFileName)){
                    SubmissionPipelineUtils.FileType fileType = SubmissionPipelineUtils.FileType.getFileTypeFromFileName(xFileName);
                    if(fileType != null)
                        return new uk.ac.ebi.pride.utilities.util.Triple<>(filePath, x, fileType);
                }
            }
            return null;
        }).filter(Objects::nonNull).collect(Collectors.toList());

        if (files.size() != spectraDataFiles.size())
            throw new IOException("The number of files provided and SpectraData referenced in the result files are different");

        return files;
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

    private Optional<ReportProtein> findProteinInReports(List<ReportProtein> proteins, String proteinKey){
        return proteins.stream().filter(x -> {
            if (x.getRepresentative().getAccession().equalsIgnoreCase(proteinKey))
                return true;
            else if (x.getSubSets().stream().anyMatch(y -> y.getAccessions().stream().map(Accession::getAccession).anyMatch(z -> z.equalsIgnoreCase(proteinKey)))) {
                return true;
            } else
                return x.getAccessions().stream().anyMatch(y -> y.getAccession().equalsIgnoreCase(proteinKey));
        }).findFirst();

    }
    public void proteinIndexStep(String fileAccession, Map<String, Object> assayObjects, String projectAccession) throws Exception {

        long initInsertPeptides = System.currentTimeMillis();

        if (assayObjects.get("modeller") != null) {
            log.info("proteinIndexStep assay file  -- " + assayObjects.get("modeller").toString());

            List<ReportProtein> proteins = (List<ReportProtein>) assayObjects.get("proteins");
            List<ReportPSM> psms = (List<ReportPSM>) assayObjects.get("psms");
            Map<String, List<PeptideSpectrumOverview>> proteinsToPsms = (Map<String, List<PeptideSpectrumOverview>>) assayObjects.get("proteinToPsms");
            System.out.println(proteinsToPsms.keySet());

            Comparator<ReportProtein> comparator = Comparator.comparing(ReportProtein::getScore);
            proteins.sort(comparator.reversed());

            Set<String> proteinIds = new HashSet<>();
            Set<String> peptideSequences = new HashSet<>();

            Set<Double> qValues = proteins.stream().map(ReportProtein::getQValue).collect(Collectors.toSet());
            Set<Double> fdrValues = proteins.stream().map(ReportProtein::getFDR).collect(Collectors.toSet());
            List<String> accessions = new ArrayList<>();

            for (String proteinKey: proteinsToPsms.keySet()) {
                Optional<ReportProtein> proteinOptional = findProteinInReports(proteins, proteinKey);

                if(!proteinOptional.isPresent()){
                    log.info("Protein from PSms not in Protein reports -- " + proteinKey);
                }else{
                    ReportProtein protein = proteinOptional.get();
                    // For some reasons the ReportProtein list is bigger that the protein keys associated with the PSMs (proteinsToPsms)
                    if (!accessions.contains(proteinKey)){

                        accessions.add(proteinKey);
                        String proteinAccession = proteinKey;
                        Set<String> proteinGroups = protein.getAccessions()
                                .stream().map(Accession::getAccession)
                                .collect(Collectors.toSet());

                        List<String> proteinPTMs = new ArrayList<>(convertProteinModifications(protein.getPeptides()));

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

                        Set<PeptideSpectrumOverview> proteinToPsms = new HashSet<>(proteinsToPsms.get(protein.getRepresentative().getAccession()));
                        proteinToPsms = proteinToPsms.stream().collect(Collectors.collectingAndThen(Collectors.toCollection(() -> new TreeSet<>(Comparator.comparing(PeptideSpectrumOverview::getPeptideSequence))),
                                        HashSet::new));

                        log.info("Protein -- " + proteinAccession + " # PSMs -- " + proteinToPsms.size());


                        proteinIds.add(proteinAccession);
                        protein.getPeptides().forEach(x -> peptideSequences.add(x.getSequence()));

                        Set<CvParam> validationMethods = (Set<CvParam>) assayObjects.get("validationMethods");

                        boolean isValid = (boolean) assayObjects.get("isValid");

                        int nPSMs = proteinToPsms.size();
                        int nPeptides = proteinToPsms.stream().map(PeptideSpectrumOverview::getPeptideSequence).collect(Collectors.toSet()).size();

                        Double coverage = protein.getCoverage(proteinAccession);
                        if(coverage.equals(Double.NaN) || coverage.equals(Double.POSITIVE_INFINITY))
                            coverage = null;

                        String proteinSequence = null;
                        if(protein.getAccessions().stream().anyMatch(y -> y.getAccession().equalsIgnoreCase(proteinAccession))){
                            proteinSequence = protein.getAccessions().stream().filter( y -> y.getAccession().equalsIgnoreCase(proteinAccession)).findAny().get().getDbSequence();
                        }

                        PrideProteinEvidence proteinEvidence = PrideProteinEvidence
                                .builder()
                                .reportedAccession(proteinAccession)
                                .isDecoy(protein.getIsDecoy())
                                .proteinGroupMembers(proteinGroups)
                                .ptms(proteinPTMs)
                                .projectAccession(projectAccession)
                                .assayAccession(fileAccession)
                                .isValid(isValid)
                                .numberPeptides(nPeptides)
                                .numberPSMs(nPSMs)
                                .scores(scores)
                                .bestSearchEngineScore(bestSearchEngineScore)
                                .sequenceCoverage(coverage)
                                .proteinSequence(proteinSequence)
                                .qualityEstimationMethods(validationMethods.stream().map(x -> new Param(x.getName(), x.getValue())).collect(Collectors.toSet()))
                                .psmAccessions(proteinToPsms)
                                .build();

                        try {
                            BackupUtil.write(proteinEvidence, (BufferedWriter) assayObjects.get("proteinEvidenceBufferedWriter"));
                            log.info(String.format("Protein %s -- Number of peptides %s", protein.getRepresentative().getAccession(), nPeptides));
                        }catch (Exception e) {
                            log.error(e.getMessage(), e);
                            throw new Exception(e);
                        }
                    }else{
                        log.info("Protein Accession already in report -- " + protein.getRepresentative().getAccession());
                    }
                }
            }
        }
    }

    private List<String> convertProteinModifications(List<ReportPeptide> peptides) {
        return new ArrayList<>(peptides.parallelStream().map(x -> x.getModifications().values().parallelStream().map(Modification::getDescription).collect(Collectors.toList())).flatMap(List::stream).collect(Collectors.toSet()));
    }

    /**
     * This method index all the highQualityPeptides that identified a protein into the mongoDB
     *  @param protein  Identified Protein
     * @param peptides Collection of identified highQualityPeptides in the experiment
     * @return
     */
//    private boolean indexPeptideByProtein(ReportProtein protein, List<ReportPeptide> peptides,
//                                          Map<String, List<PeptideSpectrumOverview>> peptideUsi,
//                                          Set<CvParam> validationMethods,
//                                          BufferedWriter peptideEvidenceBufferedWriter,
//                                          String assayAccession, boolean isValid,
//                                          String projectAccession, Set<Double> peptidesQValues,
//                                          Set<Double> peptideFDRs) {
//
//        AtomicBoolean indexPeptides = new AtomicBoolean(false);
//
//
//         protein.getPeptides().forEach( originalPeptide -> {
//             Optional<ReportPeptide> peptideOption = peptides.stream()
//                     .filter(x -> x.getStringID()
//                             .equalsIgnoreCase(originalPeptide.getStringID()))
//                     .findAny();
//
//             if(peptideOption.isPresent() && peptideUsi.containsKey(peptideOption.get().getStringID())){
//                 ReportPeptide peptide = peptideOption.get();
//
//                 Param bestSearchEngine = null;
//                 Set<Param> scores = new HashSet<>();
//                 if (!Double.isInfinite(peptide.getQValue()) && !Double.isNaN(peptide.getQValue())) {
//                     Double value = SubmissionPipelineUtils.getQValueLower(peptide.getQValue(), peptidesQValues);
//                     Param peptideScore = new Param(CvTermReference.MS_PIA_PEPTIDE_QVALUE.getName(), String.valueOf(value));
//                     scores.add(peptideScore);
//                     bestSearchEngine = peptideScore;
//                 }
//
//
//                 if (!Double.isInfinite(peptide.getScore("peptide_fdr_score"))
//                         && !Double.isNaN(peptide.getScore("peptide_fdr_score"))) {
//
//                     Double value = peptide.getScore("peptide_fdr_score");
//                     value = SubmissionPipelineUtils.getQValueLower(value, peptideFDRs);
//
//                     Param peptideScore = new Param(CvTermReference.MS_PIA_PEPTIDE_FDR.getName(), String.valueOf(value));
//                     scores.add(peptideScore);
//                 }
//
//                 Set<PeptideSpectrumOverview> usiList = null;
//                 if (peptide.getStringID() != null && peptideUsi.containsKey(peptide.getStringID()))
//                     usiList = new HashSet<>(peptideUsi.get(peptide.getStringID()));
//
//                 int startPosition = 0;
//                 int endPosition = 0;
//
//                 Optional<AccessionOccurrence> occurrence = peptide.getPeptide()
//                         .getAccessionOccurrences()
//                         .stream()
//                         .filter(x -> x.getAccession()
//                                 .getAccession()
//                                 .equalsIgnoreCase(protein.getRepresentative().getAccession()))
//                         .findFirst();
//                 if (occurrence.isPresent()) {
//                     startPosition = occurrence.get().getStart();
//                     endPosition = occurrence.get().getEnd();
//                 } else {
//                     log.info("Position of the corresponding peptide is not present -- " + protein.getRepresentative().getAccession());
//                 }
//
//                 AtomicReference<CvParam> param = new AtomicReference<>(new CvParam(PRIDETools.PrideOntologyConstants
//                         .PRIDE_SUBMITTERS_THERSHOLD.getCvLabel(),
//                         PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getAccession(),
//                         PRIDETools.PrideOntologyConstants.PRIDE_SUBMITTERS_THERSHOLD.getName(), Boolean.toString(false)));
//
//                 int misssedCleavages = peptide.getMissedCleavages();
//                 if (misssedCleavages == -1) {
//                     misssedCleavages = uk.ac.ebi.pride.utilities.mol.MoleculeUtilities.calcMissedCleavages(peptide.getSequence());
//                 }
//
//                 PrideMongoPeptideEvidence peptideEvidence = PrideMongoPeptideEvidence
//                         .builder()
//                         .assayAccession(assayAccession)
//                         .proteinAccession(protein.getRepresentative().getAccession())
//                         .isDecoy(peptide.getIsDecoy())
//                         .peptideAccession(SubmissionPipelineUtils
//                                 .encodePeptide(peptide.getSequence(), peptide.getModifications()))
//                         .peptideSequence(peptide.getSequence())
//                         .projectAccession(projectAccession)
//                         .psmAccessions(usiList)
//                         .startPosition(startPosition)
//                         .endPosition(endPosition)
//                         .missedCleavages(misssedCleavages)
//                         .scores(scores)
//                         .bestSearchEngineScore(bestSearchEngine)
//                         .ptmList(convertPeptideModifications(peptide.getModifications()))
//                         .isValid(isValid)
//                         .psmAccessions(usiList)
//                         .qualityEstimationMethods(validationMethods.stream()
//                                 .map(x -> new Param(x.getName(), x.getValue()))
//                                 .collect(Collectors.toSet()))
//                         .build();
//                 log.info(String.format("Peptide %s -- Number of psms %s", peptide.getStringID(), Objects.requireNonNull(usiList).size()));
//
//                 try {
//                     BackupUtil.write(peptideEvidence, peptideEvidenceBufferedWriter);
//                     indexPeptides.set(true);
//                 } catch (Exception e) {
//                     log.error(e.getMessage(), e);
//                 }
//             }
//         });
//         return indexPeptides.get();
//    }

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

                    if (Objects.equals(peptideEvidence.getAccession().getAccession(), proteinAccession)) {

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

//                            if (position > 0 && position < (item.getSequence().length() + 1)) {
////                                mod.addPosition(position, null);
////                                modifications.add(mod);
////                                log.info(String.valueOf(proteinPosition));
////                                log.info(ptm.getAccession());
//                            } else if (position == 0) { //n-term for protein
////                                mod.addPosition(position, null);
////                                modifications.add(mod);
////                                log.info(String.valueOf(proteinPosition));
////                                log.info(ptm.getAccession());
//
//                            }
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

}
