package uk.ac.ebi.pride.archive.indexer.services;

import lombok.extern.slf4j.Slf4j;
import org.ehcache.Cache;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.stereotype.Service;
import uk.ac.ebi.pride.archive.dataprovider.common.Triple;
import uk.ac.ebi.pride.archive.dataprovider.data.protein.PeptideSpectrumOverview;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.BinaryArchiveSpectrum;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.SummaryArchiveSpectrum;
import uk.ac.ebi.pride.archive.dataprovider.param.CvParam;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PIAModelerService;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PeptidoformClustered;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PrideJsonRandomAccess;
import uk.ac.ebi.pride.archive.indexer.utility.*;
import uk.ac.ebi.pride.utilities.term.CvTermReference;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static uk.ac.ebi.pride.archive.indexer.services.PrideAnalysisAssayService.createBackupDir;

@Configuration
@Slf4j
@Service
public class InferenceService implements Serializable{

    private PIAModelerService piaModellerInference;

    @Value("${productionPath}")
    String productionPath;

    @Value("${qValueThreshold:#{0.01}}")
    private Double qValueThreshold;

    @Value("${qFilterProteinFDR:#{0.01}}")
    private Double qFilterProteinFDR;

    @Value("${peptideLength:#{7}}")
    private Double peptideLength;

    @Autowired
    private PSMClusteringService clusterService;

    static final Map<String, PrintWriter> spectraPartitionWriters = new HashMap<>();

    public static Map<String, String> getInferenceCategories(Map<String, List<String>> peptideToProteins, Set<String> proteins) {
        Collection<List<String>> values = peptideToProteins.values();
        return proteins.stream().collect(Collectors.toMap(String::new, item -> {
              for(List<String> value: values){
                  if(value.contains(item) && value.size() == 1)
                      return "distinguishable";
              }
              return "indistinguishable";
        }));
    }

    @Bean
    PIAModelerService getPiaModellerInference() {
        piaModellerInference = new PIAModelerService();
        return piaModellerInference;
    }

    public static Map<String, Double> getBestQValue(Map<String, List<Triple<String, Double,String>>> proteins){
        Comparator<Triple<String, Double, String>> comparator = Comparator.comparing(Triple::getSecond);
        Map<String, Double> proteinScore = proteins.entrySet().stream().collect(
                Collectors.toMap(Map.Entry::getKey, x -> {
                    List<Triple<String, Double, String>> list = x.getValue();
                    list.sort(comparator);
                    return list.stream()
                            .findFirst()
                            .get()
                            .getSecond();
                }));

        proteinScore.forEach((key, value) -> System.out.println(key + "  " + value));
        return proteinScore;
    }

    public void performProteinInference(String pridePSMPath, String maraclusterResultsPath, String projectAccession,
                                        String reanalysisAccession, String folderOutput) throws Exception {

        AppCacheManager appCacheManager = AppCacheManager.getInstance();

        PrideJsonRandomAccess pridePSMJsonReader = new PrideJsonRandomAccess(pridePSMPath);
        pridePSMJsonReader.parseIndex();

        //The index of the spectrum in the reader is the same (0-based) than the key in the cluster map
        Cache<Integer, Integer> clusters = clusterService.readMaraClusterResults(maraclusterResultsPath);

//        if(clusters.size() != pridePSMJsonReader.getKeys().size())
//            throw new Exception(String.format("The number of spectra in the cluster file (%s) is different than the " +
//                    "number of spectra in the json file (%s)", clusters.size(), pridePSMJsonReader.getKeys().size()));

        Cache<Integer, List<Triple<String, PeptidoformClustered, Double>>> clusterScores = (Cache<Integer, List<Triple<String, PeptidoformClustered, Double>>>) appCacheManager.getPeptidoformCache();
        Cache<Integer, Triple<String, PeptidoformClustered, Double>> filterScores = (Cache<Integer, Triple<String, PeptidoformClustered, Double>>) appCacheManager.getFilterPeptidoformCache();

        int index = 0;
        for (Iterator<Cache.Entry<String, Long>> it = pridePSMJsonReader.getKeys(); it.hasNext(); ) {
            String usi = it.next().getKey();
            BinaryArchiveSpectrum spectrum = pridePSMJsonReader.readArchiveSpectrum(usi);
            Double pcmScore = Double.parseDouble(spectrum.getBestSearchEngineScore().getValue());
            Integer clusterId = clusters.get(index);
            List<Triple<String, PeptidoformClustered, Double>> membersCluster= new ArrayList<>();
            if(clusterScores.containsKey(clusterId))
                membersCluster = clusterScores.get(clusterId);
            membersCluster.add(new Triple<>(spectrum.getUsi(), new PeptidoformClustered(spectrum.getPeptideSequence(), spectrum.getPeptidoform(), spectrum.getIsDecoy()), pcmScore));
            clusterScores.put(clusterId, membersCluster);
            index++;
        }

        //Using cache to delete objects.
        Set<Integer> indexTobeRemoved = new HashSet<>();
        for (Iterator<Cache.Entry<Integer, List<Triple<String, PeptidoformClustered, Double>>>> it = clusterScores.iterator(); it.hasNext(); ) {
            Cache.Entry<Integer, List<Triple<String, PeptidoformClustered, Double>>> score = it.next();

            List<Triple<String, PeptidoformClustered, Double>> peptides = score.getValue();
            Set<String> isoPeptideSet = peptides.stream().map(pep -> StringUtils.makePeptideIsobaric(pep.getSecond().getSequence())).collect(Collectors.toSet());
            if(isoPeptideSet.size() > 1)
                indexTobeRemoved.add(score.getKey());

            Map<PeptidoformClustered, Long> isoPeptideForms = peptides.stream().collect(Collectors.groupingBy(Triple::getSecond, Collectors.counting()));
            Map<PeptidoformClustered, Long> validPeptideForms = new HashMap<>();
            for(Map.Entry entry: isoPeptideForms.entrySet()){
                double ratio = Math.round((Long)entry.getValue()/peptides.size() *100)/100;
                if(ratio > 0.5)
                    validPeptideForms.put((PeptidoformClustered) entry.getKey(), (Long) entry.getValue());
            }

            if (validPeptideForms.size() == 0)
                    indexTobeRemoved.add(score.getKey());

            Triple<String, PeptidoformClustered, Double> resultPeptide = null;
            for(Triple<String, PeptidoformClustered, Double> peptide: peptides){
                if(resultPeptide == null || (peptide.getSecond().equals(resultPeptide.getSecond()) && peptide.getThird() < resultPeptide.getThird()))
                    resultPeptide = peptide;
            }
            if(resultPeptide != null)
                filterScores.put(score.getKey(), resultPeptide);
        }

        String hashAssay = HashUtils.getRandomToken();
        Map<String, Object> assayObjects = new HashMap<>();
        createBackupFiles(assayObjects, folderOutput, projectAccession, hashAssay);

        Cache<String, List<PeptideSpectrumOverview>> proteinToPsms = (Cache<String, List<PeptideSpectrumOverview>>) appCacheManager.getProteinToPsmsCache();
        Map<String, List<Triple<String, Double,String>>> proteinsPSMsScores = new HashMap<>();
        Map<String, List<String>> peptideToProteins = new HashMap<>();
        Map<String, Set<String>> proteinPTMs = new HashMap<>();
        Map<String, List<Boolean>> proteinDecoys = new HashMap<>();
        int psmCount = 1;

        for (Iterator<Cache.Entry<Integer, Triple<String, PeptidoformClustered, Double>>> it = filterScores.iterator(); it.hasNext(); ) {
            boolean flush = (psmCount % 1000) == 0;
            Triple<String, PeptidoformClustered, Double> psm = it.next().getValue();
            try {
                PrintWriter batchBufferWriter = null;
                BinaryArchiveSpectrum archivePSM = pridePSMJsonReader.readArchiveSpectrum(psm.getFirst());
                SummaryArchiveSpectrum psmElastic = SummaryArchiveSpectrum
                            .builder()
                            .usi(archivePSM.getUsi())
                            .spectraUsi(archivePSM.getSpectraUsi())
                            .peptideSequence(archivePSM.getPeptideSequence())
                            .assayAccession(archivePSM.getAssayAccession())
                            .isDecoy(archivePSM.getIsDecoy())
                            .precursorCharge(archivePSM.getPrecursorCharge())
                            .isValid(archivePSM.getIsValid())
                            .projectAccession(archivePSM.getProjectAccession())
                            .reanalysisAccession(archivePSM.getReanalysisAccession())
                            .scores(archivePSM.getScores())
                            .numPeaks(archivePSM.getNumPeaks())
                            .bestSearchEngineScore(archivePSM.getBestSearchEngineScore())
                            .precursorMz(archivePSM.getPrecursorMz())
                            .proteinAccessions(archivePSM.getProteinAccessions())
                            .peptidoform(archivePSM.getPeptidoform())
                            .sampleProperties(archivePSM.getSampleProperties())
                            .build();
                    // Total number of spectrum in ArchiveSpectrum
                    BackupUtil.write(archivePSM, (PrintWriter) assayObjects.get("archiveSpectrumPrintWriter"), flush);
                    // Total number of spectrum in Elastic Search summary.
                    BackupUtil.write(psmElastic, (PrintWriter) assayObjects.get("psmSummaryEvidencePrintWriter"), flush);
                    // Writing in batches.

                    String usi = psm.getFirst();
                    String batchFile = usi.split(":")[2];
                    if(!spectraPartitionWriters.containsKey(batchFile)){
                        String prefix = (String) assayObjects.get("archiveSpectrumFilePrefix");
                        batchBufferWriter = new PrintWriter(new FileWriter(BackupUtil.getArchiveSpectrumFileBatch(prefix, batchFile), false));
                        spectraPartitionWriters.put(batchFile, batchBufferWriter);
                    }else
                        batchBufferWriter = spectraPartitionWriters.get(batchFile);

                    BackupUtil.write(archivePSM, batchBufferWriter, flush);
                    // construction of USI list.
                    PeptideSpectrumOverview psmOverview = new PeptideSpectrumOverview(archivePSM.getPrecursorCharge(),
                            archivePSM.getPrecursorMz(), usi , archivePSM.getPeptideSequence(), SubmissionPipelineUtils.removeChargeState(archivePSM.getPeptidoform()));

                    archivePSM.getProteinAccessions().forEach( x -> {
                        // For some reason for protein accessions for PSMs are not in any of the protein reported proteins.
                        List<PeptideSpectrumOverview> usis = new ArrayList<>();
                        if(proteinToPsms.containsKey(x)){
                            usis = proteinToPsms.get(x);
                        }
                        usis.add(psmOverview);
                        proteinToPsms.put(x, usis);

                        //Get the protein Score
                        List<uk.ac.ebi.pride.archive.dataprovider.common.Triple<String, Double,String>> pcms = new ArrayList<>();
                        Double pcmScore = Double.valueOf(archivePSM.getBestSearchEngineScore().getValue());

                        if(proteinsPSMsScores.containsKey(x)){
                            pcms = proteinsPSMsScores.get(x);
                        }
                        pcms.add(new uk.ac.ebi.pride.archive.dataprovider.common.Triple<>(archivePSM.getPeptidoform(), pcmScore, archivePSM.getUsi()));
                        proteinsPSMsScores.put(x, pcms);

                        // Get protein uniqueness
                        List<String> proteinIds = new ArrayList<>();
                        if(peptideToProteins.containsKey(archivePSM.getPeptidoform()))
                            proteinIds = peptideToProteins.get(archivePSM.getPeptidoform());
                        proteinIds.add(x);
                        peptideToProteins.put(archivePSM.getPeptidoform(),proteinIds);

                        // Get protein decoys
                        List<Boolean> decoys = new ArrayList<>();
                        if(proteinDecoys.containsKey(x))
                            decoys = proteinDecoys.get(x);
                        decoys.add(archivePSM.getIsDecoy());
                        proteinDecoys.put(x, decoys);

                        Set<String> ptms = new HashSet<>();
                        if(proteinPTMs.containsKey(x))
                            ptms = proteinPTMs.get(x);
                        ptms.addAll(archivePSM.getModifications().stream().map(m -> m.getModification().getName()).collect(Collectors.toList()));
                        proteinPTMs.put(x, ptms);
                    });
//                } else {
//                    System.out.println("USI with error -- " + psm.getFirst());
//                }
            } catch (Exception e) {
                log.debug("Error writing the PSMs in the files -- " + psm.getFirst());
                throw new RuntimeException(e);
            }
            psmCount++;
        };

        assayObjects.put("proteinToPsms", proteinToPsms);
        Map<String, Double> proteinScores = InferenceService.getBestQValue(proteinsPSMsScores);
        assayObjects.put("proteinScores", proteinScores);
        Map<String, String> proteinCategories = InferenceService.getInferenceCategories(peptideToProteins, proteinScores.keySet());
        assayObjects.put("proteinStatus", proteinCategories);
        // Get decoys
        Map<String, Boolean> decoyStatus = proteinDecoys.entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, e -> e.getValue().stream().anyMatch(x -> x)));
        assayObjects.put("proteinDecoys", decoyStatus);
        assayObjects.put("proteinPTMs", proteinPTMs);
        Set<CvParam> validationMethods = new HashSet<>();
        validationMethods.add(new CvParam(CvTermReference.MS_DECOY_VALIDATION_METHOD.getCvLabel(),
                CvTermReference.MS_DECOY_VALIDATION_METHOD.getAccession(), CvTermReference.MS_DECOY_VALIDATION_METHOD.getName(), String.valueOf(true)));

        assayObjects.put("validationMethods", validationMethods);
        spectraPartitionWriters.values().forEach(x -> {
            x.flush();
            x.close();
        });

        StreamSupport.stream(
                        Spliterators.spliteratorUnknownSize(proteinToPsms.iterator(), Spliterator.ORDERED), false)
                .forEach( x-> System.out.println("Proteins -- " + x.getKey() + " number of PSMs -- " + x.getValue().size() ));
        PrideAnalysisAssayService.proteinIndexStep(hashAssay, assayObjects, projectAccession, reanalysisAccession);

        for(Object object: assayObjects.values()){
            if (object instanceof PrintWriter){
                PrintWriter PrintWriter = (PrintWriter) object;
                PrintWriter.flush();
                PrintWriter.close();
            }
        }
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

    private Map<String, Object> createBackupFiles(Map<String, Object> assayObjects, String folderOutput, String projectAccession, String assayAccession) throws IOException {

        // Create first the root folder for the project
        createBackupDir(folderOutput, projectAccession);

        log.info("Creating assay file  -- " + projectAccession);

        final String proteinEvidenceFileName = BackupUtil.getProteinEvidenceFile(folderOutput, projectAccession, assayAccession);
        assayObjects.put("proteinEvidenceFileName", proteinEvidenceFileName);
        assayObjects.put("proteinEvidencePrintWriter", new PrintWriter(new FileWriter(proteinEvidenceFileName, false)));

        final String archiveSpectrumFileName = BackupUtil.getArchiveSpectrumFile(folderOutput, projectAccession, assayAccession);
        assayObjects.put("archiveSpectrumFileName", archiveSpectrumFileName);
        assayObjects.put("archiveSpectrumPrintWriter", new PrintWriter(new FileWriter(archiveSpectrumFileName, false)));

        final String archiveSpectrumFilePrefix = BackupUtil.getArchiveSpectrumFilePrefix(folderOutput, projectAccession);
        assayObjects.put("archiveSpectrumFilePrefix", archiveSpectrumFilePrefix);

        final String psmSummaryEvidenceFileName = BackupUtil.getPsmSummaryEvidenceFile(folderOutput, projectAccession, assayAccession);
        assayObjects.put("psmSummaryEvidenceFileName", psmSummaryEvidenceFileName);
        assayObjects.put("psmSummaryEvidencePrintWriter", new PrintWriter(new FileWriter(psmSummaryEvidenceFileName, false)));

        return assayObjects;

    }

}
