package uk.ac.ebi.pride.archive.indexer.services;

import lombok.extern.slf4j.Slf4j;
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
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PrideJsonRandomAccess;
import uk.ac.ebi.pride.archive.indexer.utility.BackupUtil;
import uk.ac.ebi.pride.archive.indexer.utility.HashUtils;
import uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils;
import uk.ac.ebi.pride.utilities.term.CvTermReference;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static uk.ac.ebi.pride.archive.indexer.services.PrideAnalysisAssayService.createBackupDir;

@Configuration
@Slf4j
@Service
public class InferenceService {

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

    static final Map<String, BufferedWriter> spectraPartitionWriters = new HashMap<>();

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

        PrideJsonRandomAccess pridePSMJsonReader = new PrideJsonRandomAccess(pridePSMPath);
        pridePSMJsonReader.parseIndex();

        //The index of the spectrum in the reader is the same (0-based) than the key in the cluster map
        Map<Integer, Integer> clusters = clusterService.readMaraClusterResults(maraclusterResultsPath);

        if(clusters.size() != pridePSMJsonReader.getKeys().size())
            throw new Exception(String.format("The number of spectra in the cluster file (%s) is different than the " +
                    "number of spectra in the json file (%s)", clusters.size(), pridePSMJsonReader.getKeys().size()));

        Map<Integer, List<Triple<String, PeptidoformClustered, Double>>> clusterScores = new HashMap<>();
        int index = 0;
        for(String usi: pridePSMJsonReader.getKeys()){
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

        Map<Integer, Optional<Triple<String, PeptidoformClustered, Double>>> clusteredFilter = clusterScores.entrySet()
                .parallelStream()
                .collect(Collectors.toMap(Map.Entry::getKey, e -> {
                    List<Triple<String, PeptidoformClustered, Double>> peptides = e.getValue();
                    Set<String> isoPeptideSet = peptides.stream().map(pep -> makePeptideIsobaric(pep.getSecond().getSequence())).collect(Collectors.toSet());
                    if(isoPeptideSet.size() > 1)
                        return Optional.empty();

                    Map<PeptidoformClustered, Long> isoPeptideForms = peptides.stream().collect(Collectors.groupingBy(Triple::getSecond, Collectors.counting()));
                    Map<PeptidoformClustered, Long> validPeptideForms = new HashMap<>();
                    for(Map.Entry entry: isoPeptideForms.entrySet()){
                        double ratio = Math.round((Long)entry.getValue()/peptides.size() *100)/100;
                        if(ratio > 0.5)
                            validPeptideForms.put((PeptidoformClustered) entry.getKey(), (Long) entry.getValue());
                    }
                    if (validPeptideForms.size() == 0)
                            return Optional.empty();
                    Triple<String, PeptidoformClustered, Double> resultPeptide = null;
                    for(Triple<String, PeptidoformClustered, Double> peptide: peptides){
                        if(resultPeptide == null || (peptide.getSecond().equals(resultPeptide.getSecond()) && peptide.getThird() < resultPeptide.getThird()))
                            resultPeptide = peptide;
                    }
                    return resultPeptide!= null ? Optional.of(resultPeptide):Optional.empty();
                }));
        clusteredFilter = clusteredFilter.entrySet().stream()
                .filter(e -> e.getValue().isPresent())
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        String hashAssay = HashUtils.sha1InObject(clusteredFilter);
        Map<String, Object> assayObjects = new HashMap<>();
        createBackupFiles(assayObjects, folderOutput, projectAccession, hashAssay);

        Map<String, List<PeptideSpectrumOverview>> proteinToPsms = new HashMap<>();

        Map<String, List<uk.ac.ebi.pride.archive.dataprovider.common.Triple<String, Double,String>>> proteinsPSMsScores = new HashMap<>();
        Map<String, List<String>> peptideToProteins = new HashMap<>();
        Map<String, Set<String>> proteinPTMs = new HashMap<>();
        Map<String, List<Boolean>> proteinDecoys = new HashMap<>();

        clusteredFilter.values().forEach(psmOptional -> {
            Triple<String, PeptidoformClustered, Double> psm = psmOptional.get();
            try {
                BufferedWriter batchBufferWriter = null;
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
                BackupUtil.write(archivePSM, (BufferedWriter) assayObjects.get("archiveSpectrumBufferedWriter"));
                // Total number of spectrum in Elastic Search summary.
                BackupUtil.write(psmElastic, (BufferedWriter) assayObjects.get("psmSummaryEvidenceBufferedWriter"));
                // Writing in batches.

                String usi = psm.getFirst();
                String batchFile = usi.split(":")[2];
                if(!spectraPartitionWriters.containsKey(batchFile)){
                        String prefix = (String) assayObjects.get("archiveSpectrumFilePrefix");
                        batchBufferWriter = new BufferedWriter(new FileWriter(BackupUtil.getArchiveSpectrumFileBatch(prefix, batchFile), false));
                        spectraPartitionWriters.put(batchFile, batchBufferWriter);
                }else
                    batchBufferWriter = spectraPartitionWriters.get(batchFile);

                BackupUtil.write(archivePSM, batchBufferWriter);
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
            } catch (Exception e) {
                log.debug("Error writing the PSMs in the files -- " + psm.getFirst());
                throw new RuntimeException(e);
            }
        });

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
            try {
                x.flush();
                x.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
        PrideAnalysisAssayService.proteinIndexStep(hashAssay, assayObjects, projectAccession, reanalysisAccession);
        System.out.println(clusterScores.size());

        if(assayObjects != null){
            for(Object object: assayObjects.values()){
                if (object instanceof BufferedWriter){
                    BufferedWriter bufferedWriter = (BufferedWriter) object;
                    bufferedWriter.flush();
                    bufferedWriter.close();
                }
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

    private static String makePeptideIsobaric(String peptide){
        return peptide.replace("L", "I");
    }


    private class PeptidoformClustered {
        private final String sequence;
        String peptidoform;
        Boolean isDecoy;

        public PeptidoformClustered(String sequence, String peptidoform, boolean isDecoy) {
            this.sequence = sequence;
            this.peptidoform = peptidoform;
            this.isDecoy = isDecoy;
        }

        public Boolean getDecoy() {
            return isDecoy;
        }

        public String getSequence() {
            return sequence;
        }

        @Override
        public boolean equals(Object o) {
            String oPep = makePeptideIsobaric(((PeptidoformClustered) o).peptidoform);
            String isopep = makePeptideIsobaric(peptidoform);
            if (isopep.equals(oPep)) return true;
            PeptidoformClustered that = (PeptidoformClustered) o;
            return peptidoform.equals(that.peptidoform);
        }

        @Override
        public int hashCode() {
            String isopep = makePeptideIsobaric(peptidoform);
            return Objects.hash(isopep);
        }

        @Override
        public String toString() {
            return peptidoform;
        }
    }

    private Map<String, Object> createBackupFiles(Map<String, Object> assayObjects, String folderOutput, String projectAccession, String assayAccession) throws IOException {

        // Create first the root folder for the project
        createBackupDir(folderOutput, projectAccession);

        log.info("Creating assay file  -- " + projectAccession);

        final String proteinEvidenceFileName = BackupUtil.getProteinEvidenceFile(folderOutput, projectAccession, assayAccession);
        assayObjects.put("proteinEvidenceFileName", proteinEvidenceFileName);
        assayObjects.put("proteinEvidenceBufferedWriter", new BufferedWriter(new FileWriter(proteinEvidenceFileName, false)));

        final String archiveSpectrumFileName = BackupUtil.getArchiveSpectrumFile(folderOutput, projectAccession, assayAccession);
        assayObjects.put("archiveSpectrumFileName", archiveSpectrumFileName);
        assayObjects.put("archiveSpectrumBufferedWriter", new BufferedWriter(new FileWriter(archiveSpectrumFileName, false)));

        final String archiveSpectrumFilePrefix = BackupUtil.getArchiveSpectrumFilePrefix(folderOutput, projectAccession);
        assayObjects.put("archiveSpectrumFilePrefix", archiveSpectrumFilePrefix);

        final String psmSummaryEvidenceFileName = BackupUtil.getPsmSummaryEvidenceFile(folderOutput, projectAccession, assayAccession);
        assayObjects.put("psmSummaryEvidenceFileName", psmSummaryEvidenceFileName);
        assayObjects.put("psmSummaryEvidenceBufferedWriter", new BufferedWriter(new FileWriter(psmSummaryEvidenceFileName, false)));

        return assayObjects;

    }

}
