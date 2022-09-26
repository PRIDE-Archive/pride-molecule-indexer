package uk.ac.ebi.pride.archive.indexer.services;

import lombok.extern.slf4j.Slf4j;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.stereotype.Service;
import uk.ac.ebi.pride.archive.dataprovider.common.Triple;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.BinaryArchiveSpectrum;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.MGFPRIDEWriter;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PIAModelerService;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PrideJsonRandomAccess;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.stream.Collectors;

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
                Collectors.toMap(x -> x.getKey(), x -> {
                    List<Triple<String, Double, String>> list = x.getValue();
                    Collections.sort(list, comparator);
                    return list.stream()
                            .findFirst()
                            .get()
                            .getSecond();
                }));

        proteinScore.forEach((key, value) -> System.out.println(key + "  " + value));
        return proteinScore;
    }

    public void performProteinInference(String pridePSMPath, String prideProject) throws IOException {
        PrideJsonRandomAccess pridePSMJsonReader = new PrideJsonRandomAccess(pridePSMPath);
        pridePSMJsonReader.parseIndex();
        // ProteinAccession, List<Triple<Peptidoform, Score, Usi>>
        Map<String, List<Triple<String, Double,String>>> proteins = new HashMap<>();

        for(String usi: pridePSMJsonReader.getKeys()){
            BinaryArchiveSpectrum spectrum = pridePSMJsonReader.readArchiveSpectrum(usi);
            List<String> proteinAccessions = spectrum.getProteinAccessions();
            Double pcmScore = Double.parseDouble(spectrum.getBestSearchEngineScore().getValue());
            for(String protein: proteinAccessions){
                List<Triple<String, Double,String>> pcms = new ArrayList<>();
                if(proteins.containsKey(protein)){
                    pcms = proteins.get(protein);
                }
                pcms.add(new Triple<>(spectrum.getPeptidoform(), pcmScore, spectrum.getUsi()));
                proteins.put(protein, pcms);
            }
        }
        Map<String, Double> proteinScores = getBestQValue(proteins);
        log.info(String.valueOf(proteins.size()));
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
