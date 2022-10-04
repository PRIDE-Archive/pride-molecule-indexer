package uk.ac.ebi.pride.archive.indexer.services;

import lombok.extern.slf4j.Slf4j;
import org.ehcache.Cache;
import org.springframework.stereotype.Service;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.BinaryArchiveSpectrum;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.MGFPRIDEWriter;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.PrideJsonRandomAccess;
import uk.ac.ebi.pride.archive.indexer.utility.AppCacheManager;
import uk.ac.ebi.pride.archive.indexer.utility.BackupUtil;

import java.io.*;
import java.nio.file.Files;
import java.util.Iterator;

@Slf4j
@Service
public class PSMClusteringService {

    public PSMClusteringService() {
    }

    public void convertToMgf(String prideJsonFile, String mgfOutputFile) {
        try {
            PrideJsonRandomAccess pridePSMJsonReader = new PrideJsonRandomAccess(prideJsonFile);
            pridePSMJsonReader.parseIndex();
            OutputStream outputStream       = Files.newOutputStream((new File(mgfOutputFile)).toPath());
            OutputStreamWriter outputStreamWriter = new OutputStreamWriter(outputStream);

            for (Iterator<Cache.Entry<String, Long>> it = pridePSMJsonReader.getKeys(); it.hasNext(); ) {
                String usi = it.next().getKey();
                BinaryArchiveSpectrum spec = pridePSMJsonReader.readArchiveSpectrum(usi);
                if (spec != null)
                    MGFPRIDEWriter.appendSpectrum(outputStreamWriter, spec);
            }
            outputStreamWriter.flush();
            outputStreamWriter.close();
            pridePSMJsonReader.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * MaraCluster output is a file with 3 columns, with the following structure:
     * spectra file , index of the spectra in the file (0-based), cluster index
     *
     * This method transform the file into a Map where the key is the index of the spectra in the pridejson (0-based)
     * and the cluster where it belongs.
     *
     * @param maraClusterFile MaraCluster file output
     * @return Map with the index of the spectrum and the cluster index.
     * @throws Exception
     */
    public Cache<Integer, Integer> readMaraClusterResults(String maraClusterFile) throws Exception {

        AppCacheManager cacheManager = AppCacheManager.getInstance();

        InputStream in = Files.newInputStream(new File(maraClusterFile).toPath());
        BufferedReader br = new BufferedReader(new InputStreamReader(in));

        // Clusters result will be an index of the spectrum in the pride file and the cluster where it belongs to.
        Cache<Integer, Integer> clusters = cacheManager.getClustersCache();
        String line;
        while ((line = br.readLine()) != null) {
            if(!line.trim().isEmpty()){
                String[] values = line.split("\t");
                if(values.length == 3 && !values[0].trim().isEmpty() &&
                        !values[1].trim().isEmpty() &&
                        !values[2].trim().isEmpty()){

                    if(!clusters.containsKey(Integer.parseInt(values[1].trim())))
                        clusters.put(Integer.parseInt(values[1].trim()), Integer.parseInt(values[2].trim()));
                    else
                        throw new Exception("The following spectra -- " + values[1].trim() + " belongs to more than one cluster");
                }
            }
        }
        return clusters;
    }

    public void validateJsonFile(String spectraArchiveFile, String validatedArchiveFile) {
        try {
            PrideJsonRandomAccess pridePSMJsonReader = new PrideJsonRandomAccess(spectraArchiveFile);
            pridePSMJsonReader.parseIndex();

            BufferedWriter bw = new BufferedWriter(new FileWriter(validatedArchiveFile, false));

            for (Iterator<Cache.Entry<String, Long>> it = pridePSMJsonReader.getKeys(); it.hasNext(); ) {
                String usi = it.next().getKey();
                BinaryArchiveSpectrum spec = pridePSMJsonReader.readArchiveSpectrum(usi);
                if (spec != null){
                    BackupUtil.write(spec, bw);
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
