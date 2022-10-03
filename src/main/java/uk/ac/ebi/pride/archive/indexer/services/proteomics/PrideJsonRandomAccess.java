package uk.ac.ebi.pride.archive.indexer.services.proteomics;

import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.extern.slf4j.Slf4j;
import org.ehcache.Cache;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.ArchiveSpectrum;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.BinaryArchiveSpectrum;
import uk.ac.ebi.pride.archive.indexer.utility.AppCacheManager;
import uk.ac.ebi.pride.tools.braf.BufferedRandomAccessFile;

import java.io.IOException;
import java.util.Iterator;

/**
 * {@link PrideJsonRandomAccess} is a reader of the
 * {@link uk.ac.ebi.pride.archive.dataprovider.data.spectra.SummaryArchiveSpectrum}
 * json files.
 *
 * @author ypriverol
 */
@Slf4j
public class PrideJsonRandomAccess {

    private final BufferedRandomAccessFile raf;

    private final Cache<String, Long> index;
    private final ObjectMapper objectMapper;

    public PrideJsonRandomAccess(String fileAbsolutePath) throws IOException {
        this.raf = new BufferedRandomAccessFile(fileAbsolutePath, "r", 1024 * 100);
        AppCacheManager appCacheManager = AppCacheManager.getInstance();
        this.index = appCacheManager.getPrideJsonSpectra();
        this.objectMapper = new ObjectMapper();
    }

    /**
     * Create an index of all the spectra within an {@link ArchiveSpectrum} json file.
     * The index will be a Map with usis (identifiers) as index and values the pointer
     * where the specific spectrum starts. The index is assuming that the ArchiveSpectrum
     * file is one line json representation.
     * @throws IOException
     */
    public void parseIndex() throws IOException {

        String line;
        long pos = raf.getFilePointer();

        while( (line = raf.getNextLine()) != null){
            try {
                BinaryArchiveSpectrum spectrum = BinaryArchiveSpectrum.readJson(line);
                index.put(spectrum.getUsi(), pos);
            }catch (Exception e){
                log.error("Error reading line --- " + line);
            }
            pos = raf.getFilePointer();
        }
    }

    /**
     * This class returns a specific {@link ArchiveSpectrum} for a given usi by
     * doing a lookup in the index table.
     * @param usi identifier of the spectrum
     * @return ArchiveSpectrum
     * @throws IOException
     */
    public BinaryArchiveSpectrum readArchiveSpectrum(String usi) throws IOException {
        if(index.containsKey(usi)){
            long pos = index.get(usi);
            raf.seek(pos);
            try {
                return  BinaryArchiveSpectrum.readJson(raf.readLine());
            }catch (Exception e){
                log.error("Error reading usi --- " + usi);
            }
        }
        return null;
    }

    /**
     * Return all the spectra from the Json file
     *
     * @return List of usis
     */
    public Iterator<Cache.Entry<String, Long>> getKeys(){
        return index.iterator();
    }

    public void close() throws IOException {
        raf.close();
    }
}
