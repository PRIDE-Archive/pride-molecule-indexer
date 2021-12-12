package uk.ac.ebi.pride.archive.indexer.services.proteomics;

import lombok.extern.slf4j.Slf4j;
import uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzdata_wrapper.MzMlWrapper;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;
import uk.ac.ebi.pride.tools.pkl_parser.PklFile;
import uk.ac.ebi.pride.tools.pride_wrapper.PRIDEXmlWrapper;
import uk.ac.ebi.pride.utilities.util.Tuple;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@Slf4j
public class JmzReaderSpectrumService {

    // Map of all readers containing the spectra
    final Map<String, JMzReader> readers;

    /**
     * Create a service from a list of files with the corresponding file types {@link uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils.FileType}
     * @param spectrumFileList List of Tuple of filePath and corresponding file type
     * @throws JMzReaderException
     * @throws MzXMLParsingException
     */
    private JmzReaderSpectrumService(List<Tuple<String, SubmissionPipelineUtils.FileType>> spectrumFileList) throws JMzReaderException, MzXMLParsingException {
        this.readers = new HashMap<>();
        for (Tuple<String, SubmissionPipelineUtils.FileType> entry : spectrumFileList) {
            String key = entry.getKey();
            SubmissionPipelineUtils.FileType value = entry.getValue();

            if (value == SubmissionPipelineUtils.FileType.MGF) {
                this.readers.put(key, new MgfFile(new File(key), true));
            }
            if (value == SubmissionPipelineUtils.FileType.PRIDE) {
                this.readers.put(key, new PRIDEXmlWrapper(new File(key)));
            }
            if( value == SubmissionPipelineUtils.FileType.MZML){
                this.readers.put(key, new MzMlWrapper(new File(key)));
            }
            if( value == SubmissionPipelineUtils.FileType.PKL){
                this.readers.put(key, new PklFile(new File(key)));
            }
            if( value == SubmissionPipelineUtils.FileType.MZXML){
                this.readers.put(key, new MzXMLFile(new File(key)));
            }
        }
    }

    /**
     * Return an instance that allow to read the spectra from the original file.
     *
     * @param spectrumFileList
     * @return
     * @throws JMzReaderException
     */
    public static JmzReaderSpectrumService getInstance(List<Tuple<String, SubmissionPipelineUtils.FileType>> spectrumFileList) throws JMzReaderException, MzXMLParsingException {
        return new JmzReaderSpectrumService(spectrumFileList);
    }

    public Spectrum getSpectrumById(String filePath, String id) throws JMzReaderException {
        JMzReader reader = readers.get(filePath);
        try{
            Spectrum spec = null;
            if(reader.getSpectraIds().contains(id)){
                spec = reader.getSpectrumById(id);
            }
            /***
             * In some cases the id is written wronly and the scan is annotated as controllerType=0 controllerNumber=1 scan=1 while in the
             * mzidentml is annotated as scan=1
             */
            if (spec == null){
                List<Tuple<String, String>> spectra = reader.getSpectraIds()
                        .stream().filter(x -> x.contains(id))
                        .map(x -> new Tuple<>(x, x)).collect(Collectors.toList());
                if(spectra.size() == 1)
                    spec =  reader.getSpectrumById(spectra.get(0).getValue());
                else{
                    spectra = reader.getSpectraIds().stream().map(x-> {
                        String[] scans = x.split(" ");
                        for(String idScan: scans){
                            if(idScan.contains("scan")){
                                return new Tuple<>(idScan.trim(),x);
                            }
                        }
                        return new Tuple<>(x,x);
                    }).collect(Collectors.toList()).stream().filter(x -> x.getKey().equalsIgnoreCase("scan=" + id)).collect(Collectors.toList());
                    if(spectra.size() == 1)
                        spec =  reader.getSpectrumById(spectra.get(0).getValue());
                }
            }
            if(spec != null && spec.getMsLevel() == 1)
                spec = null;
            return spec;
        }catch (NumberFormatException e){
            throw new JMzReaderException("Error parsing the following Accession -- " + id);
        }
    }

    /**
     * Read Spectrum from indexed files such as MGF or MS2
     * @param filePath FilePath of the file
     * @param id index to be search
     * @return Spectrum found in the file
     * @throws JMzReaderException
     */
    public Spectrum getSpectrumByIndex(String filePath, String id) throws JMzReaderException {
        JMzReader reader = readers.get(filePath);
        try{
            int index = Integer.parseInt(id);
            return reader.getSpectrumByIndex(index);
        }catch (Exception e){
            throw new JMzReaderException("Error parsing the following Accession -- " + id);
        }
    }


}