package uk.ac.ebi.pride.archive.pipeline.services.proteomics;

import com.amazonaws.services.s3.AmazonS3;
import lombok.extern.slf4j.Slf4j;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.pride.archive.pipeline.utility.SubmissionPipelineConstants;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzdata_wrapper.MzMlWrapper;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;
import uk.ac.ebi.pride.tools.pkl_parser.PklFile;
import uk.ac.ebi.pride.tools.pride_wrapper.PRIDEXmlWrapper;
import uk.ac.ebi.pride.utilities.util.Triple;
import uk.ac.ebi.pride.utilities.util.Tuple;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@Slf4j
public class JmzReaderSpectrumService {

    /**
     * Map of all readers containing the spectra
     */
    Map<String, JMzReader> readers = new HashMap<>();

    String clientRegion = "*** Client region ***";
    String bucketName = "*** Bucket name ***";
    String key = "*** Object key ***";

    AmazonS3 s3Client;

    private JmzReaderSpectrumService(List<Triple<String, SpectraData, SubmissionPipelineConstants.FileType>> spectrumFileList) throws JMzReaderException, MzXMLParsingException {
        this.readers = new HashMap<>();
        for (Triple<String, SpectraData, SubmissionPipelineConstants.FileType> entry : spectrumFileList) {
            String key = (String) entry.getFirst();
            SubmissionPipelineConstants.FileType value = entry.getThird();
//            if(value == null && entry.getSecond().getSpectrumIDFormat().getCvParam().getAccession().equals("MS:1000774")){
//                entry.setThird(SubmissionPipelineConstants.FileType.MGF);
//                value = SubmissionPipelineConstants.FileType.MGF;
//            }
            if (value == SubmissionPipelineConstants.FileType.MGF) {
                this.readers.put(key, new MgfFile(new File(key), true));
            }
            if (value == SubmissionPipelineConstants.FileType.PRIDE) {
                this.readers.put(key, new PRIDEXmlWrapper(new File(key)));
            }
            if( value == SubmissionPipelineConstants.FileType.MZML){
                this.readers.put(key, new MzMlWrapper(new File(key)));
            }
            if( value == SubmissionPipelineConstants.FileType.PKL){
                this.readers.put(key, new PklFile(new File(key)));
            }
            if( value == SubmissionPipelineConstants.FileType.MZXML){
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
    public static JmzReaderSpectrumService getInstance(List<Triple<String, SpectraData, SubmissionPipelineConstants.FileType>> spectrumFileList) throws JMzReaderException, MzXMLParsingException {
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
                        .map(x -> new Tuple<String, String>(x,x)).collect(Collectors.toList());
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

    public Spectrum getSpectrumByIndex(String filePath, int id) throws JMzReaderException {
        JMzReader reader = readers.get(filePath);
        try{
            return reader.getSpectrumByIndex(id);
        }catch (NumberFormatException e){
            throw new JMzReaderException("Error parsing the following Accession -- " + id);
        }
    }


}