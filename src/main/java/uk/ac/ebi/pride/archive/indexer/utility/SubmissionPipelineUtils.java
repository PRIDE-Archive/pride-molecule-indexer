package uk.ac.ebi.pride.archive.indexer.utility;

import de.mpc.PD.Spectra;
import de.mpc.pia.intermediate.Modification;
import de.mpc.pia.intermediate.PeptideSpectrumMatch;
import de.mpc.pia.modeller.psm.ReportPSM;
import org.apache.commons.io.FilenameUtils;
import uk.ac.ebi.jmzidml.model.mzidml.FileFormat;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.pride.archive.spectra.utils.Constants;
import uk.ac.ebi.pride.utilities.util.Triple;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

/**
 * This class contains a set of constants that are needed to process the data in the submission pipeline.
 *
 * @author ypriverol
 */
public class SubmissionPipelineUtils {

    /**
     * This method returns true if the filename is extension is a compression file .gz or .zip
     * @param resultFile Result file
     * @return Boolean
     */
    public static boolean isCompressedByExtension(String resultFile) {
        return resultFile.endsWith(Compress_Type.GZIP.extension) || resultFile.endsWith(Compress_Type.ZIP.extension);
    }

    /**
     * Supported id format used in the spectrum file.
     */
    public enum SpecIdFormat {
        MASCOT_QUERY_NUM,
        MULTI_PEAK_LIST_NATIVE_ID,
        SINGLE_PEAK_LIST_NATIVE_ID,
        SCAN_NUMBER_NATIVE_ID,
        MZML_ID,
        MZDATA_ID,
        WIFF_NATIVE_ID,
        SPECTRUM_NATIVE_ID,
        WIFF_MGF_TITLE,
        NONE
    }

    private static final String SIGN = "[+-]";
    public static final String INTEGER = SIGN + "?\\d+";


    public enum FileType {
        PRIDE,
        MZTAB,
        MZID,
        MGF,
        MS2,
        MZML,
        MZXML,
        DTA,
        PKL,
        APL;

        /**
         * This method returns the {@link FileType} of a result file using the name of the file
         * @param filename file name
         * @return FileType
         */
        public static FileType getFileTypeFromFileName(String filename) {
            filename = returnUnCompressPath(filename.toLowerCase());
            if (filename.toLowerCase().endsWith("mzid") || filename.toLowerCase().endsWith("mzidentml")) {
                return MZID;
            } else if (filename.toLowerCase().endsWith("mzml")) {
                return MZML;
            } else if (filename.toLowerCase().endsWith("mgf")) {
                return MGF;
            } else if (filename.toLowerCase().endsWith("mzxml")) {
                return MZXML;
            } else if (filename.toLowerCase().endsWith("mztab")) {
                return MZTAB;
            } else if (filename.toLowerCase().endsWith("apl")) {
                return APL;
            } else if (filename.toLowerCase().endsWith(".xml"))
                return PRIDE;

            return null;
        }

        public static FileType getFileTypeFromSpectraData(SpectraData spectraData) {
            FileFormat specFileFormat = spectraData.getFileFormat();
            FileType fileType = null;
            if (specFileFormat != null) {
                if (specFileFormat.getCvParam().getAccession().equals("MS:1000613")) fileType = DTA;
                if (specFileFormat.getCvParam().getAccession().equals("MS:1001062")) fileType = MGF;
                if (specFileFormat.getCvParam().getAccession().equals("MS:1000565")) fileType = PKL;
                if (specFileFormat.getCvParam().getAccession().equals("MS:1002996")) fileType = APL;
                if (specFileFormat.getCvParam().getAccession().equals("MS:1000584") || specFileFormat.getCvParam().getAccession().equals("MS:1000562"))
                    fileType = MZML;
                if (specFileFormat.getCvParam().getAccession().equals("MS:1000566")) fileType = MZXML;
                if (specFileFormat.getCvParam().getAccession().equals("MS:1001466")) fileType = MS2;
                if (specFileFormat.getCvParam().getAccession().equals("MS:1002600")) fileType = PRIDE;
            }else{
                fileType = getFileTypeFromFileName(spectraData.getLocation());
            }
            return fileType;
        }
    }

    public enum Compress_Type {
        GZIP("gz"),
        ZIP("zip");

        String extension;

        Compress_Type(String extension) {
            this.extension = extension;
        }

        public String getExtension() {
            return extension;
        }
    }

    public static String buildInternalPath(String productionPath, String projectAccession, String publicationYear, String publicationMonth) {
        return productionPath + publicationYear + "/" + publicationMonth + "/" + projectAccession + "/" + "internal/";
    }

    /**
     * This method remove from the file name the compression extension
     * @param originalPath File path
     * @return return a file name without the compression extension.
     */
    public static String returnUnCompressPath(String originalPath) {
        if (originalPath.endsWith(Compress_Type.GZIP.extension) || originalPath.endsWith(Compress_Type.ZIP.extension)) {
            return originalPath.substring(0, originalPath.length() - 3);
        }
        return originalPath;
    }

    /**
     * Check if the ms File is supported and match with some of the par of the name in the Spectra Files
     * This method should be used in high-throughput, when you add different files.
     *
     * @param msIdentMLFiles List of  the MS files related with the MZIdentML
     * @return The relation between the SpectraData and the corresponding File.
     */
    public static List<Triple<String, SpectraData, FileType>> combineSpectraControllers(String buildPath, List<String> msIdentMLFiles, List<SpectraData> spectraDataList) {

        List<Triple<String, SpectraData, FileType>> spectraFileMap = new ArrayList<>();

        for (String file : msIdentMLFiles) {
            for (SpectraData spectraData : spectraDataList) {
                if (spectraData.getLocation() != null && spectraData.getLocation().toLowerCase().contains(file.toLowerCase())) {
                    spectraFileMap.add(new Triple<>(buildPath + file, spectraData,
                            FileType.getFileTypeFromSpectraData(spectraData)));
                } else if (file.contains(spectraData.getId())
                        || (spectraData.getName() != null && file.toLowerCase().contains(spectraData.getName().toLowerCase()))) {
                    spectraFileMap.add(new Triple<>(buildPath + file, spectraData, FileType.getFileTypeFromSpectraData(spectraData)));
                }
            }
        }
        return spectraFileMap;
    }

    /**
     * The following method get for one spectrum, Triple<String, String, FileType>. The elements of the Triple are:
     *  - First: FileName provided by the user
     *  - Second: Accession of the spectrum in the file
     *  - Third: FileType
     * @param spectraDataList All the files
     * @param psm PSM
     * @return Triple
     * @throws IOException
     */
    public static Triple<String, String, SubmissionPipelineUtils.FileType> getSpectrumId(List<Triple<String, SpectraData, FileType>> spectraDataList, PeptideSpectrumMatch psm) throws IOException {

        Optional<Triple<String, SpectraData, FileType>> spectraDataOption = spectraDataList.stream().filter(x -> psm.getSpectrumIdentification()
                .getInputSpectra()
                .stream()
                .anyMatch(y -> y.getSpectraDataRef().equalsIgnoreCase(x.getSecond().getId())))
                .findAny();
        if(!spectraDataOption.isPresent())
            throw new IOException(String.format("The psm reference can't be found in the list of files -- %s", psm.toString()));

        Triple<String, SpectraData, FileType> spectraData = spectraDataOption.get();
        SpecIdFormat fileIdFormat = getSpectraDataIdFormat(spectraData.getSecond().getSpectrumIDFormat().getCvParam().getAccession());

        if (fileIdFormat == SpecIdFormat.MASCOT_QUERY_NUM || fileIdFormat == SpecIdFormat.MULTI_PEAK_LIST_NATIVE_ID) {
            String rValueStr = psm.getSourceID().replaceAll("query=", "").replaceAll("index=", "");
            String id = null;
            if (rValueStr.matches(INTEGER)) {
                id = Integer.toString(Integer.parseInt(rValueStr) + 1);
            }
            return new Triple<String, String, FileType>(spectraData.getFirst(), id, spectraData.getThird());
//        } else if (fileIdFormat == SpecIdFormat.SINGLE_PEAK_LIST_NATIVE_ID) {
//            return psm.getSourceID().replaceAll("file=", "");
//        } else if (fileIdFormat == SpecIdFormat.MZML_ID) {
//            return psm.getSourceID().replaceAll("mzMLid=", "");
//        } else if (fileIdFormat == SpecIdFormat.SCAN_NUMBER_NATIVE_ID) {
//            return psm.getSourceID().replaceAll("scan=", "");
        } else if (fileIdFormat == SpecIdFormat.SPECTRUM_NATIVE_ID) {
            String[] partsScan = psm.getSourceID().split(" ");
            Optional<String> scanOptional = Arrays.stream(partsScan).filter(x -> x.contains("scan=")).findFirst();
            String id = scanOptional.map(s -> s.replaceAll("scan=", "")).orElseGet(psm::getSourceID);
            return new Triple<String, String, FileType>(spectraData.getFirst(), id, spectraData.getThird());
        } else {
            return new Triple<>(spectraData.getFirst(), psm.getSourceID(),spectraData.getThird());
        }
    }

//    public static String buildUsi(String projectAccession, Triple<String, SpectraData, FileType> refeFile, ReportPSM psm) {
//        Constants.ScanType scanType = Constants.ScanType.INDEX;
//        SpecIdFormat fileIFormat = getSpectraDataIdFormat(refeFile.getSecond().getSpectrumIDFormat().getCvParam().getAccession());
//        String spectrumID = getSpectrumId(refeFile.getSecond(), psm);
//        if (fileIFormat == SpecIdFormat.MASCOT_QUERY_NUM || fileIFormat == SpecIdFormat.MULTI_PEAK_LIST_NATIVE_ID) {
//            scanType = Constants.ScanType.INDEX;
//        } else if (fileIFormat == SpecIdFormat.MZML_ID) {
//            scanType = Constants.ScanType.SCAN;
//            String[] scanStrings = spectrumID.split("scan=");
//            spectrumID = scanStrings[1];
//        } else if(fileIFormat == SpecIdFormat.SPECTRUM_NATIVE_ID || fileIFormat == SpecIdFormat.SCAN_NUMBER_NATIVE_ID) {
//            scanType = Constants.ScanType.SCAN;
//        }
//        Path p = Paths.get(refeFile.getFirst());
//        String fileName = p.getFileName().toString();
//        return Constants.SPECTRUM_S3_HEADER + projectAccession + ":" + fileName + ":" + scanType.getName() + ":" + spectrumID + ":" + encodePSM(psm.getSequence(), psm.getModifications(), psm.getCharge());
//    }

    public static String getSpectraUsiFromUsi(String usi){
        String spectraUsi;
        String[] usiArray = usi.split(":");
        String[] subset = Arrays.copyOfRange(usiArray, 0, 5);
        subset[2] = getFileNameNoExtension(subset[2]);
        spectraUsi = String.join(":", subset);
        return spectraUsi;
    }

    /**
     * Get the name of a file removing the following components:
     * - extension gz or zip
     * - original extension of the file .xml or .mzML
     * - parent path of the file (e.g. /home/user/)
     * @param fullPath full path
     * @return name of the file
     */
    public static String getFileNameNoExtension(String fullPath){
        String fileName = returnUnCompressPath(FilenameUtils.getName(fullPath));
        return fileName.substring(0, fileName.lastIndexOf('.'));
    }

    /**
     * build USI for PRIDE XML as spectra
     *
     * @param projectAccession project
     * @param fileName         filename with the spectra
     * @param psm              PSM
     * @return
     */
    public static String buildUsi(String projectAccession, String fileName, ReportPSM psm, String id, FileType scan) {
        Constants.ScanType scanType = Constants.ScanType.INDEX;
        if(scan == FileType.MZML)
            scanType = Constants.ScanType.SCAN;
        return Constants.SPECTRUM_S3_HEADER + projectAccession + ":" + fileName + ":" + scanType.getName() + ":" + id + ":" + encodePSM(psm.getSequence(), psm.getModifications(), psm.getCharge());
    }

    public static String encodePSM(String sequence, Map<Integer, Modification> ptms, Integer charge) {
        return encodePeptide(sequence, ptms) + "/" + charge;
    }

    public static String encodePeptide(String sequence, Map<Integer, Modification> ptms) {
        StringBuilder stringBuilder = new StringBuilder();
        String finalSequence = sequence;
        if (ptms != null && ptms.size() > 0) {
            char[] sequenceList = sequence.toCharArray();
            if (ptms.containsKey(0))
                stringBuilder.append("[").append(ptms.get(0).getAccession()).append("]");
            for (int i = 0; i < sequenceList.length; i++) {
                stringBuilder.append(sequenceList[i]);
                if (ptms.containsKey(i + 1)) {
                    stringBuilder.append("[").append(ptms.get(i + 1).getAccession()).append("]");
                }
            }

            // Add the CTerm modifications
            for (Map.Entry entry : ptms.entrySet()) {
                Integer position = (Integer) entry.getKey();
                Modification mod = (Modification) entry.getValue();
                if (position > sequence.length()) {
                    stringBuilder.append("-").append("[").append(mod.getAccession()).append("]");
                }
            }
            finalSequence = stringBuilder.toString();
        }
        return finalSequence;
    }


    /**
     * Spectrum Id format for an specific CVterm accession
     *
     * @param accession CvTerm Accession
     * @return Specific Spectrum Id Format
     */
    public static SpecIdFormat getSpectraDataIdFormat(String accession) {
        if (accession.equals("MS:1001528")) return SpecIdFormat.MASCOT_QUERY_NUM;
        if (accession.equals("MS:1000774")) return SpecIdFormat.MULTI_PEAK_LIST_NATIVE_ID;
        if (accession.equals("MS:1000775")) return SpecIdFormat.SINGLE_PEAK_LIST_NATIVE_ID;
        if (accession.equals("MS:1001530")) return SpecIdFormat.MZML_ID;
        if (accession.equals("MS:1000776")) return SpecIdFormat.SCAN_NUMBER_NATIVE_ID;
        if (accession.equals("MS:1000770")) return SpecIdFormat.WIFF_NATIVE_ID;
        if (accession.equals("MS:1000777")) return SpecIdFormat.MZDATA_ID;
        if (accession.equals(("MS:1000768"))) return SpecIdFormat.SPECTRUM_NATIVE_ID;
        if (accession.equals("MS:1000796")) return SpecIdFormat.WIFF_MGF_TITLE;
        return SpecIdFormat.NONE;
    }

    /**
     * This function aim to translate a Qvalue of 0.0 to the minor Qvalue from the distribution of qvalues.
     * @param currentQValue Current Q-value
     * @param allQValues All q-values
     * @return new q-values.
     */
    public static double getQValueLower(double currentQValue, Set<Double> allQValues){
        if(currentQValue > 0.0)
            return currentQValue;
        Double minQValue = Collections.min(allQValues.stream().filter( x -> x > 0.0).collect(Collectors.toList()));
        return new BigDecimal(minQValue/10).setScale(6, RoundingMode.HALF_UP).doubleValue();
    }

    /**
     * This method copy a gzip file from the submitted folder to uncompress internal file
     * @param compressFile Compress file mzml.gz
     * @param uncompresFile Uncompress file
     * @param productionPath Production path
     */
    public static void copyUncompressSpectraFile(String compressFile, String uncompresFile, String productionPath){

        String sourcePath = String.join("/", productionPath, "submitted", compressFile);
        String targetPath = String.join("/", productionPath,  "internal", uncompresFile);
        try (GZIPInputStream gis = new GZIPInputStream(
                new FileInputStream(sourcePath))) {
            Files.copy(gis, new File(targetPath).toPath());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
