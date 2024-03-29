package uk.ac.ebi.pride.archive.indexer.utility;

import com.fasterxml.jackson.databind.JavaType;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.module.paranamer.ParanamerModule;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.ArchiveSpectrum;
import uk.ac.ebi.pride.archive.dataprovider.data.protein.ArchiveProteinEvidence;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.BinaryArchiveSpectrum;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.SummaryArchiveSpectrum;

import java.io.*;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class BackupUtil {

    private static final ObjectMapper objectMapper;
    public static final String JSON_EXT = ".json";
    public static final String BJSON_EXT = ".bjson";

    static {
        objectMapper = new ObjectMapper();
        objectMapper.registerModule(new ParanamerModule());
    }

    public static void writeBinarySpectrum(BinaryArchiveSpectrum spec, PrintWriter bw, boolean flush) throws IOException {
        // This step has been added to do not export orginally spectra that can be binary
        try{
            String binaryLine = BinaryArchiveSpectrum.writeJson(spec);
            BinaryArchiveSpectrum binSpectrumCheck = BinaryArchiveSpectrum.readJson(binaryLine);
            if(binSpectrumCheck == null || binSpectrumCheck.getMasses().length == 0 || binSpectrumCheck.getIntensities().length == 0)
                throw new Exception(String.format("Error with spectrum -- %s", spec.toString()));
        }catch (Exception e){
            System.err.printf("Spectrum -- %s has an error%n", spec.toString());
        }
        bw.write(BinaryArchiveSpectrum.writeJson(spec));
        bw.println();
        if(flush) bw.flush();
    }

    public static void write(Object obj, PrintWriter bw, boolean flush) throws Exception {
        String s = objectMapper.writeValueAsString(obj);
        bw.write(s);
        bw.println();
        if(flush) bw.flush();
    }

    public static String getProteinEvidenceFile(String backupPath, String projectAccession, String assayAccession) {
        if (!backupPath.endsWith(File.separator)) {
            backupPath = backupPath + File.separator;
        }
        return backupPath + projectAccession + File.separator + projectAccession + "_" + assayAccession + "_" + ArchiveProteinEvidence.class.getSimpleName() + JSON_EXT;
    }

    public static String getArchiveSpectrumFileBatch(String filePrefix, String batch){
        return filePrefix + "_" + batch + "_" + ArchiveSpectrum.class.getSimpleName() + JSON_EXT;
    }

    public static String getArchiveSpectrumFile(String backupPath, String projectAccession, String assayAccession) {
        if (!backupPath.endsWith(File.separator)) {
            backupPath = backupPath + File.separator;
        }
        return backupPath + projectAccession + File.separator + projectAccession + "_" + assayAccession + "_" +ArchiveSpectrum.class.getSimpleName() + "_Total" + JSON_EXT;
    }

    public static String getArchiveSpectrumFilePrefix(String backupPath, String projectAccession) {
        if (!backupPath.endsWith(File.separator)) {
            backupPath = backupPath + File.separator;
        }
        return backupPath + projectAccession + File.separator + projectAccession;
    }

    public static String getPsmSummaryEvidenceFile(String backupPath, String projectAccession, String assayAccession) {
        if (!backupPath.endsWith(File.separator)) {
            backupPath = backupPath + File.separator;
        }
        return backupPath + projectAccession + File.separator + projectAccession + "_" + assayAccession + "_" + SummaryArchiveSpectrum.class.getSimpleName() + JSON_EXT;
    }

    public static <T> List<T> getObjectsFromFile(Path file, Class classType) throws Exception {
        List<T> list = new ArrayList<>();
        JavaType javaType = objectMapper.getTypeFactory().constructType(classType);
        BufferedReader reader;
        reader = new BufferedReader(new FileReader(file.toFile()));
        String line = reader.readLine();
        while (line != null) {
            list.add(objectMapper.readValue(line, javaType));
            line = reader.readLine();
        }
        reader.close();

        return list;
    }
}
