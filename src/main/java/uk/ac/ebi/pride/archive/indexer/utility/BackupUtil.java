package uk.ac.ebi.pride.archive.indexer.utility;

import com.fasterxml.jackson.databind.JavaType;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.module.paranamer.ParanamerModule;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.model.ArchiveSpectrum;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.model.PrideProteinEvidence;
import uk.ac.ebi.pride.archive.indexer.services.proteomics.model.PridePsmSummaryEvidence;

import java.io.*;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class BackupUtil {

    private static final ObjectMapper objectMapper;
    public static final String JSON_EXT = ".json";

    static {
        objectMapper = new ObjectMapper();
        objectMapper.registerModule(new ParanamerModule());
    }

    public static void write(Object obj, BufferedWriter bw) throws Exception {
        String s = objectMapper.writeValueAsString(obj);
        bw.write(s);
        bw.newLine();
    }

    public static String getProteinEvidenceFile(String backupPath, String projectAccession, String assayAccession) {
        if (!backupPath.endsWith(File.separator)) {
            backupPath = backupPath + File.separator;
        }
        return backupPath + projectAccession + File.separator + projectAccession + "_" + assayAccession +
                "_" + PrideProteinEvidence.class.getSimpleName() + JSON_EXT;
    }

    public static String getArchiveSpectrumFile(String backupPath, String projectAccession, String assayAccession) {
        if (!backupPath.endsWith(File.separator)) {
            backupPath = backupPath + File.separator;
        }
        return backupPath + projectAccession + File.separator + projectAccession + "_" + assayAccession +
                "_" + ArchiveSpectrum.class.getSimpleName() + JSON_EXT;
    }

    public static String getPsmSummaryEvidenceFile(String backupPath, String projectAccession, String assayAccession) {
        if (!backupPath.endsWith(File.separator)) {
            backupPath = backupPath + File.separator;
        }
        return backupPath + projectAccession + File.separator + projectAccession + "_" + assayAccession +
                "_" + PridePsmSummaryEvidence.class.getSimpleName() + JSON_EXT;
    }

    public static List<PrideProteinEvidence> getProteinEvidenceFromBackup(String backupPath, String projectAccession, String assayAccession) throws IOException {
        List<PrideProteinEvidence> prideProteinEvidences = new ArrayList<>();
        BufferedReader reader;
        reader = new BufferedReader(new FileReader(getProteinEvidenceFile(backupPath, projectAccession, assayAccession)));
        String line = reader.readLine();
        while (line != null) {
            prideProteinEvidences.add(objectMapper.readValue(line, PrideProteinEvidence.class));
            line = reader.readLine();
        }
        reader.close();

        return prideProteinEvidences;
    }

    public static List<PridePsmSummaryEvidence> getPsmSummaryEvidenceFromBackup(String backupPath, String projectAccession, String assayAccession) throws IOException {
        List<PridePsmSummaryEvidence> list = new ArrayList<>();
        BufferedReader reader;
        reader = new BufferedReader(new FileReader(getPsmSummaryEvidenceFile(backupPath, projectAccession, assayAccession)));
        String line = reader.readLine();
        while (line != null) {
            list.add(objectMapper.readValue(line, PridePsmSummaryEvidence.class));
            line = reader.readLine();
        }
        reader.close();

        return list;
    }

    public static List<ArchiveSpectrum> getArchiveSpectrumFromBackup(String backupPath, String projectAccession, String assayAccession) throws IOException {
        List<ArchiveSpectrum> list = new ArrayList<>();
        BufferedReader reader;
        reader = new BufferedReader(new FileReader(getArchiveSpectrumFile(backupPath, projectAccession, assayAccession)));
        String line = reader.readLine();
        while (line != null) {
            list.add(objectMapper.readValue(line, ArchiveSpectrum.class));
            line = reader.readLine();
        }
        reader.close();

        return list;
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
