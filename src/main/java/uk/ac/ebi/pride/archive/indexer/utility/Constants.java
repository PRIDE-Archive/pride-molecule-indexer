package uk.ac.ebi.pride.archive.indexer.utility;

public class Constants {

    public static final String SPECTRUM_S3_HEADER = "mzspec:";

    public enum ScanType{
        INDEX("index"),
        SCAN("scan"), NATIVE_ID("nativeId");

        private final String name;

        ScanType(String name){
            this.name = name;
        }

        public String getName() {
            return name;
        }
    }

    /**
     * This method build the Unified Spectrum Identifier
     * @param scanType Scan Type can be scan id or index
     * @param project ProteomeXchange project
     * @param fileName file Name
     * @param spectrumId Spectrum Id (can be scan number of index number )
     * @return
     */
    public static String buildUsi(ScanType scanType, String project, String fileName, String spectrumId){
        return SPECTRUM_S3_HEADER + project + ":" + fileName + ":" + scanType.getName() + ":" + spectrumId;
    }

}
