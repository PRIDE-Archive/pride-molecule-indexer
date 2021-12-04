package uk.ac.ebi.pride.archive.indexer.utility;

import org.apache.commons.codec.digest.DigestUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

public class HashUtils {

    public static String getSha1Checksum(File file) throws IOException {
        FileInputStream fileInputStream = new FileInputStream(file);
        String s = DigestUtils.sha1Hex(fileInputStream);
        fileInputStream.close();
        return s;
    }

    public static String getSha256Checksum(String str) {
        return DigestUtils.sha256Hex(str);
    }
}