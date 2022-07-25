package uk.ac.ebi.pride.archive.indexer.utility;

import lombok.extern.slf4j.Slf4j;


import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import javax.xml.bind.DatatypeConverter;
import java.security.NoSuchAlgorithmException;

@Slf4j
public class HashUtils {
    private static final Integer BUFFER_SIZE = 2048;

    public static MessageDigest getSha1() {
        return getHashingAlgorithm("SHA1");
    }

    private static MessageDigest getHashingAlgorithm(final String algorithm) {
        try {
            return MessageDigest.getInstance(algorithm);
        } catch (NoSuchAlgorithmException e) {
            log.error("Error in getting hash algorithm - {}", e.getMessage());
            throw new AssertionError(e);
        }
    }

    public static String normalize(MessageDigest messageDigest) {
        if ("SHA1".equals(messageDigest.getAlgorithm())) {
            return DatatypeConverter.printHexBinary(messageDigest.digest()).toLowerCase();
        }
        return new String(messageDigest.digest());
    }

    public static String calculateSha1Checksum(String filepath) throws IOException {
        byte[] bytesRead = new byte[BUFFER_SIZE];
        final MessageDigest inputStreamMessageDigest = getSha1();
        final DigestInputStream digestInputStream = new DigestInputStream(new FileInputStream(new File(filepath)), inputStreamMessageDigest);
        while (digestInputStream.read(bytesRead) != -1) ;
        return normalize(inputStreamMessageDigest);
    }
}