package uk.ac.ebi.pride.archive.indexer.utility;

import lombok.extern.slf4j.Slf4j;


import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
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
        final DigestInputStream digestInputStream = new DigestInputStream(Files.newInputStream(new File(filepath).toPath()), inputStreamMessageDigest);
        while (digestInputStream.read(bytesRead) != -1) ;
        return normalize(inputStreamMessageDigest);
    }

    public static String sha1InObject(Object object) throws Exception {
        if (object == null) {
            throw new Exception("Object is null.");
        }

        String input = String.valueOf(object);

        MessageDigest md;
        try {
            md = MessageDigest.getInstance("SHA1");
        } catch (NoSuchAlgorithmException ex) {
            return null;
        }
        md.reset();

        byte[] buffer = input.getBytes();
        md.update(buffer);

        byte[] digest = md.digest();
        String hexStr = "";
        for (int i = 0; i < digest.length; i++) {
            hexStr += Integer.toString((digest[i] & 0xff) + 0x100, 16).substring(1);
        }
        return hexStr;
    }
}