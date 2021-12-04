package uk.ac.ebi.pride.archive.indexer.utility;

/**
 * This code is licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * ==Overview==
 * <p>
 * This class
 * <p>
 * @author ypriverol
 */
public class StringUtils {


    public static final String EMPTY_STRING = "";
    public static final String URL_SEPARATOR = "/";

    // General Sample processing Capture
    public static final String CV_LABEL_ORGANISM = "NEWT";
    public static final String CV_LABEL_CELL_COMPONENT = "CL";
    public static final String CV_LABEL_CELL_TISSUE = "BTO";
    public static final String CV_LABEL_DISEASE = "DOID";

    /**
     * Get convert sentence to Capitalize Style
     * @param sentence original sentence
     * @return Capitalize sentence
     */
    public static String convertSentenceStyle(String sentence){
        sentence = sentence.toLowerCase().trim();
        return org.apache.commons.lang3.StringUtils.capitalize(sentence);
    }

}
