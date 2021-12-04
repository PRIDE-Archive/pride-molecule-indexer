package uk.ac.ebi.pride.archive.indexer.utility;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

public class DateUtils {
   private static final DateFormat DATE_FORMAT_DATE_PART = new SimpleDateFormat("yyyy-MM-dd");

    public static boolean equalsDatePartOnly(Date a, Date b) {
        String aStr = DATE_FORMAT_DATE_PART.format(a);
        String bStr = DATE_FORMAT_DATE_PART.format(b);
        return  ((a == b) || (a != null && b != null && aStr.equals(bStr)));
    }
}
