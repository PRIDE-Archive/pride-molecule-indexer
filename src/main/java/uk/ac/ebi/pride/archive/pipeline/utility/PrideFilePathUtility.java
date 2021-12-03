package uk.ac.ebi.pride.archive.pipeline.utility;

import lombok.extern.log4j.Log4j;
import uk.ac.ebi.pride.mongodb.archive.model.projects.MongoPrideProject;

import java.io.File;
import java.util.Calendar;

@Log4j
public class PrideFilePathUtility {

    public static final String INTERNAL = "internal";
    public static final String SUBMITTED = "submitted";
    public static final String SUBMISSION_PX = "submission.px";
    public static final String GENERATED = "generated";
    private static final String README_TXT = "README.txt";


    public static String getSubmissionFilePath(MongoPrideProject prideProject, String prideRepoRootPath) {
        return getPublicPath(prideProject, prideRepoRootPath) + File.separator + INTERNAL + File.separator + SUBMISSION_PX;
    }

    public static String getSubmittedFilesPath(MongoPrideProject prideProject, String prideRepoRootPath) {
        return getPublicPath(prideProject, prideRepoRootPath) + File.separator + SUBMITTED + File.separator;
    }

    private static String getPublicPath(MongoPrideProject prideProject, String prideRepoRootPath) {
        log.info("Generating public file path fragment based on the publication date and project accession");
        Calendar calendar = Calendar.getInstance();
        calendar.setTime(prideProject.getPublicationDate());
        int month = calendar.get(Calendar.MONTH) + 1; // the month are zero based, hence the correction +1
        int year = calendar.get(Calendar.YEAR);
        String datePath = year + File.separator + (month < 10 ? "0" : "") + month;
        String publicPath = datePath + File.separator + prideProject.getAccession();
        log.info("Generated public path fragment: " + publicPath);
        prideRepoRootPath = prideRepoRootPath.endsWith(File.separator) ? prideRepoRootPath : prideRepoRootPath + File.separator;
        return prideRepoRootPath + publicPath;
    }

    public static String getReadMeFilePath(MongoPrideProject prideProject, String prideRepoRootPath) {
        return getPublicPath(prideProject, prideRepoRootPath) + File.separator + GENERATED + File.separator + README_TXT;
    }
}
