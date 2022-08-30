package uk.ac.ebi.pride.archive.indexer.services.proteomics;

import de.mpc.pia.intermediate.compiler.PIACompiler;
import de.mpc.pia.intermediate.compiler.PIASimpleCompiler;
import de.mpc.pia.intermediate.compiler.parser.InputFileParserFactory;
import de.mpc.pia.modeller.PIAModeller;
import de.mpc.pia.modeller.protein.inference.OccamsRazorInference;
import de.mpc.pia.modeller.protein.scoring.AbstractScoring;
import de.mpc.pia.modeller.protein.scoring.MultiplicativeScoring;
import de.mpc.pia.modeller.protein.scoring.settings.PSMForScoring;
import de.mpc.pia.modeller.report.filter.FilterComparator;
import de.mpc.pia.modeller.report.filter.RegisteredFilters;
import de.mpc.pia.modeller.report.filter.impl.PSMScoreFilter;
import de.mpc.pia.modeller.score.FDRData;
import de.mpc.pia.modeller.score.ScoreModelEnum;
import lombok.extern.slf4j.Slf4j;
import uk.ac.ebi.pride.archive.indexer.utility.SubmissionPipelineUtils;


import java.io.File;
import java.io.IOException;
import java.util.HashMap;

@Slf4j
public class PIAModelerService {

    public PIAModelerService() {
    }

    /**
     * Perform the protein inference for a file. Including all the thershold.
     * @param filePath assay file path, pride xml, mzidentml
     * @param psmQThreshold q-value threshold
     * @param proteinQThreshold q-value threshold
     */
    public PIAModeller performProteinInference(String assayId,
                                               String filePath, SubmissionPipelineUtils.FileType fileType,
                                               double psmQThreshold, double proteinQThreshold) throws IOException {

        PIAModeller piaModeller = computeFDRPSMLevel(assayId, filePath, fileType);
        return performFilteringInference(piaModeller, psmQThreshold, proteinQThreshold);


    }

    public PIAModeller performFilteringInference(PIAModeller piaModeller, double psmQThreshold, double proteinQThreshold){

        if (piaModeller != null){

            piaModeller.setCreatePSMSets(false);

            piaModeller.getPSMModeller().setAllDecoyPattern("searchengine");
            piaModeller.getPSMModeller().setAllTopIdentifications(0);

            piaModeller.getPSMModeller().addFilter(0L, new PSMScoreFilter(FilterComparator.less_equal,
                    false, psmQThreshold, ScoreModelEnum.PSM_LEVEL_FDR_SCORE.getShortName()));

            piaModeller.getPSMModeller().addFilter(0L, RegisteredFilters.PSM_SOURCE_ID_FILTER
                    .newInstanceOf(FilterComparator.equal,"index=null" ,true));

            piaModeller.getPSMModeller().calculateAllFDR();
            piaModeller.getPSMModeller().calculateCombinedFDRScore();
            piaModeller.setConsiderModifications(true);

            // protein level
            OccamsRazorInference seInference = new OccamsRazorInference();

            seInference.addFilter(new PSMScoreFilter(FilterComparator.less_equal,
                    false, psmQThreshold, ScoreModelEnum.PSM_LEVEL_FDR_SCORE.getShortName()));
            seInference.addFilter(RegisteredFilters.PSM_SOURCE_ID_FILTER
                    .newInstanceOf(FilterComparator.equal,"index=null" ,true));

            seInference.setScoring(new MultiplicativeScoring(new HashMap<>()));
            seInference.getScoring()
                    .setSetting(AbstractScoring.SCORING_SETTING_ID,
                            ScoreModelEnum.PSM_LEVEL_FDR_SCORE.getShortName());
            seInference.getScoring()
                    .setSetting(AbstractScoring.SCORING_SPECTRA_SETTING_ID,
                            PSMForScoring.ONLY_BEST.getShortName());


            piaModeller.getProteinModeller().infereProteins(seInference);

            piaModeller.getProteinModeller()
                    .updateFDRData(FDRData.DecoyStrategy.SEARCHENGINE, "searchengine", proteinQThreshold);
            piaModeller.getProteinModeller().updateDecoyStates();
            piaModeller.getProteinModeller().calculateFDR();
//            piaModeller.getPeptideModeller().calculateFDR(0L);

        }

        return piaModeller;
    }

    /**
     * Compute the PSM FDR as PSM and Protein level
     * @param assayKey Assay Key
     * @param filePath File path of the assay
     * @param fileType File type, it can be mzTab, mzidentml or PRIDE xml
     * @return PIAModeller
     * @throws IOException
     */
    private PIAModeller computeFDRPSMLevel(String assayKey, String filePath, SubmissionPipelineUtils.FileType fileType) throws IOException {

        PIAModeller piaModeller = null;
        PIACompiler piaCompiler = new PIASimpleCompiler();

        String type = InputFileParserFactory.InputFileTypes.MZTAB_INPUT.getFileTypeShort();
        if(fileType == SubmissionPipelineUtils.FileType.MZID)
           type = InputFileParserFactory.InputFileTypes.MZIDENTML_INPUT.getFileTypeShort();

        piaCompiler.getDataFromFile(assayKey, filePath, null, type);

        piaCompiler.buildClusterList();
        piaCompiler.buildIntermediateStructure();

        if (piaCompiler.getAllPeptideSpectrumMatcheIDs() != null
                && !piaCompiler.getAllPeptideSpectrumMatcheIDs().isEmpty()) {

            File inferenceTempFile = File.createTempFile(assayKey, ".tmp");
            piaCompiler.writeOutXML(inferenceTempFile);
            piaCompiler.finish();
            piaModeller = new PIAModeller(inferenceTempFile.getAbsolutePath());

            if (inferenceTempFile.exists()) {
                inferenceTempFile.deleteOnExit();
            }
        }
        return piaModeller;
    }

}
