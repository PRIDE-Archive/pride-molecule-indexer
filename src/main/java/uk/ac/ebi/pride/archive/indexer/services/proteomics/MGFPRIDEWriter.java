package uk.ac.ebi.pride.archive.indexer.services.proteomics;

import lombok.extern.slf4j.Slf4j;
import uk.ac.ebi.pride.archive.dataprovider.data.spectra.BinaryArchiveSpectrum;

import java.io.IOException;
import java.io.Serializable;

@Slf4j
public class MGFPRIDEWriter implements Serializable {

    public static void appendSpectrum(Appendable out, BinaryArchiveSpectrum spectrum) {
        try {
            out.append("BEGIN IONS");
            out.append("\n");

            appendTitle(spectrum, out);
            out.append("\n");

            double precursorCharge = spectrum.getPrecursorCharge();
            double massChargeRatio = spectrum.getPrecursorMz();

            out.append("PEPMASS=").append(String.valueOf(massChargeRatio));
            out.append("\n");

            out.append("CHARGE=").append(String.valueOf(precursorCharge));
            if (precursorCharge > 0)
                out.append("+");
            out.append("\n");

//            appendProperties(spectrum, out);

            appendPeaks(spectrum, out);

            out.append("END IONS");
            out.append("\n");
        } catch (IOException e) {
            log.error("Error to write the spectrum --- " + spectrum.getUsi());
        }
    }

//    public void appendProperties(ISpectrum spectrum, Appendable out) {
//        final Properties properties = spectrum.getProperties();
//        try {
//            for (String s : properties.stringPropertyNames()) {
//                final String property = properties.getProperty(s);
//                final String line = KnownProperties.toMGFLine(s, property);
//                out.append(line);
//                out.append("\n");
//            }
//        } catch (IOException e) {
//            throw new UnsupportedOperationException(e);
//        }
//    }

    /**
     * override to add peptide later
     *
     * @param out
     * @throws IOException
     */
    public static void appendTitle(final BinaryArchiveSpectrum spectrum, final Appendable out) throws IOException {
        out.append("TITLE=id=").append(spectrum.getUsi());
        final String peptide = spectrum.getPeptidoform();
        if (peptide != null && peptide.length() > 0)
            out.append(",sequence=").append(peptide);
    }

    protected static void appendPeaks(final BinaryArchiveSpectrum spectrum, final Appendable out) throws IOException {
        for(int i=0; i < spectrum.getIntensities().length; i++){
            String line = String.format("%10.3f", spectrum.getMasses()[i]) + "\t" +
                    String.format("%10.3f", spectrum.getIntensities()[i]).trim();
            out.append(line);
            out.append("\n");
        }
    }
}
