package uk.ac.ebi.pride.archive.indexer.services.proteomics;

import uk.ac.ebi.pride.archive.indexer.utility.StringUtils;

import java.io.Serializable;
import java.util.Objects;

public class PeptidoformClustered implements Serializable {

    private static final long serialVersionUID = -5553406563748295904L;
    private final String sequence;
    String peptidoform;
    Boolean isDecoy;

    public PeptidoformClustered(String sequence, String peptidoform, boolean isDecoy) {
        this.sequence = sequence;
        this.peptidoform = peptidoform;
        this.isDecoy = isDecoy;
    }

    public Boolean getDecoy() {
        return isDecoy;
    }

    public String getSequence() {
        return sequence;
    }

    @Override
    public boolean equals(Object o) {
        String oPep = StringUtils.makePeptideIsobaric(((PeptidoformClustered) o).peptidoform);
        String isopep = StringUtils.makePeptideIsobaric(peptidoform);
        if (isopep.equals(oPep)) return true;
        PeptidoformClustered that = (PeptidoformClustered) o;
        return peptidoform.equals(that.peptidoform);
    }

    @Override
    public int hashCode() {
        String isopep = StringUtils.makePeptideIsobaric(peptidoform);
        return Objects.hash(isopep);
    }

    @Override
    public String toString() {
        return peptidoform;
    }
}
