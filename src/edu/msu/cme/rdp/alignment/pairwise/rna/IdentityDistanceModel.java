/*
 * Created on Sep 29, 2005
 */
package edu.msu.cme.rdp.alignment.pairwise.rna;

import edu.msu.cme.rdp.readseq.utils.SeqUtils;

/**
 * @author farrisry
 */
public class IdentityDistanceModel extends DistanceModel {

    private static final byte DASH = (byte) '-';
    private static final byte DOT = (byte) '.';
    private boolean metric = false;

    public IdentityDistanceModel() {
    }

    public IdentityDistanceModel(boolean metric) {
        this.metric = metric;
    }

    /* (non-Javadoc)
     * @see edu.msu.cme.rdp.pairwisedistance.DistanceModel#getDistance(edu.msu.cme.rdp.pairwisedistance.DivergenceMatrix)
     */
    public double getDistance(byte[] seqX, byte[] seqY, int overlapMin) throws OverlapCheckFailedException {
        int match = 0;
        int compPositions = 0;
        for (int index = 0; index < seqX.length; index++) {
            if (!metric
                    && (seqX[index] == DASH || seqX[index] == SeqUtils.GAP || seqX[index] == DOT
                    || seqY[index] == DASH || seqY[index] == SeqUtils.GAP || seqY[index] == DOT)) {
                continue;
            }

            compPositions++;
            if (seqX[index] == seqY[index]) {
                match++;
            }
        }

        if (compPositions < overlapMin) {
            throw new OverlapCheckFailedException("");
        }

        return 1 - ((float)match / compPositions);
    }

    public static void main(String [] args) throws Exception {
        String s1 = "--------KIAVYENETPLFVSNIKHSVEELSAFPEVIDQFEFRKNLVLQELENNKIPF-SFDAIIGRGGLVKPIPGGVYEVNEAMKRDTVHAMR-THACNLGGLIASELASTLPCPAFIADPGVVDELEDIARITGSPLMPKIT--------------------";
        String s2 = "---------IAVYENETPLFVSNIKHSVEELSAFPEVIDQFEFRKNLVLQELENNKIPF-SFDAIIGRGGLVKPIPGGVYEVNEAMKRDTVHAMR-THACNLGGLIASELASTLPCPAFIADPGVVDELEDIARITGSPLMPKITI-------------------";

        System.out.println(new IdentityDistanceModel().getDistance(s1.getBytes(), s2.getBytes(), 0));
        System.out.println(new IdentityDistanceModel(true).getDistance(s1.getBytes(), s2.getBytes(), 0));
    }
}
