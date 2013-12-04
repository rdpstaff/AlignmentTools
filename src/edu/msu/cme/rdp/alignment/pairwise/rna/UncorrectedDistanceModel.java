/*
 * Created on Sep 29, 2005
 */
package edu.msu.cme.rdp.alignment.pairwise.rna;

/**
 * @author farrisry
 */
public class UncorrectedDistanceModel extends DistanceModel {

    /* (non-Javadoc)
     * @see edu.msu.cme.rdp.pairwisedistance.DistanceModel#getDistance(edu.msu.cme.rdp.pairwisedistance.DivergenceMatrix)
     */
    public double getDistance( byte[] seqX, byte[] seqY, int overlapLimit) throws OverlapCheckFailedException {
    	DivergenceMatrix dm = new DivergenceMatrix( seqX, seqY, overlapLimit );
        double distance = 1.0 - ( dm.getFrequency( dm.A ) + dm.getFrequency( dm.F ) + dm.getFrequency( dm.K ) + dm.getFrequency( dm.P ));
        return distance;
    }

}
