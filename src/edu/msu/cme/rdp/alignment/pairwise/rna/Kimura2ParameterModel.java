package edu.msu.cme.rdp.alignment.pairwise.rna;
//  Kimura2ParameterModel.java
//  SequenceRelatedness
//
//  Created by Ryan Farris on Fri Jun 07 2002.
//  Copyright (c) 2002 Michigan State University Board of Trustees. All rights reserved.
//


public class Kimura2ParameterModel extends DistanceModel {

    // this could be static, but I want this to be part of a interface
    public double getDistance( byte[] seqX, byte[] seqY, int overlapLimit) throws OverlapCheckFailedException {
    	DivergenceMatrix dm = new DivergenceMatrix( seqX, seqY, overlapLimit );
        double P = dm.getFrequency( dm.C ) + dm.getFrequency( dm.H ) + dm.getFrequency( dm.I ) + dm.getFrequency( dm.N );
        double Q = dm.getFrequency( dm.B ) + dm.getFrequency( dm.D ) + dm.getFrequency( dm.E ) + dm.getFrequency( dm.G )
            + dm.getFrequency( dm.J ) + dm.getFrequency( dm.L ) + dm.getFrequency( dm.M ) + dm.getFrequency( dm.O );

        double correctedDistance = (.5) * Math.log( 1.0 / ( 1.0 - (2.0 * P) - Q )) + (.25) * Math.log( 1.0 / ( 1.0 - (2.0 * Q) ));

        return correctedDistance;
    }

}
