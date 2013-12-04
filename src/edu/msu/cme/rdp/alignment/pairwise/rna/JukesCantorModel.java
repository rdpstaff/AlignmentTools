package edu.msu.cme.rdp.alignment.pairwise.rna;

//
//  JukesCantorModel.java
//  SequenceRelatedness
//
//  Created by Ryan Farris on Fri Jun 07 2002.
//  Copyright (c) 2002 Michigan State University Board of Trustees. All rights reserved.
//

public class JukesCantorModel extends DistanceModel {

    public final double getDistance( byte[] seqX, byte[] seqY, int overlapLimit) throws OverlapCheckFailedException {
    	DivergenceMatrix dm = new DivergenceMatrix( seqX, seqY, overlapLimit );
        double D = 1.0 - ( dm.getFrequency( dm.A ) + dm.getFrequency( dm.F ) + dm.getFrequency( dm.K ) + dm.getFrequency( dm.P ));
        double correctedDistance = -(3.0/4.0) * Math.log( 1.0 - (4.0/3.0) * D );
        return correctedDistance;
    }
}
