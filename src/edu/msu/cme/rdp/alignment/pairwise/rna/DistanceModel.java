package edu.msu.cme.rdp.alignment.pairwise.rna;
//  CorrectedDistanceModel.java
//  SequenceRelatedness
//
//  Created by Ryan Farris on Fri Jun 07 2002.
//  Copyright (c) 2002 Michigan State University Board of Trustees. All rights reserved.
//

public abstract class DistanceModel {

    public abstract double getDistance( byte[] seqX, byte[] seqY, int overlapLimit ) throws OverlapCheckFailedException;
}
