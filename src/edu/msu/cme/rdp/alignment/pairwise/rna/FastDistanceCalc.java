package edu.msu.cme.rdp.alignment.pairwise.rna;
//
//  DivergenceMatrix.java

import edu.msu.cme.rdp.readseq.utils.SeqUtils;

//
//  Created by Ryan Farris on Fri Jun 07 2002.
//  Copyright (c) 2002 Michigan State University Board of Trustees. All rights reserved.
//
public class FastDistanceCalc {
    private static final float nan = -1.0f;

    public static float FastDistanceCalc(byte[] seqX, byte[] seqY, int overlapLimit) {
        int matchable = 0;
        int matches = 0;
        int x, y;

        for (int i = 0; i < seqX.length; i++) {            
            x = seqX[i];
            y = seqY[i];
            
            if(((x | y) & SeqUtils.NON_COMPAREABLE) == 0) {
                matchable++;
                if(x == y) {
                    matches++;
                }
            }
        }
        
        if(matchable == 0 || matchable < overlapLimit) {
            return nan;
        }
        return 1 - matches / (float)matchable;
    }
}
