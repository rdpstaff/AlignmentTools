package edu.msu.cme.rdp.alignment.pairwise.rna;
//
//  DivergenceMatrix.java

import edu.msu.cme.rdp.readseq.utils.SeqUtils;

//
//  Created by Ryan Farris on Fri Jun 07 2002.
//  Copyright (c) 2002 Michigan State University Board of Trustees. All rights reserved.
//

public class DivergenceMatrix {
    // constants other classes will use when accessing positions on the matrix.
    protected final int[] A = {0, 0}; protected final int[] B = {1, 0}; protected final int[] C = {2, 0}; protected final int[] D = {4, 0};
    protected final int[] E = {0, 1}; protected final int[] F = {1, 1}; protected final int[] G = {2, 1}; protected final int[] H = {4, 1};
    protected final int[] I = {0, 2}; protected final int[] J = {1, 2}; protected final int[] K = {2, 2}; protected final int[] L = {4, 2};
    protected final int[] M = {0, 4}; protected final int[] N = {1, 4}; protected final int[] O = {2, 4}; protected final int[] P = {4, 4};

    // useful local constants
    private static final int x = 0;
    private static final int y = 1;

    //This matrix is 5 by 5 so that the actual value from the byte array (shifted right 1) can be used to index
    //in to the frequency matrix, this however does mean that row 3 is not used
    private double[][] matrix = new double[5][5];
    
    public DivergenceMatrix( byte[] seqX, byte[] seqY, int overlapLimit ) throws OverlapCheckFailedException {
        int matchablePositions = 0;

        double [][] localMatrix = new double[matrix.length][matrix.length];

        // alignment position by alignment position, count the number of matches and mismatches.
        // the alignments had better be the same length, 'cause I'm not checking here.
        for ( int i=0; i < seqX.length; i++ ) {
            int x = seqX[i], y = seqY[i];

            if(((x & SeqUtils.NON_COMPAREABLE) == 0) && ((y & SeqUtils.NON_COMPAREABLE) == 0)) {
                matchablePositions++;
                localMatrix[x >> 1][y >> 1]++;
            }
        }
        
        if (matchablePositions < overlapLimit) {
        		throw new OverlapCheckFailedException("did not pass the overlap limit");
        }

        // turn our matrix of counts into a matrix of frequencies.
        for ( int i=0; i < localMatrix.length; i++ ) {
            for ( int j=0; j < localMatrix.length; j++ ) {
                matrix[i][j] = localMatrix[i][j] / (double)matchablePositions;
            }
        }
    }


    public double getFrequency( int[] position ) {
        return this.matrix[position[x]][position[y]];
    }
}
