/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.alignment.pairwise;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author fishjord
 */
public class PairwiseAligner {

    public static void reverse(StringBuffer s) {
        int l = s.length() - 1;
        for (int index = 0; index < (l + 1) / 2; index++) {
            char c = s.charAt(index);
            s.setCharAt(index, s.charAt(l - index));
            s.setCharAt(l - index, c);
        }
    }
    private static final int match = 0;
    private static final int gap = 1;
    private static final int left = 0;
    private static final int up = 1;
    private static final int diag = 2;
    private static final int trace = 2; //trace 0 = -1, -1 | trace 1 = -1, 0 | trace 2 = 0, -1

    private static int[][][] populateMatrix(char[] seqi, char[] seqj, ScoringMatrix matrix, AlignmentMode mode) {
        int[][][] scoringMatrix = new int[seqi.length + 1][seqj.length + 1][3];

        scoringMatrix[0][0][match] = 0;
        scoringMatrix[0][0][gap] = Integer.MIN_VALUE;

        for (int index = 1; index < seqi.length + 1; index++) {
            int initScore = (mode == AlignmentMode.global/* || mode == AlignmentMode.glocal*/) ? (index - 1) * matrix.getGapExtend() + matrix.getGapOpen() : 0;

            scoringMatrix[index][0][match] = initScore;
            scoringMatrix[index][0][gap] = initScore + matrix.getGapExtend();
            scoringMatrix[index][0][trace] = left;
        }
        for (int index = 1; index < seqj.length + 1; index++) {
            int initScore = (mode == AlignmentMode.global || mode == AlignmentMode.glocal) ? (index - 1) * matrix.getGapExtend() + matrix.getGapOpen() : 0;

            scoringMatrix[0][index][match] = initScore;
            scoringMatrix[0][index][gap] = initScore + matrix.getGapExtend();
            scoringMatrix[0][index][trace] = up;
        }

        for (int i = 1; i < seqi.length + 1; i++) {

            for (int j = 1; j < seqj.length + 1; j++) {

                int sxy = matrix.score(seqi[i - 1], seqj[j - 1]);

                int scoreUp = Math.max(scoringMatrix[i - 1][j][gap] + matrix.getGapExtend(), scoringMatrix[i - 1][j][match] + matrix.getGapOpen());
                int scoreLeft = Math.max(scoringMatrix[i][j - 1][gap] + matrix.getGapExtend(), scoringMatrix[i][j - 1][match] + matrix.getGapOpen());
                int m = Math.max(scoringMatrix[i - 1][j - 1][match], scoringMatrix[i - 1][j - 1][gap]) + sxy;


                if (mode == AlignmentMode.local) {
                    scoringMatrix[i][j][match] = (m < 0)? 0 : m;
                    scoringMatrix[i][j][gap] = Math.min(0, Math.max(scoreLeft, scoreUp));
                } else {
                    scoringMatrix[i][j][match] = m;
                    scoringMatrix[i][j][gap] = Math.max(scoreLeft, scoreUp);
                }

                if (m >= scoreLeft && m >= scoreUp) {
                    scoringMatrix[i][j][trace] = diag;
                } else if (scoreLeft >= m && scoreLeft >= scoreUp) {
                    scoringMatrix[i][j][trace] = up;
                } else {
                    scoringMatrix[i][j][trace] = left;
                }
            }
        }

        return scoringMatrix;
    }

    private static PairwiseAlignment traceback(int[][][] scoringMatrix, char[] seqi, char[] seqj, AlignmentMode mode) {
        StringBuffer alignedSeqi = new StringBuffer();
        StringBuffer alignedSeqj = new StringBuffer();

        List<Integer> scores = new ArrayList();

        int i = seqi.length, j = seqj.length;
        boolean fillFromJ = false;

        switch (mode) {
            case overlap:
            case overlap_trim: //For overlap we find the best score on either edge
            {
                int bestEdge = Integer.MIN_VALUE;
                for (int index = 1; index < seqi.length + 1; index++) {
                    if (scoringMatrix[index][seqj.length][match] > bestEdge) {
                        bestEdge = scoringMatrix[index][seqj.length][match];
                        i = index;
                    }
                }

                for (int index = 1; index < seqj.length + 1; index++) {
                    if (scoringMatrix[seqi.length][index][match] > bestEdge) {
                        bestEdge = scoringMatrix[seqi.length][index][match];
                        i = seqi.length;
                        j = index;
                        fillFromJ = true;
                    }
                }

            }
            break;

            case local: //For local, find the best score anywhere in the scoringMatrix
            {
                int bestEdge = Integer.MIN_VALUE;
                for (int row = 1; row < seqi.length + 1; row++) {
                    for (int col = 1; col < seqj.length + 1; col++) {
                        if (scoringMatrix[row][col][match] > bestEdge) {
                            bestEdge = scoringMatrix[row][col][match];
                            i = row;
                            j = col;
                        }
                    }
                }

            }
            break;

            case glocal: //For glocal we find the best score anywhere on the right edge of the scoringMatrix
            {
                int bestEdge = Integer.MIN_VALUE;
                for (int index = 1; index < seqi.length + 1; index++) {
                    if (scoringMatrix[index][seqj.length][match] > bestEdge) {
                        bestEdge = scoringMatrix[index][seqj.length][match];
                        i = index;
                    }
                }
                /*for (int index = 1; index < seqj.length + 1; index++) {
                    if (scoringMatrix[seqi.length][index][match] > bestEdge) {
                        bestEdge = scoringMatrix[seqi.length][index][match];
                        j = index;
                    }
                }*/
            }
            break;
        }

        if (mode == AlignmentMode.overlap || mode == AlignmentMode.glocal) {
            int x;
            char[] appendBases;
            StringBuffer appendBasesSeq, appendGaps;

            if (fillFromJ) {
                x = j;
                appendBases = seqj;
                appendBasesSeq = alignedSeqj;
                appendGaps = alignedSeqi;
            } else {
                x = i;
                appendBases = seqi;
                appendBasesSeq = alignedSeqi;
                appendGaps = alignedSeqj;
            }


            for (int index = appendBases.length - 1; index >= x; index--) {
                appendBasesSeq.append(appendBases[index]);
                appendGaps.append('-');
                if (appendBases == seqi) {
                    scores.add(scoringMatrix[i][seqj.length][match]);
                } else {
                    scores.add(scoringMatrix[seqi.length][j][match]);
                }
            }
        }

        int endi = i, endj = j;
        boolean done = false;
        while (!done) {
            int traceVal = scoringMatrix[i][j][trace];
            if (traceVal == diag) {
                alignedSeqi.append(seqi[i - 1]);
                alignedSeqj.append(seqj[j - 1]);

                i--;
                j--;
            } else if (traceVal == up) {
                alignedSeqi.append('-');
                alignedSeqj.append(seqj[j - 1]);

                j--;
            } else if (traceVal == left) {
                alignedSeqi.append(seqi[i - 1]);
                alignedSeqj.append('-');

                i--;
            } else {
                throw new IllegalArgumentException("Unknown trace value " + traceVal);
            }

            if (mode == AlignmentMode.local && scoringMatrix[i][j][match] < 0) {
                scores.add(0);
            } else {
                scores.add(scoringMatrix[i][j][match]);
            }

            switch (mode) {
                case global:
                    done = (i == 0) && (j == 0);
                    break;
                case local:
                    done = ((i == 0) && (j == 0)) || (scoringMatrix[i][j][match] <= 0);
                    break;
                case overlap:
                case overlap_trim:
                    done = (i == 0) || (j == 0);
                    break;
                case glocal:
                    //done = j == 0;
                    done = i == 0;
                    break;
                default:
                    throw new IllegalArgumentException("Unknown alignment mode " + mode);
            }
        }
        int starti = i;
        int startj = j;

        if (mode == AlignmentMode.overlap){// || mode == AlignmentMode.glocal) {
            int x;
            char[] appendBases;
            StringBuffer appendBasesSeq, appendGaps;

            if (i == 0) {
                x = j;
                appendBases = seqj;
                appendBasesSeq = alignedSeqj;
                appendGaps = alignedSeqi;
            } else {
                x = i;
                appendBases = seqi;
                appendBasesSeq = alignedSeqi;
                appendGaps = alignedSeqj;
            }

            while (x > 0) {
                appendGaps.append('-');
                appendBasesSeq.append(appendBases[x - 1]);
                scores.add(0);
                x--;
            }
        }

        reverse(alignedSeqi);

        reverse(alignedSeqj);

        Collections.reverse(scores);
        return new PairwiseAlignment(alignedSeqi.toString().toUpperCase(), alignedSeqj.toString().toUpperCase(), scores, starti, endi, startj, endj);
    }

    /**
     *
     * @param seq1 For glocal and overlap modes this sequence is assumed to be the REFERENCE
     * @param seq2 For glocal and overlap modes this seuqence is assumed to be the QUERY
     * @param scoringMatrix Scoring matrix to use, the ScoringMatrix class has defaults for Nucleotide and Protein scoring matrices
     * @param mode Alignment mode, global, local, global-local, overlap, or overlap trim
     * @return
     */
    public static PairwiseAlignment align(String seq1, String seq2, ScoringMatrix scoringMatrix, AlignmentMode mode) {
        int[][][] scores = populateMatrix(seq1.toCharArray(), seq2.toCharArray(), scoringMatrix, mode);
        return traceback(scores, seq1.toCharArray(), seq2.toCharArray(), mode);
    }
}
