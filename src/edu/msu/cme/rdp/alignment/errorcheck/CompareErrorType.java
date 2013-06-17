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
package edu.msu.cme.rdp.alignment.errorcheck;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAligner;
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAlignment;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

public class CompareErrorType {

    public static class CountResult {

        int backCount;
        int forwardCount;
    }
    private static final char gapChar = '-';
    private PrintStream misMatch_writer = null;
    private PrintStream indel_writer = null;
    private PrintStream qualOut = null;
    private IndexedSeqReader qualReader = null;

    public CompareErrorType(File mismatch_out, File indel_out) throws IOException {
        this(mismatch_out, indel_out, null, null);
    }

    public CompareErrorType(File mismatch_out, File indel_out, File qualFile, File qualOutFile) throws IOException {
        misMatch_writer = new PrintStream(mismatch_out);
        indel_writer = new PrintStream(indel_out);

        if (qualFile != null && qualFile.exists()) {
            qualOut = new PrintStream(qualOutFile);
            qualReader = new IndexedSeqReader(qualFile, false, false);
        }

    }

    public void processSequence(String querySeqId, String querySeq, String refSeqId, String refSeq, int refStart, boolean reversed) throws IOException {
        int[][] positionMapping = mapSequences(refSeq, refStart, querySeq);

        String[] qualInfo = null;
        if (qualOut != null) {
            Object[] tmp = mapQualScores(qualReader, querySeqId, querySeq, reversed);
            qualInfo = (String[]) tmp[0];
            float avgQScore = (Float) tmp[1];
            float avgEQScore = (Float) tmp[2];

            qualOut.println(querySeqId + "\t" + avgQScore + "\t" + avgEQScore);
        }

        compare(querySeqId, refSeqId, refSeq, querySeq, qualInfo, positionMapping[0], positionMapping[1]);
    }

    public void close() throws IOException {
        misMatch_writer.close();
        indel_writer.close();

        if (qualOut != null) {
            qualOut.close();
        }

        if (qualReader != null) {
            qualReader.close();
        }
    }

    private int[][] mapSequences(String alignedRefSeq, int refStart, String alignedQuerySeq) {
        int[][] ret = new int[2][alignedRefSeq.length()];
        int refLoc = 0;
        int queryLoc = 0;

        char[] ref = alignedRefSeq.toCharArray();
        char[] query = alignedQuerySeq.toCharArray();

        for (int index = 0; index < ref.length; index++) {
            if (ref[index] != '-') {
                refLoc++;
            }

            if (query[index] != '-') {
                queryLoc++;
            }

            ret[0][index] = queryLoc;
            ret[1][index] = refLoc + refStart;
        }

        return ret;
    }

    private Object[] mapQualScores(IndexedSeqReader qualReader, String seqName, String alignedSeqString, boolean reversed) throws IOException {
        String[] ret = new String[alignedSeqString.length()];

        String dealignedSeqString = SeqUtils.getUnalignedSeqString(alignedSeqString);
        String[] qualScores = qualReader.readSeq(seqName).getSeqString().trim().split("\\s+");

        if(qualScores.length != dealignedSeqString.length()) {
            throw new IOException(seqName + "'s quality string [" + qualScores.length + "] and seq string lengths [" + dealignedSeqString.length() + "] do not match");
        }

        if(reversed) {
            int start = 0;
            int end = qualScores.length - 1;

            while(start < end) {
                String tmp = qualScores[start];
                qualScores[start] = qualScores[end];
                qualScores[end] = tmp;
                start++;
                end--;
            }
        }

        float avgEQual = 0;
        int numQualScores = 0;

        int scoreIndex = 0;
        for (int index = 0; index < dealignedSeqString.length(); index++) {
            while (alignedSeqString.charAt(scoreIndex) == '-') {
                ret[scoreIndex++] = "-1";
            }

            ret[scoreIndex++] = qualScores[index];
            numQualScores++;
            avgEQual += Math.pow(10, -1 * (Double.valueOf(qualScores[index]) / 10.0));
        }

        avgEQual /= numQualScores;
        float avgQual = -10 * (float) Math.log10(avgEQual);

        return new Object[]{ret, avgQual, avgEQual};
    }

    private String compare(String querySeqName, String refSeqName, String r, String q, String[] qualInfo, int[] queryPosMapping, int[] refPosMapping) throws IOException {
        StringBuilder buf = new StringBuilder();
        char badChar = gapChar;

        for (int i = 0; i < r.length();) {
            char rc = r.charAt(i);
            char qc = q.charAt(i);
            if (rc != gapChar && qc != gapChar) {
                if (rc != qc) {
                    misMatch_writer.println(querySeqName + "\t" + refSeqName + "\t" + (i + 1) + "\t" + rc + "\t" + qc + "\t" + queryPosMapping[i] + "\t" + refPosMapping[i] + "\t" + ((qualInfo == null) ? "null" : qualInfo[i]));
                }
                i++;
                continue;
            } else {
                if (rc == gapChar) {
                    badChar = qc;
                } else {
                    badChar = rc;
                }


                // refseq
                CountResult refCount = count(r, i, badChar);
                CountResult queryCount = count(q, i, badChar);

                indel_writer.println(querySeqName + "\t" + refSeqName + "\t" + (i + 1) + "\t" + (refCount.forwardCount + refCount.backCount) + "\t" + (queryCount.forwardCount + queryCount.backCount) + "\t" + badChar + "\t" + queryPosMapping[i] + "\t" + refPosMapping[i] + "\t" + ((qualInfo == null) ? "null" : qualInfo[i]));

                i++;
                
            }

        }
        return buf.toString();
    }

    private static CountResult count(String r, int i, char badChar) {
        // going back
        int backCount = 0;
        for (int j = i; j >= 0; j--) {
            if ( r.charAt(j) == gapChar){  // if it's gap, find the next non-gap character
                continue;
            }
            if (r.charAt(j) == badChar) {
                backCount++;
            } else {
                break;
            }
        }
        // going forward
        int forwardCount = 0;
        for (int j = i + 1; j < r.length(); j++) {
            if ( r.charAt(j) == gapChar){  // if it's gap, find the previous non-gap character
                continue;
            }
            if (r.charAt(j) == badChar) {
                forwardCount++;
            } else {
                break;
            }
        }

        CountResult result = new CountResult();
        result.backCount = backCount;
        result.forwardCount = forwardCount;
        return result;
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 5 && args.length != 7) {
            System.err.println("USAGE: CompareErrorType <ref_nucl.fasta> <query_nucl.fasta> <pairwise alignment out> <mismatch.out> <indel.out> [<unaligned_qual_file> <qual.out>]");
            return;
        }

        File refSeqFile = new File(args[0]);
        File querySeqFile = new File(args[1]);
        File alignOutFile = new File(args[2]);
        File mismatchOutFile = new File(args[3]);
        File indelOutFile = new File(args[4]);

        File qualInFile = null;
        File qualOutFile = null;

        if (args.length == 7) {
            qualInFile = new File(args[5]);
            if(!qualInFile.exists()) {
                throw new IllegalArgumentException("Quality [" + args[5] + "] file doesn't exist");
            }
            qualOutFile = new File(args[6]);
        }

        //ScoringMatrix scoringMatrix = ScoringMatrix.getDefaultNuclMatrix();
        // use a simple scoring function, match score 0, mismatch -1, gap opening -1, gap extension -1.
        ScoringMatrix scoringMatrix = new ScoringMatrix(ScoringMatrix.class.getResourceAsStream("/data/simple_scoringmatrix.txt"), -1, -1);

        List<Sequence> refSeqList = SequenceReader.readFully(refSeqFile);
        SequenceReader queryReader = new SequenceReader(querySeqFile);
        PrintStream alignOutStream = new PrintStream(alignOutFile);
        CompareErrorType errorProcessor = new CompareErrorType(mismatchOutFile, indelOutFile, qualInFile, qualOutFile);

        Sequence seq;

        while ((seq = queryReader.readNextSequence()) != null) {
            try {
                PairwiseAlignment bestResult = null;
                Sequence bestSeq = null;
                boolean bestReversed = false;
                String querySeqStr = seq.getSeqString().toLowerCase();
                String reversedQuery = IUBUtilities.reverseComplement(querySeqStr);


                for (Sequence refSeq : refSeqList) {
                    String refSeqStr = refSeq.getSeqString().toLowerCase();
                    PairwiseAlignment result = PairwiseAligner.align( refSeqStr, querySeqStr, scoringMatrix, AlignmentMode.global);
                    PairwiseAlignment reversedResult = PairwiseAligner.align(refSeqStr, IUBUtilities.reverseComplement(querySeqStr), scoringMatrix, AlignmentMode.global);


                    PairwiseAlignment currBest = (result.getScore() > reversedResult.getScore()) ? result : reversedResult;

                    if (bestResult == null || currBest.getScore() > bestResult.getScore()) {
                        bestResult = currBest;
                        bestSeq = refSeq;
                        if (currBest == reversedResult) {
                            bestReversed = true;
                        } else {
                            bestReversed = false;
                        }
                    }
                }

                int refStart = bestResult.getStarti();
                int refEnd = bestResult.getEndi();

                alignOutStream.println(">\t" + seq.getSeqName() + "\t" + bestSeq.getSeqName() + "\t" + refStart + "\t" + refEnd + "\t" + bestResult.getScore() + "\t" + ((bestReversed) ? "\treversed" : ""));
                alignOutStream.println(bestResult.getAlignedSeqj() + "\n");
                alignOutStream.println(bestResult.getAlignedSeqi() + "\n");

                //seqi is reference seq, seqj is the refseq
                errorProcessor.processSequence(seq.getSeqName(), bestResult.getAlignedSeqj(), bestSeq.getSeqName(), bestResult.getAlignedSeqi(), refStart, bestReversed);
            } catch (Exception e) {
                throw new RuntimeException("Failed while processing seq " + seq.getSeqName(), e);
            }
        }
        queryReader.close();
        alignOutStream.close();
        errorProcessor.close();
    }
}
