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
import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;
import edu.msu.cme.rdp.readseq.readers.QSeqReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

public class CompareErrorType {

    public static class CountResult {

        int backCount;
        int forwardCount;
    }
    
    public static class PAObject {
        PairwiseAlignment PA;
        boolean reversed;
        Sequence refSeq;
        
        public PAObject(PairwiseAlignment pa, boolean reverse, Sequence seq){
            this.PA = pa;
            this.reversed = reverse;
            this.refSeq = seq;
        }
        public PairwiseAlignment getPA(){
            return PA;
        }
        public boolean getReversed(){
            return reversed;
        }
        public Sequence getRefSeq(){
            return refSeq;
        }
    }
    private static final char gapChar = '-';
    private PrintStream misMatch_writer = null;
    private PrintStream indel_writer = null;
    private PrintStream qualOut = null;

    public CompareErrorType(File mismatch_out, File indel_out) throws IOException {
        this(mismatch_out, indel_out, null);
    }

    public CompareErrorType(File mismatch_out, File indel_out, File qualOutFile) throws IOException {
        misMatch_writer = new PrintStream(mismatch_out);
        indel_writer = new PrintStream(indel_out);

        if (qualOutFile != null) {
            qualOut = new PrintStream(qualOutFile);
        }

    }

    public void processSequence(Sequence querySeq, String alignedQuery, String refSeqId, String refSeq, int refStart, boolean reversed) throws IOException {
        int[][] positionMapping = mapSequences(refSeq, refStart, alignedQuery);

        String[] qualInfo = null;
        if (qualOut != null) {
            Object[] tmp = mapQualScores((QSequence)querySeq, alignedQuery, reversed);
            qualInfo = (String[]) tmp[0];
            float avgQScore = (Float) tmp[1];
            float avgEQScore = (Float) tmp[2];

            qualOut.println(querySeq.getSeqName() + "\t" + avgQScore + "\t" + avgEQScore);
        }

        compare(querySeq.getSeqName(), refSeqId, refSeq, alignedQuery, qualInfo, positionMapping[0], positionMapping[1]);
    }

    public void close() throws IOException {
        misMatch_writer.close();
        indel_writer.close();

        if (qualOut != null) {
            qualOut.close();
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

    private Object[] mapQualScores(QSequence seq, String alignedSeqString, boolean reversed) throws IOException {
        String[] ret = new String[alignedSeqString.length()];

        String dealignedSeqString = SeqUtils.getUnalignedSeqString(alignedSeqString);
        byte[] qualScores = seq.getQuality();

        if (qualScores.length != dealignedSeqString.length()) {
            throw new IOException(seq.getSeqName() + "'s quality string [" + qualScores.length + "] and seq string lengths [" + dealignedSeqString.length() + "] do not match");
        }

        if (reversed) {
            int start = 0;
            int end = qualScores.length - 1;

            while (start < end) {
                byte tmp = qualScores[start];
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

            ret[scoreIndex++] = qualScores[index] + "";
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
            if (r.charAt(j) == gapChar) {  // if it's gap, find the next non-gap character
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
            if (r.charAt(j) == gapChar) {  // if it's gap, find the previous non-gap character
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
        Options options = new Options();
        options.addOption("s", "stem", true, "Output stem (default <query_nucl.fasta>)");

        final SeqReader queryReader;
        final List<Sequence> refSeqList;
        final PrintStream alignOutStream;
        final CompareErrorType errorProcessor;
        Sequence seq;
        Map<String, PAObject> matchMap = new HashMap();


        try {
            CommandLine line = new PosixParser().parse(options, args);

            String stem;

            args = line.getArgs();
            if (args.length != 2 && args.length != 3) {
                throw new Exception("Unexpected number of arguments");
            }

            File refFile = new File(args[0]);
            File queryFile = new File(args[1]);

            if (line.hasOption("stem")) {
                stem = line.getOptionValue("stem");
            } else {
                stem = queryFile.getName();
            }

            File alignOutFile = new File(stem + "_alignments.txt");
            File mismatchOutFile = new File(stem + "_mismatches.txt");
            File indelOutFile = new File(stem + "_indels.txt");
            File qualOutFile = null;

            refSeqList = SequenceReader.readFully(refFile);
            if (args.length == 3) {
                queryReader = new QSeqReader(queryFile, new File(args[2]));
            } else {
                queryReader = new SequenceReader(queryFile);
            }

            seq = queryReader.readNextSequence();
            if(seq instanceof QSequence) {
                qualOutFile = new File(stem + "_qual.txt");
            }

            errorProcessor = new CompareErrorType(mismatchOutFile, indelOutFile, qualOutFile);
            alignOutStream = new PrintStream(alignOutFile);

            System.err.println("Starting CompareErrorType");
            System.err.println("*  Time:              " + new Date());
            System.err.println("*  Reference File:    " + refFile);
            System.err.println("*  Query File:        " + queryFile);
            if(args.length == 3) {
                System.err.println("*  Qual File:         " + args[2]);
            }
            System.err.println("*  Query format:      " + queryReader.getFormat());
            System.err.println("*  Alignment Output:  " + alignOutFile);
            System.err.println("*  Mismatches Output: " + mismatchOutFile);
            System.err.println("*  Alignment Output:  " + indelOutFile);
            if(qualOutFile != null) {
                System.err.println("*  Quality Output:    " + qualOutFile);
            }

        } catch (Exception e) {
            new HelpFormatter().printHelp("CompareErrorType [options] <ref_nucl> (<query_nucl> | <query_nucl.fasta> <query_nucl.qual>)", options);
            System.err.println("ERROR: " + e.getMessage());
            throw new RuntimeException(e);
            //System.exit(1);
            //return;
        }

        //ScoringMatrix scoringMatrix = ScoringMatrix.getDefaultNuclMatrix();
        // use a simple scoring function, match score 0, mismatch -1, gap opening -1, gap extension -1.
        ScoringMatrix scoringMatrix = new ScoringMatrix(ScoringMatrix.class.getResourceAsStream("/data/simple_scoringmatrix.txt"), -1, -1);

        do {
            try {
                PairwiseAlignment bestResult = null;
                Sequence bestSeq = null;
                boolean bestReversed = false;
                String querySeqStr = seq.getSeqString().toLowerCase();
                String reversedQuery = IUBUtilities.reverseComplement(querySeqStr);
                PAObject bestMatch = null;
                
                //checking if sequence has been seen before
                if(matchMap.containsKey(seq.getSeqString())){
                    bestMatch = matchMap.get(seq.getSeqString());
                }
                else{
                    for (Sequence refSeq : refSeqList) {
                            String refSeqStr = refSeq.getSeqString().toLowerCase();
                            PairwiseAlignment result = PairwiseAligner.align(refSeqStr, querySeqStr, scoringMatrix, AlignmentMode.global);
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
                     
                     //Since this is a new sequence, make a new PAObject to put into the map to compare against later
                     bestMatch = new PAObject(bestResult, bestReversed, bestSeq);
                     matchMap.put(seq.getSeqString(), bestMatch);
                    }
                }
                int refStart = bestMatch.getPA().getStarti();
                int refEnd = bestMatch.getPA().getEndi();
                bestSeq = bestMatch.getRefSeq();
                bestReversed = bestMatch.getReversed();
                bestResult = bestMatch.getPA();
                
                //output information
                alignOutStream.println(">\t" + seq.getSeqName() + "\t" + bestSeq.getSeqName() + "\t" +  seq.getSeqString().length() + "\t" + refStart + "\t" + refEnd + "\t" + bestResult.getScore() + "\t" + ((bestReversed) ? "\treversed" : ""));
                alignOutStream.println(bestResult.getAlignedSeqj() + "\n");
                alignOutStream.println(bestResult.getAlignedSeqi() + "\n");

                //seqi is reference seq, seqj is the refseq
                errorProcessor.processSequence(seq, bestResult.getAlignedSeqj(), bestSeq.getSeqName(), bestResult.getAlignedSeqi(), refStart, bestReversed);
                    
            } catch (Exception e) {
                throw new RuntimeException("Failed while processing seq " + seq.getSeqName(), e);
            }
        } while ((seq = queryReader.readNextSequence()) != null);
        queryReader.close();
        alignOutStream.close();
        errorProcessor.close();     
    }
}