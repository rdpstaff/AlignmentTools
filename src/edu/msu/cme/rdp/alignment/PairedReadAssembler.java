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

package edu.msu.cme.rdp.alignment;


import edu.msu.cme.rdp.alignment.pairwise.PairwiseAligner;
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAlignment;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.File;
import java.io.IOException;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class PairedReadAssembler {
    
    public static class AssemblyResult {
        private Sequence assembledSeq;
        private int overlapStart;
        private int overlapEnd;
        private int overlapErrors;

        public AssemblyResult(Sequence assembledSeq, int overlapStart, int overlapEnd, int overlapErrors) {
            this.assembledSeq = assembledSeq;
            this.overlapStart = overlapStart;
            this.overlapEnd = overlapEnd;
            this.overlapErrors = overlapErrors;
        }

        public Sequence getAssembledSeq() {
            return assembledSeq;
        }

        public int getOverlap() {
            return overlapEnd - overlapStart;
        }

        public int getOverlapEnd() {
            return overlapEnd;
        }

        public int getOverlapStart() {
            return overlapStart;
        }

        public int getOverlapErrors() {
            return overlapErrors;
        }
    }
    
    private static enum OverlapMode {left, overlap, right};

    private ScoringMatrix scoringMatrix;

    
    public PairedReadAssembler(ScoringMatrix scoringMatrix) {
        this.scoringMatrix = scoringMatrix;        
    }

    /**
     * Care must be taken when using this, the left read MUST be the forward sense read and the right read MUST be the reverse sense read
     * The right read will be reverse and overlap pairwise aligned with the left read, and overlap determined from that
     *
     * @param left
     * @param right
     * @return
     */
    public AssemblyResult assemble(Sequence left, Sequence right) {
        PairwiseAlignment align = PairwiseAligner.align(left.getSeqString(), IUBUtilities.reverseComplement(right.getSeqString()), scoringMatrix, AlignmentMode.overlap);

        char[] alignedLeft = align.getAlignedSeqi().toCharArray();
        char[] alignedRight = align.getAlignedSeqj().toCharArray();

        int overlapStart = -1;
        int overlapEnd = -1;
        int overlapErrors = 0;

        OverlapMode overlap = OverlapMode.left;
        StringBuilder assembledSeqStr = new StringBuilder();

        for(int index = 0;index < alignedLeft.length;index++) {
            if(alignedLeft[index] == '-' && overlapEnd == -1) {
                overlapEnd = index;
                overlap = OverlapMode.right;
            }

            if(alignedRight[index] != '-' && overlapStart == -1) {
                overlapStart = index;
                overlap = OverlapMode.overlap;
            }

            if(overlap == OverlapMode.overlap && alignedLeft[index] != alignedRight[index]) {
                overlapErrors++;
                assembledSeqStr.append("N");
                continue;
            }

            switch(overlap) {
                case left:
                case overlap:
                    assembledSeqStr.append(alignedLeft[index]);
                    break;
                case right:
                    assembledSeqStr.append(alignedRight[index]);
                    break;
            }
        }

        String assembledId = left.getSeqName();
        if(assembledId.contains(".")) {
            assembledId = assembledId.substring(0, assembledId.lastIndexOf("."));
        }

        /*System.out.println(overlapStart + "\t" + overlapEnd + "\t" + (overlapEnd - overlapStart));
        System.out.println(align.getAlignedSeqi());
        System.out.println(align.getAlignedSeqj());
        System.out.println(assembledSeqStr);
        System.out.println();*/

        return new AssemblyResult(new Sequence(assembledId, left.getDesc(), assembledSeqStr.toString()), overlapStart, overlapEnd, overlapErrors);
    }

    public static void main(String [] args) throws Exception {

        if(args.length != 5) {
            System.err.println("USAGE: PairedReadAssembler <fseq_file> <rseq_file> <min_ident> <min_overlap> <assembled_out_file>");
            return;
        }

        File fSeqFile = new File(args[0]);
        File rSeqFile = new File(args[1]);
        float maxError = 1 - Float.parseFloat(args[2]);
        int minOverlap = Integer.parseInt(args[3]);
        FastaWriter seqOut = new FastaWriter(new File(args[4]));

        PairedReadAssembler assembler = new PairedReadAssembler(ScoringMatrix.getDefaultNuclMatrix());

        SequenceReader fSeqReader = new SequenceReader(fSeqFile);
        IndexedSeqReader rSeqReader = new IndexedSeqReader(rSeqFile);
        Set<String> rseqIds = rSeqReader.getSeqIdSet();

        System.out.println("fseq_id\trseq_id\tassembled_length\toverlap_start\toverlap_end\ttotal_overlap\tmismatches\terror_ratio\tassembled?");

        Sequence fseq;
        while((fseq = fSeqReader.readNextSequence()) != null) {
            if(!rseqIds.contains(fseq.getSeqName())) {
                System.out.println(fseq.getSeqName() + "\tnot found in rseq file");
                continue;
            }

            Sequence rseq = rSeqReader.readSeq(fseq.getSeqName());

            if(!fseq.getSeqName().equals(rseq.getSeqName())) {
                throw new IOException("Read " + fseq.getSeqName() + " from fseq file and " + rseq.getSeqName() + " from rseq file, can't match sequences with different names");
            }

            AssemblyResult result = assembler.assemble(fseq, rseq);

            float errorRatio = (float)result.getOverlapErrors() / result.getOverlap();
            boolean assembled = false;
            if(result.getOverlap() > minOverlap && errorRatio < maxError) {
                seqOut.writeSeq(result.getAssembledSeq());
                assembled = true;
            }

            System.out.println(fseq.getSeqName() + "\t" + result.getAssembledSeq().getSeqString().length() + "\t" + result.getOverlapStart() + "\t" + result.getOverlapEnd() + "\t" + result.getOverlap() + "\t" + result.getOverlapErrors() + "\t" + errorRatio + "\t" + assembled);
        }
        seqOut.close();
        fSeqReader.close();
        rSeqReader.close();
    }

}
