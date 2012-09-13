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
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAligner;
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAlignment;
import edu.msu.cme.rdp.alignment.pairwise.rna.DistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.IdentityDistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.io.File;
import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author fishjord
 */
public class PairwiseMatcher {

    public static void main(String[] args) throws Exception {
        if(args.length != 3) {
            System.err.println("USAGE: PairwiseMatcher <ref_file> <query_file> <out_file>");
            return;
        }

        File inFile = new File(args[1]);
        File refFile = new File(args[0]);

        List<Sequence> refSeqs = SequenceReader.readFully(refFile);
        ScoringMatrix matrix = ScoringMatrix.getDefaultProteinMatrix();
        DistanceModel model = new IdentityDistanceModel();

        SeqReader reader = new SequenceReader(inFile);
        Sequence seq;

        Map<String, Integer> hits = new LinkedHashMap();

        for (Sequence refSeq : refSeqs) {
            hits.put(refSeq.getSeqName(), 0);
        }

        PrintStream out = new PrintStream(args[2]);

        while ((seq = reader.readNextSequence()) != null) {
            long startTime = System.currentTimeMillis();
            PairwiseAlignment bestAlignment = null;
            Sequence bestRef = null;

            for (Sequence refSeq : refSeqs) {
                PairwiseAlignment alignment = PairwiseAligner.align(refSeq.getSeqString(), seq.getSeqString(), matrix, AlignmentMode.glocal);

                if (bestAlignment == null || bestAlignment.getScore() < alignment.getScore()) {
                    bestAlignment = alignment;
                    bestRef = refSeq;
                }
            }

            int start = -1, end = -1;
            int refPos = 0;

            char[] seqi = bestAlignment.getAlignedSeqi().toCharArray();
            char[] seqj = bestAlignment.getAlignedSeqj().toCharArray();

            for (int index = 0; index < bestAlignment.getAlignedSeqi().length(); index++) {

                if (seqj[index] != '-') {
                    if (start == -1) {
                        start = refPos;
                    }

                    end = refPos + 1;
                }

                if (seqi[index] != '-') {
                    refPos++;
                }
            }

            double ident = model.getDistance(bestAlignment.getAlignedSeqi().getBytes(), bestAlignment.getAlignedSeqj().getBytes(), 0);

            out.println(">\t" + seq.getSeqName() + "\t" + bestRef.getSeqName() + "\t" + bestAlignment.getScore() + "\t" + start + "\t" + end + "\t" + ident);

            out.println(bestAlignment.getAlignedSeqi());
            out.println(bestAlignment.getAlignedSeqj());

            hits.put(bestRef.getSeqName(), hits.get(bestRef.getSeqName()) + 1);

            System.out.println("Processed " + seq.getSeqName() + " in " + (System.currentTimeMillis() - startTime) + "ms");
        }
        out.close();

        for (String seqid : hits.keySet()) {
            int hitCount = hits.get(seqid);
            if (hitCount == 0) {
                continue;
            }

            System.out.println(seqid + "\t" + hitCount);
        }
    }
}
