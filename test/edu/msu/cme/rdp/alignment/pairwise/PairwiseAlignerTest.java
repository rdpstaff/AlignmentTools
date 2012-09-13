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
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author wangqion
 */
public class PairwiseAlignerTest {

    public PairwiseAlignerTest() {
    }

    /**
     * Test of reverse method, of class PairwiseAligner.
     */
    @Test
    public void testReverse() {
        StringBuffer refSeq = new StringBuffer("CACGCAGGCG");
        PairwiseAligner.reverse(refSeq);
        assertEquals(refSeq.toString(), "GCGGACGCAC");
    }

    /**
     * Test of align method, of class PairwiseAligner.
     */
    @Test
    public void testAlignProtSeq() {
        String refSeq = "TRLILNSKAQTTVMDLARERGTVEDLELEDVLVEGHLGVRCAESGGPEPGVGCAGRGVITAINFLEENGAYTEDTDYVFYDVLGDVVCGGFAMPIRENKAKEIYIVT";
        String querySeq = "ergedleledvlveghlgvrcaesggpepgvgcagrgvitainfleengayt";

        ScoringMatrix scoringMatrix = ScoringMatrix.getDefaultProteinMatrix();

        PairwiseAlignment result = PairwiseAligner.align(refSeq, querySeq, scoringMatrix, AlignmentMode.overlap_trim);
        assertEquals("ERGTVEDLELEDVLVEGHLGVRCAESGGPEPGVGCAGRGVITAINFLEENGAYT", result.getAlignedSeqi());
        assertEquals("ERG--EDLELEDVLVEGHLGVRCAESGGPEPGVGCAGRGVITAINFLEENGAYT", result.getAlignedSeqj());

        result = PairwiseAligner.align(refSeq, querySeq, scoringMatrix, AlignmentMode.glocal);
        assertEquals("TRLILNSKAQTTVMDLARERGTVEDLELEDVLVEGHLGVRCAESGGPEPGVGCAGRGVITAINFLEENGAYTEDTDYVFYDVLGDVVCGGFAMPIRENKAKEIYIVT", result.getAlignedSeqi());
        assertEquals("------------------ERG--EDLELEDVLVEGHLGVRCAESGGPEPGVGCAGRGVITAINFLEENGAYT-----------------------------------", result.getAlignedSeqj());

        result = PairwiseAligner.align(refSeq, querySeq, scoringMatrix, AlignmentMode.local);
        assertEquals("EDLELEDVLVEGHLGVRCAESGGPEPGVGCAGRGVITAINFLEENGAYT", result.getAlignedSeqi());
        assertEquals("EDLELEDVLVEGHLGVRCAESGGPEPGVGCAGRGVITAINFLEENGAYT", result.getAlignedSeqj());

        result = PairwiseAligner.align(refSeq, querySeq, scoringMatrix, AlignmentMode.global);
    }

    public void testAlignNuclSeq() {
        ScoringMatrix scoringMatrix = ScoringMatrix.getDefaultNuclMatrix();
        String refSeq = "TGTGCGGCGGCTTCGCCATGCCGATTTCGCGAAACAAGGCGCAGGAAAATCTACATCGTGATA";
        String querySeq = "TGCGCCATGCCGATTCGCGAAAACAAGGCGCAGGAAATCTACATC";

        PairwiseAlignment result = PairwiseAligner.align(refSeq, querySeq, scoringMatrix, AlignmentMode.glocal);
        assertEquals("TGTGCGGCGGCTTCGCCATGCCGATTTCGCG-AAACAAGGCGCAGGAAAATCTACATCGTGATA", result.getAlignedSeqi());
        assertEquals("-----------TGCGCCATGCCGA-TTCGCGAAAACAAGGCGCAGG-AAATCTACATC------", result.getAlignedSeqj());

        result = PairwiseAligner.align(refSeq, querySeq, scoringMatrix, AlignmentMode.overlap_trim);
        assertEquals("TTCGCCATGCCGATTTCGCG-AAACAAGGCGCAGGAAAATCTACATC", result.getAlignedSeqi());
        assertEquals("TGCGCCATGCCGA-TTCGCGAAAACAAGGCGCAGG-AAATCTACATC", result.getAlignedSeqj());

        result = PairwiseAligner.align(refSeq, querySeq, scoringMatrix, AlignmentMode.global);
        assertEquals("TGTGCGGCGGCTTCGCCATGCCGATTTCGCG-AAACAAGGCGCAGGAAAATCTACATCGTGATA", result.getAlignedSeqi());
        assertEquals("-----------TGCGCCATGCCGA-TTCGCGAAAACAAGGCGCAGG-AAATCTACATC------", result.getAlignedSeqj());
    }
}
