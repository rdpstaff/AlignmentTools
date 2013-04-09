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
package edu.msu.cme.rdp.alignment.hmm.scoring;

import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.alignment.hmm.XSTATES;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.util.ArrayList;
import java.util.List;


import static java.lang.StrictMath.*;

/**
 *
 * @author fishjord
 */
public abstract class HMMScorer {
    /*
    final public void consume(char b) {
    extend();
    seq.append(b);
    float sc;
    mmx(i, 0, Float.NEGATIVE_INFINITY);
    imx(i, 0, Float.NEGATIVE_INFINITY);
    dmx(i, 0, Float.NEGATIVE_INFINITY);
    xmx(i, XSTATES.E, Float.NEGATIVE_INFINITY);
    for (int k = 1; k < hmm.M(); k++) {
    sc = max(mmx(i - 1, k - 1) + hmm.tsc(k - 1, MM),
    imx(i - 1, k - 1) + hmm.tsc(k - 1, IM));
    sc = max(sc, dmx(i - 1, k - 1) + hmm.tsc(k - 1, DM));
    //sc = max(sc, xmx(i - 1, B) + hmm.tsc(k - 1, BM));     //So there is a small problem with this since we don't know the length at the beginning...
    mmx(i, k, sc + hmm.msc(k, b));
    //xmx[i][XSTATES.E.ordinal()] = max(xmx[i][XSTATES.E.ordinal()], matrix[M][i][k] + esc); //So there is a small problem with this since we don't know the length at the beginning...but since esc is inf for glocal we're cool
    sc = max(mmx(i - 1, k) + hmm.tsc(k, MI),
    imx(i - 1, k) + hmm.tsc(k, II));
    imx(i, k, sc + hmm.isc(k, b));
    dmx(i, k, max(mmx(i, k - 1) + hmm.tsc(k - 1, MD),
    dmx(i, k - 1) + hmm.tsc(k - 1, DD)));
    }
    sc = max(mmx(i - 1, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, MM),
    imx(i - 1, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, IM));
    sc = max(sc, dmx(i - 1, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, DM));
    //sc = max(sc, xmx(i - 1, B) + hmm.tsc(hmm.M() - 1, BM));   //So there is a small problem with this since we don't know the length at the beginning...
    mmx(i, hmm.M(), sc + hmm.msc(hmm.M(), b));
    imx(i, hmm.M(), Float.NEGATIVE_INFINITY);
    dmx(i, hmm.M(), max(mmx(i, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, MD),
    dmx(i, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, DD)));
    //E
    sc = max(xmx(i, XSTATES.E), mmx(i, hmm.M()));
    xmx(i, XSTATES.E, max(sc, dmx(i, hmm.M())));
     *
     * In theory these should not be used for unihit
     *
    //J
    sc = xmx(i - 1, J) + hmm.xsc(J, XSC.LOOP);
    xmx(i, J, max(sc, xmx(i, XSTATES.E) + hmm.xsc(XSTATES.E, XSC.LOOP)));
    //C
    sc = xmx(i - 1, C) + hmm.xsc(C, XSC.LOOP);
    xmx(i, C, max(sc, xmx(i, XSTATES.E) + hmm.xsc(XSTATES.E, XSC.MOVE)));
    //N
    xmx(i, N, xmx(i - 1, N) + hmm.xsc(N, XSC.LOOP));
    //B
     *
     * In theory J is always going to be infinite in unihit configuration
     *
     *
    sc = xmx(i, N) + hmm.xsc(N, XSC.MOVE);
    xmx(i, B, max(sc, xmx(i, J) + hmm.xsc(J, XSC.MOVE)));
    }
     */
    //protected StringBuilder seq = new StringBuilder();
    protected final ProfileHMM hmm;
    protected final int startingState;
    protected final int modelLength;
    protected final int modelStates;

    private List<double[][]> matrix;
    private static final int M = 0, I = 1, D = 2;
    protected int i;  //curr seq position
    public final static float ln2 = (float)log(2);
    private double[] bestPossibleScore;
    private List<double[]> xmx;
    private static final double[] nullModels = new double[5000];

    protected List<Double> maxScores = new ArrayList();

    static {
        for (int i = 0; i < 5000; i++) {
            nullModels[i] = genNull1(i);
        }
    }

    public HMMScorer(ProfileHMM hmm) {
        this(hmm, -1);
    }

    public HMMScorer(ProfileHMM hmm, int startingState) {
        this.hmm = hmm;
        this.modelLength = hmm.M();
        this.modelStates = hmm.M() - 1;
        //reset();

        /*hmm.reconfigureLength(475);

        /*
         * Configure the best possible remaining scores
         */
        /*String bestSeq = HMMEmit.mostProbableSeq(hmm);
        bestPossibleScore = new float[hmm.M() + 1];

        for (char c : bestSeq.toCharArray()) {
            consume(c);
        }

        if (bestSeq.length() != hmm.M()) {
            throw new IllegalStateException("Sequence length doesn't match model length!");
        }

        float score = mmx(hmm.M(), hmm.M());
        float startingScore = 0;//(startingState != -1)? mmx(startingState, startingState) : 0;

        for (int k = 0; k <= hmm.M(); k++) {
            bestPossibleScore[k] = (float)(score - mmx(k, k) - startingScore);
        }*/

        this.startingState = startingState;
        reset();
    }

    private static double genNull1(int L) {
        float p1 = (float) L / (L + 1);
        float null1 = (float)( L * log(p1) + log(1 - p1));

        return null1;
    }

    public static double getNull1(int L) {
        if (L < nullModels.length) {
            return nullModels[L];
        } else {
            return genNull1(L);
        }
    }

    protected abstract void initializeSM();
    public abstract void consume(char b);

    protected final void reset() {
        i = -1;
        matrix = new ArrayList();
        xmx = new ArrayList();
        //seq = new StringBuilder();

        maxScores = new ArrayList();
        maxScores.add(Double.NEGATIVE_INFINITY);

        extend();
        initializeSM();
    }

    protected void extend() {
        i++;
        double[][] nextRow = new double[3][hmm.M() + 1];
        matrix.add(nextRow);
        xmx.add(new double[XSTATES.values().length]);
    }

    public Object[] getBestScore() {
        int bestK = -1;
        double score = Float.NEGATIVE_INFINITY;
        double[] mx = mmx(i);
        double[] dx = dmx(i);

        for (int k = 2; k <= modelLength; k++) {
            double sc = max(mx[k], dx[k]);
            if (sc > score) {
                score = sc;
                bestK = k;
            }
            score = max(score, sc);
        }

        return new Object[]{bestK, score};
    }

    public Object[] getBestScore(int mp) {
        int bestK = mp;
        double[] mx = mmx(i);
        double[] dx = dmx(i);

        return new Object[]{bestK, max(mx[mp], dx[mp])};
    }

    public double getMaxScore() {
/*        float score = Float.NEGATIVE_INFINITY;

        float[] mx = mmx(i);
        float[] dx = dmx(i);

        for (int k = 2; k <= modelLength; k++) {
            float sc = max(mx[k], dx[k]);
            if (sc > score) {
                score = sc;
            }
            score = max(score, sc);
        }*/

        double score = maxScores.get(i);
        return (score - getNull1(i)) / ln2;
        //return new Object[] { i, (score - getNull1(i)) / ln2, score, bestK, state, hmm.tsc(bestK - 1, II), hmm.tsc(bestK - 1, MI) };
    }

    public double getMaxScorePrint() {
        double score = Float.NEGATIVE_INFINITY;
        int maxState = -1;

        double[] mx = mmx(i);
        double[] dx = dmx(i);
        boolean m = false;

        for (int k = startingState - 1; k <= modelLength; k++) {
            double sc = max(mx[k], dx[k]);
            if (sc > score) {
                m = sc == mx[k];
                score = sc;
                maxState = k;
            }
            score = max(score, sc);
        }

        System.out.println("Max state: " + maxState + " seq pos: " + i + " match state?: " + m + " score: " + score + " scores[startingState - 1]: " + mx[startingState - 1] + ", " + dx[startingState - 1] + "scores[startingState - 2]: " + mx[startingState - 2] + ", " + dx[startingState - 2]);
        System.out.println(mmx(0, 0) + ", " + dmx(0, 0) + "\t" + mmx(0, 1) + ", " + dmx(0, 1));
        System.out.println(mmx(1, 0) + ", " + dmx(1, 0) + "\t" + mmx(1, 1) + ", " + dmx(1, 1));
        System.out.println();

        return (score - getNull1(i)) / ln2;
    }

    public double getScore() {
        double score = Math.max(mmx(i, hmm.M()), dmx(i, hmm.M()));//xmx(i, XSTATES.E);

        return (score - getNull1(i)) / ln2;
    }

    public void contract() {
        if (i >= 0) {
            matrix.remove(i);
            xmx.remove(i);
            //seq.deleteCharAt(i);
            maxScores.remove(i);
            i--;
        }
    }

    /*final public String getSeq() {
        return seq.toString();
    }*/

    final public int getQueryLength() {
        return i;
    }

    final protected double mmx(int i, int k) {
        return matrix.get(i)[M][k];
    }

    final protected double imx(int i, int k) {
        return matrix.get(i)[I][k];
    }

    final protected double dmx(int i, int k) {
        return matrix.get(i)[D][k];
    }

    final protected double xmx(int i, XSTATES s) {
        return xmx.get(i)[s.ordinal()];
    }

    final protected double[] mmx(int i) {
        return matrix.get(i)[M];
    }

    final protected double[] imx(int i) {
        return matrix.get(i)[I];
    }

    final protected double[] dmx(int i) {
        return matrix.get(i)[D];
    }

    final protected double[] xmx(int i) {
        return xmx.get(i);
    }

    final protected void mmx(int i, int k, double val) {
        matrix.get(i)[M][k] = val;
    }

    final protected void imx(int i, int k, double val) {
        matrix.get(i)[I][k] = val;
    }

    final protected void dmx(int i, int k, double val) {
        matrix.get(i)[D][k] = val;
    }

    final protected void xmx(int i, XSTATES s, double val) {
        xmx.get(i)[s.ordinal()] = val;
    }

    public static double scoreSequence(ProfileHMM hmm, Sequence seq) {
        return scoreSequence(hmm, seq.getSeqString());
    }

    public static double scoreSequence(ProfileHMM hmm, String seq) {
        ViterbiScorer scorer = new ViterbiScorer(hmm);
        char[] bases = seq.toCharArray();
        for (char base : bases) {
            scorer.consume(base);
        }
        return scorer.getMaxScore();
    }
}
