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
package edu.msu.cme.rdp.alignment.hmm;

import edu.msu.cme.rdp.readseq.SequenceType;
import java.util.Arrays;

import static java.lang.StrictMath.*;

import static edu.msu.cme.rdp.alignment.hmm.TSC.*;
import static edu.msu.cme.rdp.alignment.hmm.XSC.*;

/**
 *
 * This is a trimmed down representation of a HMMER3/b ascii formatted model
 *
 * It disregards several fields and simply checked to make sure we're using an
 * amino alphabet, checks the format against HMMER3/b and the length
 *
 * Everything else besides the model is ignored
 *
 * @author fishjord
 */
public class ProfileHMM {

    static final int NUM_EMISSION_STATES = 2;
    static final int NUM_TRANSITIONS = 7;
    static final int NXSTATES = 5;
    static final int NXTRANS = 2;
    static final int msc = 0;
    static final int isc = 1;
    String version;
    String name;
    SequenceType alphabet;
    double[][] transitions;
    double[][][] emissions;
    double[] compo;
    double[][] xsc = new double[NXSTATES][NXTRANS];
    int[] alphaMapping;
    double[] maxMatchEmissions;
    Integer m;          //model length
    int k;              //alphabet length
    int l;
    private double nj;
    private boolean multihit = false;
    private boolean local = false;
    private boolean normalized;
    private MostProbableHCostHMM hcost;
    private static final double[] aa_bg;

    static {
        double[] bg = new double[20];

        bg[0] = 0.0787945;             /* A */
        bg[1] = 0.0151600;             /* C */
        bg[2] = 0.0535222;             /* D */
        bg[3] = 0.0668298;             /* E */
        bg[4] = 0.0397062;             /* F */
        bg[5] = 0.0695071;             /* G */
        bg[6] = 0.0229198;             /* H */
        bg[7] = 0.0590092;             /* I */
        bg[8] = 0.0594422;             /* K */
        bg[9] = 0.0963728;             /* L */
        bg[10] = 0.0237718;             /* M */
        bg[11] = 0.0414386;             /* N */
        bg[12] = 0.0482904;             /* P */
        bg[13] = 0.0395639;             /* Q */
        bg[14] = 0.0540978;             /* R */
        bg[15] = 0.0683364;             /* S */
        bg[16] = 0.0540687;             /* T */
        bg[17] = 0.0673417;             /* V */
        bg[18] = 0.0114135;             /* W */
        bg[19] = 0.0304133;             /* Y */


        aa_bg = bg;
    }

    @Override
    public ProfileHMM clone() {
        ProfileHMM copy = new ProfileHMM();
        copy.version = version;
        copy.name = name;
        copy.alphabet = alphabet;
        copy.m = m;
        copy.k = k;
        copy.l = l;
        copy.nj = nj;
        copy.multihit = multihit;
        copy.local = local;

        copy.alphaMapping = Arrays.copyOf(alphaMapping, alphaMapping.length);

        if (compo != null) {
            copy.compo = Arrays.copyOf(compo, compo.length);
        }

        copy.transitions = new double[transitions.length][];
        for (int index = 0; index < transitions.length; index++) {
            copy.transitions[index] = Arrays.copyOf(transitions[index], transitions[index].length);
        }

        copy.xsc = new double[xsc.length][];
        for (int index = 0; index < xsc.length; index++) {
            copy.xsc[index] = Arrays.copyOf(xsc[index], xsc[index].length);
        }

        copy.emissions = new double[emissions.length][][];
        for (int i = 0; i < emissions.length; i++) {
            copy.emissions[i] = new double[emissions[i].length][];
            for (int j = 0; j < emissions[i].length; j++) {
                copy.emissions[i][j] = Arrays.copyOf(emissions[i][j], emissions[i][j].length);
            }
        }

        copy.maxMatchEmissions = Arrays.copyOf(maxMatchEmissions, maxMatchEmissions.length);

        return copy;
    }

    public ProfileHMM() {
    }

    public String getName() {
        return name;
    }

    public int M() {
        return m;
    }

    public int K() {
        return k;
    }

    public double msc(int k, char b) {
        if (alphaMapping[b] == -1) {
            throw new IllegalArgumentException("No mapping for " + b);
        }
        return emissions[k][alphaMapping[b]][msc];
    }

    public void rescaleMatchEmission(int k, char b, double scale) {
        msc(k, alphaMapping[b], msc(k, b) * scale);
    }

    public double isc(int k, char b) {
        return emissions[k][alphaMapping[b]][isc];
    }

    public double msc(int k, int b) {
        if (k == 0) {
            return Double.NEGATIVE_INFINITY;
        }
        return emissions[k][b][msc];
    }

    public double isc(int k, int b) {
        return emissions[k][b][isc];
    }

    public double tsc(int k, TSC trans) {
        return transitions[trans.ordinal()][k];
    }

    public double[] tsc(TSC trans) {
        return transitions[trans.ordinal()];
    }

    public double xsc(XSTATES xstate, XSC trans) {
        return xsc[xstate.ordinal()][trans.ordinal()];
    }

    public double getMaxMatchEmission(int i) {
        if (normalized) {
            return maxMatchEmissions[i];
        } else {
            return 0;
        }
    }

    void msc(int k, int i, double val) {
        emissions[k][i][msc] = val;
        if (val > maxMatchEmissions[k]) {
            maxMatchEmissions[k] = val;
        }
    }

    void isc(int k, int i, double val) {
        emissions[k][i][isc] = val;
    }

    void tsc(int k, TSC trans, double val) {
        transitions[trans.ordinal()][k] = val;
    }

    public void xsc(XSTATES xstate, XSC trans, double val) {
        xsc[xstate.ordinal()][trans.ordinal()] = val;
    }

    void configureGlocal(boolean normalize) {
        this.normalized = normalize;
        double Z = log(tsc(0, MD));

        Arrays.fill(maxMatchEmissions, 0);

        tsc(0, BM, log(1 - tsc(0, MD)));
        //System.out.println("Z: " + Z + " BM: " + tsc(0, BM));

        for (int k = 1; k < m; k++) {
            tsc(k, BM, Z + log(tsc(k, DM)));
            Z += log(tsc(k, DD));
        }
        tsc(m, BM, Double.NEGATIVE_INFINITY);

        if (multihit) {
            throw new RuntimeException("I hate you");
        } else {
            xsc(XSTATES.E, MOVE, 0);
            xsc(XSTATES.E, LOOP, Double.NEGATIVE_INFINITY);
        }

        for (int k = 1; k <= m; k++) {
            for (int x = 0; x < this.k; x++) {
                double val = msc(k, x);
                if (normalize) {
                    //val /= aa_bg[x];
                    val /= compo[x];
                }

                //System.err.println("k" + k + " x" + x + ": " + log(val));

                msc(k, x, log(val));
                if (normalize) {
                    isc(k, x, 0);
                    //isc(k, x, log(isc(k, x) / compo[x]));
                } else {
                    isc(k, x, log(isc(k, x)));
                }
            }
        }

        for (int k = 0; k <= m; k++) {
            for (TSC trans : TSC.values()) {
                if (trans == BM) {
                    continue;
                }
                //if(trans == TSC.MM) {
                //    System.err.println("k" + k + " MM: " +log(tsc(k, trans)));
                //}
                tsc(k, trans, log(tsc(k, trans)));
            }
        }

        for (int x = 0; x < this.k; x++) {
            isc(m, x, Double.NEGATIVE_INFINITY);
        }
    }

    public void reconfigureLength(int L) {
        double ploop, pmove;

        pmove = ((2.0 + nj) / (L + 2 + nj));
        ploop = 1 - pmove;

        xsc(XSTATES.N, LOOP, log(ploop));
        xsc(XSTATES.C, LOOP, log(ploop));
        xsc(XSTATES.J, LOOP, log(ploop));

        xsc(XSTATES.N, MOVE, log(pmove));
        xsc(XSTATES.C, MOVE, log(pmove));
        xsc(XSTATES.J, MOVE, log(pmove));

        this.l = L;
    }

    public SequenceType getAlphabet() {
        return alphabet;
    }

    public MostProbableHCostHMM getHCost() {
        if (hcost == null) {
            hcost = new MostProbableHCostHMM(this);
        }

        return hcost;
    }
}
