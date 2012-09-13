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
import edu.msu.cme.rdp.alignment.hmm.XSC;
import edu.msu.cme.rdp.alignment.hmm.XSTATES;

import static edu.msu.cme.rdp.alignment.hmm.TSC.*;
import static edu.msu.cme.rdp.alignment.hmm.XSTATES.*;

import static java.lang.StrictMath.*;

/**
 *
 * @author fishjord
 */
public final class ViterbiScorer extends HMMScorer {

    public ViterbiScorer(ProfileHMM hmm) {
        this(hmm, -1);
    }

    public ViterbiScorer(ProfileHMM hmm, int startingState) {
        super(hmm, startingState);
    }

    protected final void initializeSM() {
        for (int k = 0; k <= modelLength; k++) {
            if (k == (startingState - 1) || startingState == -1) {
                mmx(0, k, 0);
                imx(0, k, 0);
                dmx(0, k, 0);
            } else {
                mmx(0, k, Double.NEGATIVE_INFINITY);
                imx(0, k, Double.NEGATIVE_INFINITY);
                dmx(0, k, Double.NEGATIVE_INFINITY);
            }
        }

        xmx(0, N, 0);
        xmx(0, B, hmm.xsc(N, XSC.MOVE));
        xmx(0, XSTATES.E, Double.NEGATIVE_INFINITY);
        xmx(0, C, Double.NEGATIVE_INFINITY);
        xmx(0, J, Double.NEGATIVE_INFINITY);
    }

    final public void consume(char b) {
        extend();
        //seq.append(b);

        double lastMaxScore = Double.NEGATIVE_INFINITY;

        double sc;
        double[] mx = mmx(i);
        double[] ix = imx(i);
        double[] dx = dmx(i);

        double[] mx_m1 = mmx(i - 1);
        double[] ix_m1 = imx(i - 1);
        double[] dx_m1 = dmx(i - 1);
        
        double[] tscMM = hmm.tsc(MM);
        double[] tscIM = hmm.tsc(IM);
        double[] tscDM = hmm.tsc(DM);
        
        double[] tscMI = hmm.tsc(MI);
        double[] tscII = hmm.tsc(II);
        
        double[] tscMD = hmm.tsc(MD);
        double[] tscDD = hmm.tsc(DD);
        
        if(startingState == 0) {
            System.err.println("It's really weird for the starting state to be zero...");
        }
        
        if(startingState > 0) {
            mx[startingState - 1] = ix[startingState - 1] = dx[startingState - 1] = Double.NEGATIVE_INFINITY;
        } else {
            mx[0] = ix[0] = dx[0] = Double.NEGATIVE_INFINITY;            
        }

        //xmx(i, XSTATES.E, Float.NEGATIVE_INFINITY);

        for (int k = max(1, startingState); k < modelLength; k++) {
            sc = max(mx_m1[k - 1] + tscMM[k - 1],
                    ix_m1[k - 1] + tscIM[k - 1]);
            sc = max(sc, dx_m1[k - 1] + tscDM[k - 1]);
            //sc = max(sc, xmx(i - 1, B) + hmm.tsc(k - 1, BM));     //So there is a small problem with this since we don't know the length at the beginning...

            mx[k] = sc + hmm.msc(k, b);

            //xmx[i][XSTATES.E.ordinal()] = max(xmx[i][XSTATES.E.ordinal()], matrix[M][i][k] + esc); //So there is a small problem with this since we don't know the length at the beginning...but since esc is inf for glocal we're cool

            sc = max(mx_m1[k] + tscMI[k],
                    ix_m1[k] + tscII[k]);

            ix[k] = sc + hmm.isc(k, b);

            dx[k] = max(mx[k - 1] + tscMD[k - 1],
                    dx[k - 1] + tscDD[k - 1]);

            lastMaxScore = max(max(mx[k], dx[k]), lastMaxScore);
        }
        sc = max(mx_m1[modelStates] + tscMM[modelStates],
                ix_m1[modelStates] + tscIM[modelStates]);

        sc = max(sc, dx_m1[modelStates] + tscDM[modelStates]);
        //sc = max(sc, xmx(i - 1, B) + hmm.tsc(hmm.M() - 1, BM));   //So there is a small problem with this since we don't know the length at the beginning...

        mx[modelLength] = sc + hmm.msc(modelLength, b);

        ix[modelLength] = Double.NEGATIVE_INFINITY;

        dx[modelLength] = max(mmx(i, modelStates) + tscMD[modelStates],
                dx[modelStates] + tscDD[modelStates]);

        lastMaxScore = max(max(mx[modelStates], dx[modelStates]), lastMaxScore);

        maxScores.add(lastMaxScore);

        //E
        sc = max(xmx(i, XSTATES.E), mmx(i, hmm.M()));
        xmx(i, XSTATES.E, max(sc, dmx(i, hmm.M())));

        /*
         * In theory these should not be used for unihit
         */
        //J
        sc = xmx(i - 1, J) + hmm.xsc(J, XSC.LOOP);
        xmx(i, J, max(sc, xmx(i, XSTATES.E) + hmm.xsc(XSTATES.E, XSC.LOOP)));

        //C
        sc = xmx(i - 1, C) + hmm.xsc(C, XSC.LOOP);
        xmx(i, C, max(sc, xmx(i, XSTATES.E) + hmm.xsc(XSTATES.E, XSC.MOVE)));

        //N
        xmx(i, N, xmx(i - 1, N) + hmm.xsc(N, XSC.LOOP));

        //B
        /*
         * In theory J is always going to be infinite in unihit configuration
         * 
         */
        //sc = xmx(i, N) + hmm.xsc(N, XSC.MOVE);
        xmx(i, B, max(sc, xmx(i, J) + hmm.xsc(J, XSC.MOVE)));
    }
}
