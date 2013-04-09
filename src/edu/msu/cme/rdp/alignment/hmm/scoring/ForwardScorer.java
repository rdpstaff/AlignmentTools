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
public class ForwardScorer extends HMMScorer {

    private static final double[] logsumLookup = new double[16000];

    static {
        for (int i = 0; i < 16000; i++) {
            logsumLookup[i] = log(1 + exp((double) i / 1000));
        }
    }

    private static final double logsum(double a, double b) {
        double max = max(a, b);
        double min = min(a, b);

        if (min == Double.NEGATIVE_INFINITY) {
            return max;
        }

        return max + log(1 + exp(min - max));

        //return (min == Double.NEGATIVE_INFINITY || (max - min) >= 15.7) ? max : max + logsumLookup[(int) ((max - min) * 1000)];
    }

    public ForwardScorer(ProfileHMM hmm) {
        this(hmm, -1);
    }

    public ForwardScorer(ProfileHMM hmm, int startingState) {
        super(hmm, startingState);
    }

    protected final void initializeSM() {
        for (int k = 0; k <= hmm.M(); k++) {
            if (k == startingState || startingState == -1) {
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

        double sc;

        mmx(i, 0, Double.NEGATIVE_INFINITY);
        imx(i, 0, Double.NEGATIVE_INFINITY);
        dmx(i, 0, Double.NEGATIVE_INFINITY);

        xmx(i, XSTATES.E, Double.NEGATIVE_INFINITY);

        double lastMaxScore;
        if (i == 0) {
            lastMaxScore = Double.NEGATIVE_INFINITY;
        } else {
            lastMaxScore = maxScores.get(i - 1);
        }

        for (int k = 1; k < hmm.M(); k++) {
            sc = logsum(
                    logsum(
                    mmx(i - 1, k - 1) + hmm.tsc(k - 1, MM),
                    imx(i - 1, k - 1) + hmm.tsc(k - 1, IM)),
                    dmx(i - 1, k - 1) + hmm.tsc(k - 1, DM) //Simplified by removing an extra logsum(xmx(i-1, begin state) + tsc(BM)) which will always be neg infinity, we're doing begining transitions differently
                    );

            mmx(i, k, sc + hmm.msc(k, b));

            sc = logsum(
                    mmx(i - 1, k) + hmm.tsc(k, MI),
                    imx(i - 1, k) + hmm.tsc(k, II));

            imx(i, k, sc + hmm.isc(k, b));

            sc = logsum(
                    mmx(i, k - 1) + hmm.tsc(k - 1, MD),
                    dmx(i, k - 1) + hmm.tsc(k - 1, DD));

            dmx(i, k, sc);
            lastMaxScore = max(max(mmx(i, k), dmx(i, k)), lastMaxScore);
        }

        sc = logsum(
                logsum(
                mmx(i - 1, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, MM),
                imx(i - 1, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, IM)),
                dmx(i - 1, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, DM) //Simplified by removing an extra logsum(xmx(i-1, begin state) + tsc(BM)) which will always be neg infinity, we're doing begining transitions differently
                );

        mmx(i, hmm.M(), sc + hmm.msc(hmm.M(), b));

        imx(i, hmm.M(), Double.NEGATIVE_INFINITY);

        sc = logsum(
                mmx(i, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, MD),
                dmx(i, hmm.M() - 1) + hmm.tsc(hmm.M() - 1, DD));

        dmx(i, hmm.M(), sc);

        lastMaxScore = max(max(mmx(i, hmm.M()), dmx(i, hmm.M())), lastMaxScore);

        maxScores.add(lastMaxScore);

        /*
         * Everything below here hasn't been updated for forward...but it isn't
         * used right now so get over it
         */

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
        sc = xmx(i, N) + hmm.xsc(N, XSC.MOVE);
        xmx(i, B, max(sc, xmx(i, J) + hmm.xsc(J, XSC.MOVE)));
    }
}
