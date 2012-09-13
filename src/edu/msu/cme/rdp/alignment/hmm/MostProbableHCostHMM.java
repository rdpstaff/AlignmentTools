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

import static edu.msu.cme.rdp.alignment.hmm.TSC.*;

/**
 *
 * @author fishjord
 */
public class MostProbableHCostHMM {

    private ProfileHMM hmm;
    private static final char[] mapping = new char[]{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
    private final double[][] mostProbFromState;

    public MostProbableHCostHMM(ProfileHMM hmm) {
        this.hmm = hmm;
        mostProbFromState = new double[128][hmm.M() + 1];

        for (int index = 0; index <= hmm.M(); index++) {
            mostProbFromState['m'][index] = computeCostInternal('m', index);
            mostProbFromState['i'][index] = computeCostInternal('i', index);
            mostProbFromState['d'][index] = computeCostInternal('d', index);
        }
    }

    public double computeHeuristicCost(char state, int stateNo) {
        return mostProbFromState[state][stateNo];
    }

    private double computeCostInternal(char prevState, int stateNo) {

        double h = 0;
        double matchTrans = Double.NEGATIVE_INFINITY;
        double insTrans = Double.NEGATIVE_INFINITY;
        double delTrans = Double.NEGATIVE_INFINITY;
        double bestMatchProb = Double.NEGATIVE_INFINITY;
        double bestInsProb = Double.NEGATIVE_INFINITY;

        for (int state = stateNo + 1; state <= hmm.M(); state++) {
            //System.out.println(state + " " + layerScale);

            switch (prevState) {
                case 'm':
                    matchTrans = hmm.tsc(state - 1, MM);
                    insTrans = hmm.tsc(state - 1, MI);
                    delTrans = hmm.tsc(state - 1, MD);
                    break;
                case 'd':
                    matchTrans = hmm.tsc(state - 1, DM);
                    insTrans = Double.NEGATIVE_INFINITY;
                    delTrans = hmm.tsc(state - 1, DD);
                    break;
                case 'i':
                    matchTrans = hmm.tsc(state - 1, IM);
                    insTrans = hmm.tsc(state - 1, II);
                    delTrans = Double.NEGATIVE_INFINITY;
                    break;
                default:
                    throw new RuntimeException("I hate you.");
            }

            bestMatchProb = Float.NEGATIVE_INFINITY;
            bestInsProb = Float.NEGATIVE_INFINITY;

            for (int k = 0; k < hmm.K(); k++) {
                if (hmm.msc(state, k) > bestMatchProb) {
                    bestMatchProb = hmm.msc(state, k);
                }

                if (hmm.isc(state, k) > bestInsProb) {
                    bestInsProb = hmm.isc(state, k);
                }
            }

            matchTrans += bestMatchProb - hmm.getMaxMatchEmission(state);
            delTrans -= hmm.getMaxMatchEmission(state);
            insTrans += bestInsProb;
	    insTrans = Float.NEGATIVE_INFINITY;

            if (insTrans > matchTrans && insTrans > delTrans) {
                h += insTrans;
                prevState = 'i';
                state--;
            } else if (delTrans > matchTrans && delTrans > insTrans) {
                h += delTrans;
                prevState = 'd';
            } else {
                h += matchTrans;
                prevState = 'm';
            }
        }

        return h;
    }
}
