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

import java.util.List;
import java.util.Collections;

/**
 *
 * @author fishjord
 */
public class PairwiseAlignment {

    private String alignedSeqi;
    private String alignedSeqj;
    private List<Integer> scores;
    private int starti, endi;
    private int startj, endj;
    private double ident = Double.NaN;

    public PairwiseAlignment(String alignedSeqi, String alignedSeqj, List<Integer> scores, int starti, int endi, int startj, int endj) {
        this.alignedSeqi = alignedSeqi;
        this.alignedSeqj = alignedSeqj;
        this.scores = Collections.unmodifiableList(scores);

        this.starti = starti;
        this.startj = startj;
        this.endi = endi;
        this.endj = endj;
    }

    public String getAlignedSeqi() {
        return alignedSeqi;
    }

    public String getAlignedSeqj() {
        return alignedSeqj;
    }

    public int getScore() {
        return scores.get(scores.size() - 1);
    }

    public List<Integer> getScores() {
        return scores;
    }

    public int getEndi() {
        return endi;
    }

    public int getEndj() {
        return endj;
    }

    public int getStarti() {
        return starti;
    }

    public int getStartj() {
        return startj;
    }
    
    public double getIdent(){
        return ident;
    }
    
    public void setIdent(double i){
        ident = i;
    }
}
