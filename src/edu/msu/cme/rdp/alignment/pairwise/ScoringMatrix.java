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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;

/**
 *
 * @author fishjord
 */
@XmlAccessorType(XmlAccessType.FIELD)
public class ScoringMatrix {

    private static ScoringMatrix defaultNuclMatrix = null;
    private static ScoringMatrix defaultProtMatrix = null;
    private static ScoringMatrix defaultProtMetricMatrix = null;

    
    public static final int DEFAULT_GAP_OPEN_PEALTY = -10;
    public static final int DEFAULT_GAP_EXT_PENALTY = -1;
    public static final int DEFAULT_FRAME_SHIFT_PENALTY = -10;
    public static final int DEFAULT_METRIC_GAP_OPEN_PEALTY = -13;
    public static final int DEFAULT_METRIC_GAP_EXT_PENALTY = -4;
    
    private int[][] scoringMatrix;
    private int[] reverseLookup = new int [127];
    private int gapPenalty;
    private int gapExtend;
    private int frameshiftPenalty = DEFAULT_FRAME_SHIFT_PENALTY;

    private ScoringMatrix() {}

    public ScoringMatrix(File f, int gapOpenPenalty, int gapExtendPenalty) throws IOException {
        this(new FileInputStream(f), gapOpenPenalty, gapExtendPenalty);
    }
    
    public ScoringMatrix(InputStream is, int gapOpenPenalty, int gapExtendPenalty, int frameshiftPenalty) throws IOException {
        this(is, gapOpenPenalty, gapExtendPenalty);
        this.frameshiftPenalty = frameshiftPenalty;
    }
        
    public ScoringMatrix(InputStream is, int gapOpenPenalty, int gapExtendPenalty) throws IOException {
        this.gapPenalty = gapOpenPenalty;
        this.gapExtend = gapExtendPenalty;

        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        String line;

        for(int index = 0;index < 127;index++) {
            reverseLookup[index] = -1;
        }

        int curRow = -1;
        List<Character> validChars = new ArrayList();
        while ((line = reader.readLine()) != null) {
            line = line.trim().toLowerCase();
            if (line.startsWith("#") || line.equals("")) {
                continue;
            }

            if (curRow == -1) {
                for (String s : line.split("\\s+")) {
                    if (s.length() > 1) {
                        throw new IOException("Scoring matrix header " + s + " is longer than 1 character, cannot parse");
                    }
                    validChars.add(s.charAt(0));
                }
                curRow = 0;
                scoringMatrix = new int[validChars.size()][validChars.size()];
            } else {
                String[] values = line.split("\\s+");
                if (values.length != validChars.size() + 1) {
                    throw new IOException("Expected " + validChars.size() + " columns in matrix row " + curRow + " but instead found " + (values.length - 1));
                }

                if (values[0].charAt(0) != validChars.get(curRow) || values[0].length() > 1) {
                    throw new IOException("Scoring matrix must be symetric (expected " + validChars.get(curRow) + " but saw " + values[0] + ")");
                }

                for (int index = 1; index < values.length; index++) {
                    scoringMatrix[curRow][index - 1] = Integer.valueOf(values[index]);
                }
                curRow++;
            }
        }

        for(int index = 0;index < validChars.size();index++) {
            reverseLookup[Character.toLowerCase(validChars.get(index))] = reverseLookup[Character.toUpperCase(validChars.get(index))] = index;
        }
    }

    public static ScoringMatrix getDefaultNuclMatrix() {
        if(defaultNuclMatrix == null) {
            try {
                defaultNuclMatrix = new ScoringMatrix(ScoringMatrix.class.getResourceAsStream("/data/NUC.4.4"), DEFAULT_GAP_OPEN_PEALTY, DEFAULT_GAP_EXT_PENALTY);
            } catch(Exception e) {
                throw new RuntimeException("Failed to get default nucl matrix...something is very wrong!", e);
            }
        }

        return defaultNuclMatrix;
    }

    public static ScoringMatrix getSimpleScoringMatrix(int match, int mismatch) {
        ScoringMatrix ret = new ScoringMatrix();
        ret.gapExtend = mismatch;
        ret.gapPenalty = mismatch;

        for(int index = 0;index < 127;index++) {
            ret.reverseLookup[index] = -1;
        }

        for(char c = 'a';c <= 'z';c++) {
            ret.reverseLookup[c] = ret.reverseLookup[Character.toUpperCase(c)] = c - 'a';
        }

        ret.scoringMatrix = new int['z' - 'a']['z' - 'a'];
        for(int row = 0;row < ret.scoringMatrix.length;row++) {
            ret.scoringMatrix[row][row] = match;
            for(int col = row + 1;col < ret.scoringMatrix.length;col++) {
                ret.scoringMatrix[row][col] = ret.scoringMatrix[col][row] = mismatch;
            }
        }

        return ret;
    }

    public static ScoringMatrix getDefaultProteinMatrix() {
        if(defaultProtMatrix == null) {
            try {
                defaultProtMatrix = new ScoringMatrix(ScoringMatrix.class.getResourceAsStream("/data/blosum62.txt"), DEFAULT_GAP_OPEN_PEALTY, DEFAULT_GAP_EXT_PENALTY);               
            } catch(Exception e) {
                throw new RuntimeException("Failed to get default protein scoring matrix...something is very wrong!", e);
            }
        }

        return defaultProtMatrix;
    }
    
    public static ScoringMatrix getDefaultProteinMetricMatrix() {
        if(defaultProtMetricMatrix == null) {
            try {
                defaultProtMetricMatrix = new ScoringMatrix(ScoringMatrix.class.getResourceAsStream("/data/blosum62_metric.txt"), DEFAULT_METRIC_GAP_OPEN_PEALTY, DEFAULT_METRIC_GAP_EXT_PENALTY);
            } catch(Exception e) {
                throw new RuntimeException("Failed to get default protein scoring matrix blosum62_metric...something is very wrong!", e);
            }
        }

        return defaultProtMetricMatrix;
    }

    public int score(Character b1, Character b2) {
        int i1 = reverseLookup[b1];
        int i2 = reverseLookup[b2];

        if (i1 == -1 || i2 == -1) {
            throw new IllegalArgumentException("Cannot score " + b1 + ", " + b2);
        }

        return scoringMatrix[i1][i2];
    }

    public int getGapOpen() {
        return gapPenalty;
    }

    public int getGapExtend() {
        return gapExtend;
    }
    
    public int getFrameshiftPenalty(){
        return frameshiftPenalty;
    }
    
    public int getIndelPenalty(){
        return gapPenalty + gapExtend ;
    }
    
    public static InputStream getDefaultProteinMatrixStream(){
        return ScoringMatrix.class.getResourceAsStream("/data/blosum62.txt");
    }
    
    public static InputStream getDefaultProteinMatrixMetricStream(){
        return ScoringMatrix.class.getResourceAsStream("/data/blosum62_metric.txt");
    }
}
