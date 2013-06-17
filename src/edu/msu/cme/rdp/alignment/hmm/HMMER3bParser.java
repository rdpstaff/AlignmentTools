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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;

import static java.lang.StrictMath.*;

/**
 *
 * @author fishjord
 */
public class HMMER3bParser {

    public static ProfileHMM readModel(File f) throws IOException {
        return readModel(new FileInputStream(f));
    }

    public static ProfileHMM readModel(InputStream is) throws IOException {
        ProfileHMM hmm = readModelInternal(is);

        hmm.configureGlocal(true);

        return hmm;
    }

    public static ProfileHMM readUnnormalized(File f) throws IOException {
        return readUnnormalized(new FileInputStream(f));
    }

    public static ProfileHMM readUnnormalized(InputStream is) throws IOException {
        ProfileHMM hmm = readModelInternal(is);

        hmm.configureGlocal(false);

        return hmm;
    }

    private static ProfileHMM readModelInternal(InputStream is) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        ProfileHMM hmm = new ProfileHMM();

        hmm.version = reader.readLine().trim().split("\\s+")[0];

        if (!hmm.version.equals("HMMER3/b") && !hmm.version.equals("HMMER3/f")) {
            throw new IOException("Cannot parse " + hmm.version + " version hmmer model, only supports HMMER3/b");
        }

        String line;
        String alpha = null;
        while ((line = reader.readLine()) != null) {
            String[] lexemes = line.trim().split("\\s+");
            if (lexemes[0].equals("NAME")) {
                hmm.name = lexemes[1];
            } else if (lexemes[0].equals("LENG")) {
                hmm.m = Integer.valueOf(lexemes[1]);
            } else if (lexemes[0].equals("ALPH")) {
                alpha = lexemes[1];

                if (alpha.equals("amino")) {
                    hmm.alphabet = SequenceType.Protein;
                } else if (alpha.equals("RNA")) {
                    hmm.alphabet = SequenceType.Nucleotide;
                } else if (alpha.equals("DNA")) {
                    hmm.alphabet = SequenceType.Nucleotide;
                } else {
                    throw new IOException("Unknown sequence alphabet \"" + alpha + "\"");
                }
            } else if (lexemes[0].equals("HMM")) {
                parseAlpha(hmm, line);
                break;
            }
        }

        if (hmm.m == null) {
            throw new IOException("Model length not found");
        }

        if (alpha == null) {
            throw new IOException("No alphabet specified in HMM");
        }

        reader.readLine(); // discard the transition labels (we expect them to be in the order m->m     m->i     m->d     i->m     i->i     d->m     d->d

        hmm.transitions = new double[ProfileHMM.NUM_TRANSITIONS + 1][hmm.m + 1]; // Leave room for the implicit B->M transition that is computed on load
        hmm.emissions = new double[hmm.m + 1][hmm.k][ProfileHMM.NUM_EMISSION_STATES];

        //So we're gonna configure this here...but the model sets these
        //in the msc(state, alpha, val) (checking if val > maxEmission[state])
        //Two biproducts of this, first if the probabilities are tweaked at runtime
        //maxEmissions will be updated accordingly, second, you can't use maxEmissions
        //til the model is calibrated, which you couldn't anyway so it doesn't matter
        hmm.maxMatchEmissions = new double[hmm.m + 1];
        Arrays.fill(hmm.maxMatchEmissions, Double.NEGATIVE_INFINITY);

        line = reader.readLine().trim();
        if (line.startsWith("COMPO")) {
            hmm.compo = readProbabilities(line, 1, hmm.k);
        }

        double[] probabilities;

        for (int index = 0; index <= hmm.m; index++) {
            if (index > 0) {
                probabilities = readProbabilities(reader.readLine(), 1, hmm.k);
                for (int emission = 0; emission < hmm.k; emission++) {
                    hmm.emissions[index][emission][ProfileHMM.msc] = probabilities[emission];
                }
            }

            probabilities = readProbabilities(reader.readLine(), 0, hmm.k);
            for (int emission = 0; emission < hmm.k; emission++) {
                hmm.emissions[index][emission][ProfileHMM.isc] = probabilities[emission];
            }

            probabilities = readProbabilities(reader.readLine(), 0, ProfileHMM.NUM_TRANSITIONS);
            for (int trans = 0; trans < ProfileHMM.NUM_TRANSITIONS; trans++) {
                hmm.transitions[trans][index] = probabilities[trans];
            }
        }

        reader.close();

        return hmm;
    }

    private static void parseAlpha(ProfileHMM hmm, String line) throws IOException {
        String[] lexemes = line.trim().split("\\s+");

        hmm.alphaMapping = new int[127];
        Arrays.fill(hmm.alphaMapping, -1);

        for (int index = 1; index < lexemes.length; index++) {
            if (lexemes[index].length() != 1) {
                throw new IOException("Alphabet symbol " + lexemes[index] + " too long");
            }

            char sym = Character.toUpperCase(lexemes[index].charAt(0));
            hmm.alphaMapping[(int) sym] = index - 1;
            hmm.alphaMapping[(int) Character.toLowerCase(sym)] = index - 1;

            if (sym == 'U' && hmm.alphabet == SequenceType.Nucleotide) {
                hmm.alphaMapping[(int) 'T'] = index - 1;
                hmm.alphaMapping[(int) 't'] = index - 1;
            }
        }

        hmm.k = lexemes.length - 1;
    }

    /**
     * Parses the probailities for a model from the provided string
     *
     * @param line
     * @param offset
     * @param expectedValues
     * @return
     */
    private static double[] readProbabilities(String line, int offset, int expectedValues) throws IOException {
        if (line == null || line.equals("//")) {
            throw new IOException("Unexpected end of model");
        }
        String[] lexemes = line.trim().split("\\s+");
        double[] ret = new double[expectedValues];

        if (lexemes.length < expectedValues + offset) {
            throw new IOException("Too few probabilities on line \"" + line + "\"");
        }

        for (int index = offset; index < offset + expectedValues; index++) {
            if (lexemes[index].equals("*")) {
                ret[index - offset] = 0;
            } else {
                double v = Double.valueOf(lexemes[index]);
                ret[index - offset] = exp(-1 * v);
            }
        }

        return ret;
    }

    public static void main(String[] args) throws IOException {
        ProfileHMM hmm = HMMER3bParser.readModel(new File("/work/fishjord/other_projects/hmmgs/models/rplb/rplb_for.hmm"));

        System.out.print("double expected_compo[] = {");
        for (int index = 0; index < hmm.compo.length; index++) {
            System.out.print(hmm.compo[index]);
            if (index + 1 != hmm.compo.length) {
                System.out.print(",");
            }
        }
        System.out.println("};");

        System.out.println("double expected_msc[][20] = {");
        for (int state = 0; state < hmm.M() + 1; state++) {
            System.out.print("\t{");
            for (int index = 0; index < hmm.compo.length; index++) {
                if (hmm.msc(state, index) == Double.NEGATIVE_INFINITY) {
                    System.out.print("-std::numeric_limits<double>::infinity()");
                } else {
                    System.out.print(hmm.msc(state, index));
                }
                if (index + 1 != hmm.compo.length) {
                    System.out.print(",");
                }
            }
            System.out.print("}");
            if (state + 1 != hmm.M() + 1) {
                System.out.println(",");
            }
        }
        System.out.println("};");

        System.out.println("double expected_tsc[][7] = {");
        for (int state = 0; state < hmm.M() + 1; state++) {
            System.out.print("\t{");
            for (TSC tsc : new TSC[]{TSC.MM, TSC.MI, TSC.MD, TSC.IM, TSC.II, TSC.DM, TSC.DD}) {

                if (hmm.tsc(state, tsc) == Double.NEGATIVE_INFINITY) {
                    System.out.print("-std::numeric_limits<double>::infinity()");
                } else {
                    System.out.print(hmm.tsc(state, tsc));
                }
                if (tsc != TSC.DD) {
                    System.out.print(",");
                }
            }
            System.out.print("}");
            if (state + 1 != hmm.M() + 1) {
                System.out.println(",");
            }
        }
        System.out.println("};");

        System.out.print("double expected_max_emission[] = {");
        for(int state = 0;state <= hmm.M();state++) {
            System.out.print(hmm.getMaxMatchEmission(state));
            if(state != hmm.M()) {
                System.out.println(", ");
            }
        }
        System.out.println("};");

        System.out.println("double expected_m_hcost = " + hmm.getHCost().computeHeuristicCost('m', 0) + ";");
        System.out.println("double expected_d_hcost = " + hmm.getHCost().computeHeuristicCost('d', 0) + ";");
        System.out.println("double expected_i_hcost = " + hmm.getHCost().computeHeuristicCost('i', 0) + ";");
    }
}
