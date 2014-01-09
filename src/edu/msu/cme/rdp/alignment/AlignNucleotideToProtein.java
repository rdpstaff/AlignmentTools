/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.alignment;

import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author fishjord
 */
public class AlignNucleotideToProtein {

    public static void main(String[] args) throws Exception {

        if (args.length != 4) {
            System.err.println("USAGE: AlignNucleotideToProtein <aligned prot seqs> <unaligned_nucl_seqs> <aligned nucl out> <stats out>");
            return;
        }

        File alignProtSeqFile = new File(args[0]);
        File nuclSeqFile = new File(args[1]);
        File alignedNuclOutFile = new File(args[2]);
        File statsOut = new File(args[3]);

        Map<String, Sequence> protTemplateMap = new HashMap();
        for (Sequence seq : SequenceReader.readFully(alignProtSeqFile)) {
            protTemplateMap.put(seq.getSeqName(), seq);
        }

        SequenceReader nuclSeqReader = new SequenceReader(nuclSeqFile);
        Sequence seq;

        PrintStream out = new PrintStream(statsOut);
        FastaWriter alignedNuclOut = new FastaWriter(alignedNuclOutFile);
        boolean seqHasStop = false;
        while ((seq = nuclSeqReader.readNextSequence()) != null) {
            Sequence protTemplate = protTemplateMap.get(seq.getSeqName());

            if (protTemplate == null) {
                System.err.println("Failed to find template for " + seq.getSeqName());
                continue;
            }

            String alignedProtSeq = protTemplate.getSeqString();
            String unalignedProtSeq = SeqUtils.getUnalignedSeqString(alignedProtSeq);

            String unalignedNuclSeq = seq.getSeqString();
            //Just let it truncate, for whatever reason these aren't always aligned on a 3 base boundary in genbank, but they're in frame (that I've seen)
            //Note this does assume that it is the trailing characters that are not translated, not the leading
            int nuclLengthInProt = unalignedNuclSeq.length() / 3;

            //We want to check to see if this is a recoverable error, ie
            //if the nucleotide sequence has a stop codon in the last place
            //but the protein sequence does not (many tools will remove/omit the stop codon, specifically genbank translations leave out the *)
            if (nuclLengthInProt - 1 == unalignedProtSeq.length()) {
                String lastCodon = unalignedNuclSeq.substring(unalignedNuclSeq.length() - 3);
                char lastAA = ProteinUtils.getInstance().translateToProtein(lastCodon, false, 11).charAt(0);
                if (lastAA == '*' && !unalignedProtSeq.endsWith("*")) {
                    unalignedNuclSeq = unalignedNuclSeq.substring(0, unalignedNuclSeq.length() - 3);
                    nuclLengthInProt = unalignedNuclSeq.length() / 3;
                }
            }

            if (unalignedProtSeq.length() != nuclLengthInProt) {
                System.err.println("Length mismatch between " + seq.getSeqName() + " protein [" + unalignedProtSeq.length() + "] and nucleotide sequences [" + unalignedNuclSeq.length() + "], not translating");
                continue;
            }

            int starIndex = unalignedProtSeq.indexOf('*');
            if (starIndex != -1 && starIndex != unalignedProtSeq.length() - 1) {
                System.err.println("Sequence " + seq.getSeqName() + " has a stop codon in that isn't at the end, this could be very bad");
            }

            if (starIndex != -1) {
                seqHasStop = true;
            }

            String alignedNuclSeqString = ProteinUtils.getInstance().getAlignedNucSeq(alignedProtSeq, unalignedNuclSeq, true, 11);
            String sanityCheckSeqString = SeqUtils.getUnalignedSeqString(alignedNuclSeqString);

            if (!sanityCheckSeqString.equalsIgnoreCase(unalignedNuclSeq)) {
                System.err.println(sanityCheckSeqString);
                System.err.println(unalignedNuclSeq);
                System.err.println(alignedNuclSeqString);
                System.err.println("Huge error with " + seq.getSeqName() + ", aligned nucleotide sequence doesn't match the unaligned one");
                continue;
            }

            float translScore = ProteinUtils.getInstance().getTranslScore(unalignedProtSeq, unalignedNuclSeq, true, 11);

            out.println("@ " + seq.getSeqName() + "\t" + translScore);
            for (char c : unalignedProtSeq.toCharArray()) {
                out.print(" " + c + "  | ");
            }
            out.println();

            for (int index = 0; index < unalignedNuclSeq.length(); index += 3) {
                int end = index + 3;
                if (end > seq.getSeqString().length()) {
                    end = seq.getSeqString().length();
                }
                out.print(seq.getSeqString().subSequence(index, end) + " | ");
            }

            out.println();

            String predictedProt = ProteinUtils.getInstance().translateToProtein(unalignedNuclSeq, true, 11);
            for (char c : predictedProt.toCharArray()) {
                out.print(" " + c + "  | ");
            }
            out.println();
            alignedNuclOut.writeSeq(seq.getSeqName(), seq.getDesc(), alignedNuclSeqString);
        }

        if (protTemplateMap.containsKey("#=GC_RF")) {
            StringBuilder mask = new StringBuilder();

            for (char c : protTemplateMap.get("#=GC_RF").getSeqString().toCharArray()) {
                if (c == 'x') {
                    mask.append("xxx");
                } else if (c == '.') {
                    mask.append("...");
                } else {
                    throw new IOException("I hate you. " + c);
                }
            }

            alignedNuclOut.writeSeq("#=GC_RF", mask.toString());
        }

        out.close();
        alignedNuclOut.close();

        if (seqHasStop) {
            System.err.println("One or more protein sequences contained a stop codon, this could cause trouble if HMMER ever aligns a stop codon to a model position");
        }
    }
}
