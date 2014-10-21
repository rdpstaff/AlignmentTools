package edu.msu.cme.rdp.alignment;

import edu.msu.cme.rdp.readseq.SequenceFormat;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.*;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class AlignmentMerger {

    public Set parseIgnorePositions(File maskFileName) throws IOException {
        //	 read the positions to ignore
        Set<Integer> ignorePositions = new HashSet<Integer>();
        BufferedReader reader = new BufferedReader(new FileReader(maskFileName));
        String line;
        while ((line = reader.readLine()) != null) {
            if (!line.startsWith("#")) {
                ignorePositions.add(Integer.parseInt(line));
            }
        }
        return ignorePositions;
    }

    /**
     *
     * @param alignmentFiles
     * @param writer
     * @param ignoredPositions Set of model positions that should be masked
     * (expected starting at ONE not ZERO)
     * @return Map of files not merged to the reason they weren't
     * @throws IOException
     */
    public static Map<File, String> mergeAlignment(List<File> alignmentFiles, File mergeToFile, Set<Integer> ignoredPositions) throws IOException {
        OutputStream out = new BufferedOutputStream(new FileOutputStream(mergeToFile));
        try {
            return mergeAlignment(alignmentFiles, out, ignoredPositions);
        } finally {
            out.close();
        }
    }

    /**
     *
     * @param alignmentFiles
     * @param writer
     * @param ignoredPositions Set of model positions that should be masked
     * (expected starting at ONE not ZERO)
     * @return Map of files not merged to the reason they weren't
     * @throws IOException
     */
    public static Map<File, String> mergeAlignment(List<File> alignmentFiles, OutputStream os, Set<Integer> ignoredPositions) throws IOException {
        Map<File, String> skippedFiles = new HashMap();
        Map<File, char[]> refSeqMap = getReferenceSeqs(alignmentFiles, ignoredPositions, skippedFiles);
        int[] insertLengths = getInserts(refSeqMap.values());

        FastaWriter mergeStream = new FastaWriter(os);

        String mergedRefSeq = null;
        for (File f : refSeqMap.keySet()) {
            System.err.println("Merging file " + f);
            mergedRefSeq = merge(f, mergeStream, refSeqMap.get(f), insertLengths);
        }

        mergeStream.writeSeq("#=GC_RF", mergedRefSeq);
        mergeStream.close();

        return skippedFiles;
    }

    private static Map<File, char[]> getReferenceSeqs(List<File> alignmentFiles, Set<Integer> ignoredPositions, Map<File, String> skippedFiles) throws IOException {
        Map<File, char[]> refSeqMap = new LinkedHashMap();
        int realModelLength = -1;

        for (File alignmentFile : alignmentFiles) {
            if (alignmentFile.isFile()) {
                SequenceFormat format = SeqUtils.guessFileFormat(alignmentFile);

                String refSeqId;
                if (format == SequenceFormat.FASTA) {
                    refSeqId = "#=GC_RF";
                } else if (format == SequenceFormat.STK) {
                    refSeqId = "#=GC RF";
                } else {
                    skippedFiles.put(alignmentFile, "Unprocessable sequence format " + format);
                    continue;
                }

                char[] refSeq = null;

                try {
                    if (format == SequenceFormat.FASTA) {
                        BufferedInputStream is = new BufferedInputStream(new FileInputStream(alignmentFile));
                       
                        SeqReader reader = new SequenceReader(is);

                        Sequence seq;
                        while ((seq = reader.readNextSequence()) != null) {
                            if (seq.getSeqName().equals(refSeqId)) {
                                refSeq = seq.getSeqString().toCharArray();
                                break;
                            }
                        }
                        reader.close();
                    } else {
                        IndexedSeqReader reader = new IndexedSeqReader(alignmentFile);
                        refSeq = reader.readSeq(refSeqId).getSeqString().toCharArray();
                        reader.close();
                    }

                    if (refSeq == null) {
                        throw new IOException();
                    }
                    for (int index = 0; index < refSeq.length; index++) {
                        if (refSeq[index] == '~') {
                            refSeq[index] = '.';
                        }
                    }

                } catch (IOException ignore) {
                    skippedFiles.put(alignmentFile, "Reference sequence \"" + refSeqId + "\" not found");
                    continue;
                }

                int modelLength = 0;
                for (int index = 0; index < refSeq.length; index++) {
                    if (refSeq[index] != '.') {
                        modelLength++;
                    }
                    //Increase the model position before we check if it's ignored
                    //since the model positions start at ONE
                    if (ignoredPositions.contains(modelLength)) {
                        refSeq[index] = '.';
                    }
                }

                if (realModelLength == -1) {
                    realModelLength = modelLength;
                } else if (realModelLength != modelLength) {
                    skippedFiles.put(alignmentFile, "Model length [" + modelLength + "] doesn't match expected [" + realModelLength + "]");
                    continue;
                }

                refSeqMap.put(alignmentFile, refSeq);
            } else {
                skippedFiles.put(alignmentFile, "Not a file");
            }
        }

        return refSeqMap;
    }

    private static int[] getInserts(Collection<char[]> refSeqs) throws IOException {
        //We keep track of the maximum insert length at every model position
        //The convention is the number of inserts AFTER the model position
        //With model positions starting at 1
        //So inserts before the first model position are in insertLengths[0]
        //and inserts after the last are at insertLengths[insertLengths.length - 1]
        int[] insertLengths = null;

        for (char[] refSeq : refSeqs) {
            if (insertLengths == null) {
                int modelLength = new String(refSeq).replace(".", "").length();
                //We can have inserts at the begining AND end
                insertLengths = new int[modelLength + 1];
            }

            int modelPosition = 0;
            int insertLength = 0;
            for (int index = 0; index < refSeq.length; index++) {
                if (refSeq[index] == '.') {
                    insertLength++;
                } else {
                    if (insertLength > insertLengths[modelPosition]) {
                        insertLengths[modelPosition] = insertLength;
                    }
                    modelPosition++;
                    insertLength = 0;
                }
            }

            //Since we won't always have a model position following the last
            //insert, gotta do it like this
            if (insertLength > insertLengths[modelPosition]) {
                insertLengths[modelPosition] = insertLength;
            }
        }

        return insertLengths;
    }

    /**
     *
     * @param alignmentFile
     * @param mergeStream
     * @param refSeq
     * @param insertLengths
     * @return Merged reference sequence
     * @throws IOException
     */
    private static String merge(File alignmentFile, FastaWriter mergeStream, char[] refSeq, final int[] insertLengths) throws IOException {
        SequenceReader reader = new SequenceReader(alignmentFile);
        Sequence seq;
        StringBuilder mergedRef = new StringBuilder();

        String[] paddings = new String[insertLengths.length];

        int modelPosition = 0;
        int insertLength = 0;
        for (int index = 0; index < refSeq.length; index++) {
            if (refSeq[index] == '.') {
                insertLength++;
                mergedRef.append(".");
            } else {
                StringBuilder insertPadding = new StringBuilder();
                for (int ins = insertLength; ins < insertLengths[modelPosition]; ins++) {
                    insertPadding.append('.');
                }

                paddings[modelPosition] = insertPadding.toString();
                mergedRef.append(insertPadding);
                mergedRef.append("x");

                modelPosition++;
                insertLength = 0;
            }
        }

        StringBuilder insertPadding = new StringBuilder();
        for (int ins = insertLength; ins < insertLengths[modelPosition]; ins++) {
            insertPadding.append('.');
        }
        paddings[modelPosition] = insertPadding.toString();
        mergedRef.append(insertPadding);

        while ((seq = reader.readNextSequence()) != null) {
            //Skip meta seqs, we'll add our own dummy ones at the end
            if (seq.getSeqName().startsWith("#")) {
                continue;
            }

            StringBuilder modifiedSeq = new StringBuilder();
            char[] bases = seq.getSeqString().toCharArray();
            modelPosition = 0;

            if (bases.length != refSeq.length) {
                throw new IOException("Sequence " + seq.getSeqName() + " in file " + alignmentFile.getName() + "'s length [" + bases.length + "] isn't the same as it's ref seq length [" + refSeq.length + "]");
            }

            for (int b = 0; b < bases.length; b++) {
                char base = bases[b];
                if (base == 'u' || base == 'U') {
                    base = 't';
                }

                if (refSeq[b] == '.') {
                    modifiedSeq.append(Character.toLowerCase(base));
                } else {
                    modifiedSeq.append(paddings[modelPosition]).append(Character.toUpperCase(base));
                    modelPosition++;
                }
            }

            modifiedSeq.append(paddings[modelPosition]);

            mergeStream.writeSeq(seq.getSeqName(), seq.getDesc(), modifiedSeq.toString());
        }

        reader.close();
        return mergedRef.toString();
    }

    public static void main(String[] args) throws Exception {

        String usage = "Usage: java AlignmentMerger <alignfiledir> <outfile.fasta> [ <mask_file>] \n"
                + "  This program reads in all the files from the input directory and merges the alignment into one single file\n"
                + "  stkdir contains a list of aligned stk files to be merged \n"
                + "  maskfile contains the model positions to be ignored";

        if (args.length != 2 && args.length != 3) {
            throw new IllegalArgumentException(usage);
        }

        String stkdir = args[0];
        String outfile = args[1];

        AlignmentMerger parser = new AlignmentMerger();
        Set ignorePositions = new HashSet();
        if (args.length == 3) {
            ignorePositions = parser.parseIgnorePositions(new File(args[2]));
        }
        Map<File, String> notMerged = AlignmentMerger.mergeAlignment(Arrays.asList(new File(stkdir).listFiles()), new File(outfile), ignorePositions);
        for (File f : notMerged.keySet()) {
            System.err.println("Did not merge " + f + " beacse " + notMerged.get(f));
        }
    }
}
