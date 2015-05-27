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

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.rna.DistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.IdentityDistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.OverlapCheckFailedException;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.utils.kmermatch.KmerMatchCore;
import edu.msu.cme.rdp.readseq.utils.kmermatch.NuclSeqMatch;
import edu.msu.cme.rdp.readseq.utils.kmermatch.ProteinSeqMatch;
import edu.msu.cme.rdp.readseq.utils.orientation.GoodWordIterator;
import edu.msu.cme.rdp.readseq.utils.orientation.ProteinWordGenerator;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class PairwiseKNN {

    private File queryFile;
    private File refFile;
    private int k;  
    private int prefilter = 0;
    private int wordSize;
    private  AlignmentMode mode;
    private List<Sequence> dbSeqs;
    private PrintStream out;
    private static final String dformat = "%1$.3f";
    
    public static class Neighbor {

        PairwiseAlignment alignment;
        boolean reverse;
        Sequence dbSeq;
    }

    private static <T> void insert(T n, List<T> list, Comparator<T> comp, int k) {
        int i = list.size();
        list.add(n);

        while (i > 0 && comp.compare(list.get(i), list.get(i - 1)) > 0) {
            Collections.swap(list, i, i - 1);
            i--;
        }

        if (list.size() > k) {
            list.remove(k);
        }
    }

    public static List<Neighbor> getKNN(Sequence query, List<Sequence> dbSeqs, AlignmentMode mode, int k, int wordSize, int prefilter) throws IOException {
        List<Neighbor> ret = new ArrayList();
        Neighbor n;
        Comparator c = new Comparator<Neighbor>() {
            public int compare(Neighbor t, Neighbor t1) {
                return t.alignment.getScore() - t1.alignment.getScore();
            }
        };

        SequenceType seqType = SeqUtils.guessSequenceType(query);
        ScoringMatrix matrix;
        KmerMatchCore kerMatchCore;
        if (seqType == SequenceType.Nucleotide) {
            matrix = ScoringMatrix.getDefaultNuclMatrix();
            kerMatchCore = new NuclSeqMatch(dbSeqs, wordSize); 
        } else {
            matrix = ScoringMatrix.getDefaultProteinMatrix();
            kerMatchCore = new ProteinSeqMatch(dbSeqs, wordSize); 
        }

        List<Sequence> refList;
       
        if ( prefilter == 0) {   // do not pre-filter the reference seqs
            refList = dbSeqs;
        }else {
            refList = new ArrayList<Sequence>();
            ArrayList<ProteinSeqMatch.BestMatch> topKMatches= kerMatchCore.findTopKMatch(query, prefilter);
            for (KmerMatchCore.BestMatch bestTarget : topKMatches) {
                refList.add(bestTarget.getBestMatch());
            }
        }
        
        for (Sequence dbSeq : refList) {
            n = new Neighbor();
            n.dbSeq = dbSeq;
            PairwiseAlignment fwd = PairwiseAligner.align(n.dbSeq.getSeqString(), query.getSeqString(), matrix, mode);
            if (seqType == SequenceType.Nucleotide) {
                PairwiseAlignment rc = PairwiseAligner.align(n.dbSeq.getSeqString(), IUBUtilities.reverseComplement(query.getSeqString()), matrix, mode);

                if (rc.getScore() > fwd.getScore()) {
                    n.alignment = rc;
                    n.reverse = true;
                } else {
                    n.alignment = fwd;
                    n.reverse = false;
                }
            } else {
                n.alignment = fwd;
                n.reverse = false;
            }

            insert(n, ret, c, k);
        }
                    
        return ret;
    }

    public PairwiseKNN(File queryFile, File refFile, PrintStream out, AlignmentMode mode, int k, int wordSize, int prefilter) throws IOException{
        this.queryFile = queryFile;
        this.refFile = refFile;
        this.out = out;
        this.mode = mode;
        this.k = k;
        this.prefilter = prefilter;
        this.wordSize = wordSize;
        SequenceType querySeqType = SeqUtils.guessSequenceType(queryFile);
        SequenceType refSeqType = SeqUtils.guessSequenceType(refFile);

        if ( querySeqType !=  refSeqType) {
            throw new RuntimeException("reference seqs and query seqs must be the same type, either protein or nucleotide. " );
        }
        if ( wordSize == 0 ){
            if ( refSeqType == SequenceType.Protein){
                this.wordSize = ProteinWordGenerator.WORDSIZE;
            } else {
                this.wordSize = GoodWordIterator.DEFAULT_WORDSIZE ;
            }
        }
            
        dbSeqs = SequenceReader.readFully(refFile);
    }
    
    public void match() throws IOException, OverlapCheckFailedException {
        match(false);
    }
    
    public void match(boolean removeBaseN) throws IOException, OverlapCheckFailedException {
        DistanceModel dist = new IdentityDistanceModel();

        out.println("#query file: " + queryFile.getName() + " db file: " + refFile.getName() + " k: " + k + " mode: " + mode + " usePrefilter: " + prefilter);
        out.println("#seqname\tk\tref seqid\tref desc\torientation\tscore\tident\tquery start\tquery end\tquery length\tref start\tref end");
        Sequence seq;
        List<Neighbor> alignments;
        Neighbor n;
        PairwiseAlignment alignment;
        SequenceReader queryReader = new SequenceReader(queryFile);
        while ((seq = queryReader.readNextSequence()) != null) {
            // remove the N's
            if ( removeBaseN){
                Sequence temp = new Sequence(seq.getSeqName(), seq.getDesc(), seq.getSeqString().toUpperCase().replace("N", ""));
                seq = temp;
            }
            System.err.println("seq " + seq.getSeqString());
            alignments = getKNN(seq, dbSeqs, mode, k, wordSize, prefilter);

            for (int index = 0; index < alignments.size(); index++) {
                n = alignments.get(index);
                alignment = n.alignment;
                double ident = 1 - dist.getDistance(alignment.getAlignedSeqi().getBytes(), alignment.getAlignedSeqj().getBytes(), 0);

                out.println("@" + seq.getSeqName()
                        + "\t" + (index + 1)
                        + "\t" + n.dbSeq.getSeqName()
                        + "\t" + n.dbSeq.getDesc()
                        + "\t" + (n.reverse ? "-" : "+")
                        + "\t" + alignment.getScore()
                        + "\t" + String.format(dformat,ident)
                        + "\t" + alignment.getStartj()
                        + "\t" + alignment.getEndj()
                        + "\t" + seq.getSeqString().length()
                        + "\t" + alignment.getStarti()
                        + "\t" + alignment.getEndi());

                out.println(">" + alignment.getAlignedSeqj());
                out.println(">" + alignment.getAlignedSeqi());
            }
        }
        queryReader.close();
        out.close();
    }
    
    public static void main(String[] args) throws Exception { 
        File queryFile;
        File refFile;
        AlignmentMode mode = AlignmentMode.glocal;
        int k = 1;
        int wordSize = 0 ;
        int prefilter = 10 ;  //  The top p closest protein targets
        boolean removeBaseN = false;
        PrintStream out = new PrintStream(System.out);

        Options options = new Options();
        options.addOption("m", "mode", true, "Alignment mode {global, glocal, local, overlap, overlap_trimmed} (default= glocal)");
        options.addOption("k", true, "K-nearest neighbors to return. (default = 1)");
        options.addOption("o", "out", true, "Redirect output to file instead of stdout");
        options.addOption("p", "prefilter", true, "The top p closest targets from kmer prefilter step. Set p=0 to disable the prefilter step. (default = 10) ");
        options.addOption("w", "word-size", true, "The word size used to find closest targets during prefilter. (default " + ProteinWordGenerator.WORDSIZE 
                + " for protein, " + GoodWordIterator.DEFAULT_WORDSIZE  + " for nucleotide)");
        options.addOption("n", false, "Remove Ns from the query. Default is false");

        try {
            CommandLine line = new PosixParser().parse(options, args);

            if (line.hasOption("mode")) {
                mode = AlignmentMode.valueOf(line.getOptionValue("mode"));
            }

            if (line.hasOption('k')) {
                k = Integer.valueOf(line.getOptionValue('k'));
                if ( k < 1 ){
                    throw new Exception("k must be at least 1");
                }
            }
           
            if (line.hasOption("word-size")) {
                wordSize = Integer.parseInt(line.getOptionValue("word-size"));
                if ( wordSize < 3 ){
                    throw new Exception("Word size must be at least 3");
                }
            }
            if (line.hasOption("prefilter")) {
                prefilter = Integer.parseInt(line.getOptionValue("prefilter"));
                // prefilter == 0 means no prefilter
                if ( prefilter > 0 && prefilter < k ){
                    throw new Exception("prefilter must be at least as big as k " + k);
                }
            }
            
            if (line.hasOption("out")) {
                out = new PrintStream(line.getOptionValue("out"));
            }
            if (line.hasOption('n')) {
                removeBaseN = true;
            }
            args = line.getArgs();

            if (args.length != 2) {
                throw new Exception("Unexpected number of command line arguments");
            }

            queryFile = new File(args[0]);
            refFile = new File(args[1]);            

        } catch (Exception e) {
            new HelpFormatter().printHelp("PairwiseKNN <options> <queryFile> <dbFile>", options);
            System.err.println("ERROR: " + e.getMessage());
            return;
        }
        
        PairwiseKNN theObj = new PairwiseKNN(queryFile, refFile, out, mode, k, wordSize, prefilter);
        theObj.match(removeBaseN);
        

    }
}
