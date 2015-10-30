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
import edu.msu.cme.rdp.readseq.utils.orientation.OrientationChecker;
import edu.msu.cme.rdp.readseq.utils.orientation.ProteinWordGenerator;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class PairwiseKNN {

    private final File refFile;
    private final int k;  // final number of hits to return
    private final int prefilter;  // number of hits to keep from prefilter stages
    private int wordSize;
    private final AlignmentMode mode;
    private final HashMap<String, Sequence> dbSeqsMap = new HashMap(); // keep all the refseq in memory for pairwise alignment
    private final ScoringMatrix matrix;
    private KmerMatchCore kerMatchCore;
    private static final String dformat = "%1$.3f";
    private static final DistanceModel dist = new IdentityDistanceModel();
    private static final Comparator c = new ResultComparator();
    private final SequenceType refSeqType ;
    
    public static class Neighbor {

        PairwiseAlignment alignment;
        boolean reverse;
        Sequence dbSeq;
        
        public boolean isReverse(){
            return reverse;
        }
        
        public PairwiseAlignment getAlignment(){
            return alignment;
        }
        
        public Sequence getDbSeq(){
            return dbSeq;
        }
    }
    
    public static class ResultComparator implements Comparator<Neighbor> {
        public int compare(Neighbor t, Neighbor t1) {
            return t.alignment.getScore() - t1.alignment.getScore();
        }
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

   
    public PairwiseKNN( File refFile, AlignmentMode mode, int k, int ws, int prefilter) throws IOException{
        this.refFile = refFile;
        this.mode = mode;
        this.k = k;
        this.prefilter = prefilter;
        this.wordSize = ws;
        
        refSeqType = SeqUtils.guessSequenceType(refFile);
        parseRefSeq(refFile);
        if ( refSeqType == SequenceType.Protein){
            matrix = ScoringMatrix.getDefaultProteinMatrix();
            if ( wordSize == 0 ){
                this.wordSize = ProteinWordGenerator.WORDSIZE;
            }
            if ( prefilter > 0){
                kerMatchCore = new ProteinSeqMatch(new ArrayList<Sequence>(dbSeqsMap.values()), wordSize);
            }
        } else {
            matrix = ScoringMatrix.getDefaultNuclMatrix();
            if ( wordSize == 0 ){
                this.wordSize = GoodWordIterator.DEFAULT_WORDSIZE ;
            }
            if ( prefilter > 0){
                kerMatchCore = new NuclSeqMatch(new ArrayList<Sequence>(dbSeqsMap.values()), wordSize);
            }
        }
    }
     
    private synchronized void parseRefSeq(File file) throws IOException{
        SeqReader reader = new SequenceReader(file);
        Sequence seq;
        while ((seq = reader.readNextSequence()) != null) {
            dbSeqsMap.put(seq.getSeqName(), seq);
        }
        reader.close();
    }
    
    public Sequence getRefSeq(String seqName){
        return this.dbSeqsMap.get(seqName);
    }
    
    public String getRefFilename(){
        return this.refFile.getName();
    }
    
    public int getK(){
        return k;
    }
    /**
     * 
     * @param seq
     * @param refList, allow different reference set for flexibility
     * @param isSeqReversed indicates the orientation of the sequence compared to the original seq
     * @param checkReverse for bacteria and archaea, we can use OrientationChecker, for other sequences, we should provide option
     * @return
     * @throws IOException
     * @throws OverlapCheckFailedException 
     */
    public List<Neighbor> getKNN(Sequence seq, Collection<Sequence> refList, boolean removeBaseN, boolean isSeqReversed, boolean checkReverse) throws IOException, OverlapCheckFailedException {
        List<Neighbor> ret = new ArrayList();
        Neighbor n;
        if ( removeBaseN){
            Sequence temp = new Sequence(seq.getSeqName(), seq.getDesc(), seq.getSeqString().toUpperCase().replace("N", ""));
            seq = temp;
        }
        
        for (Sequence dbSeq : refList) {
            n = new Neighbor();
            n.dbSeq = dbSeq;
            PairwiseAlignment fwd = PairwiseAligner.align(n.dbSeq.getSeqString(), seq.getSeqString(), matrix, mode);
            if (refSeqType == SequenceType.Nucleotide &&  checkReverse) {
                            
                PairwiseAlignment rc = PairwiseAligner.align(n.dbSeq.getSeqString(), IUBUtilities.reverseComplement(seq.getSeqString()), matrix, mode);
                if (rc.getScore() > fwd.getScore()) {
                    n.alignment = rc;
                    n.reverse = isSeqReversed ? false :true;
                    double ident = 1 - dist.getDistance(rc.getAlignedSeqi().getBytes(), rc.getAlignedSeqj().getBytes(), 0);
                    rc.setIdent(ident);
                } else {
                    n.alignment = fwd;
                    n.reverse = isSeqReversed ? true : false;
                    double ident = 1 - dist.getDistance(fwd.getAlignedSeqi().getBytes(), fwd.getAlignedSeqj().getBytes(), 0);
                    fwd.setIdent(ident);
                }
               
            } else {
                double ident = 1 - dist.getDistance(fwd.getAlignedSeqi().getBytes(), fwd.getAlignedSeqj().getBytes(), 0);
                fwd.setIdent(ident);
                n.alignment = fwd;
                n.reverse = isSeqReversed;
            }

            insert(n, ret, c, k);
        }
                    
        return ret;
    }  
    
    public List<Neighbor> findMatch(Sequence seq, boolean removeBaseN) throws IOException, OverlapCheckFailedException {
        boolean isReversed = false;
        if ( this.refSeqType == SequenceType.Nucleotide ){
            //check orientation
            isReversed = OrientationChecker.getChecker().isSeqReversed(seq.getSeqString());
            if ( isReversed ){
                seq = new Sequence(seq.getSeqName(), seq.getDesc(), IUBUtilities.reverseComplement(seq.getSeqString()));
            }
        }
        if ( prefilter == 0) {   // do not pre-filter the reference seqs, so need to check both orientation
            return getKNN(seq, dbSeqsMap.values(), removeBaseN, isReversed, true);
        }else {
            List<Sequence> refList = new ArrayList<Sequence>();
            ArrayList<ProteinSeqMatch.BestMatch> topKMatches= kerMatchCore.findTopKMatch(seq, prefilter);
                        
            for (KmerMatchCore.BestMatch bestTarget : topKMatches) {
                refList.add(bestTarget.getBestMatch());
            }
            return getKNN(seq, refList, removeBaseN, isReversed, false);
        }        
    }
    
    private synchronized void printAlignment(Sequence seq, List<Neighbor> alignments, PrintStream out) throws IOException{
        Neighbor n;
        PairwiseAlignment alignment;
        for (int index = 0; index < alignments.size(); index++) {
            n = alignments.get(index);
            alignment = n.alignment;

            out.println("@" + seq.getSeqName()
                    + "\t" + (index + 1)
                    + "\t" + (n.reverse ? "-" : "+")
                    + "\t" + alignment.getScore()
                    + "\t" + String.format(dformat,alignment.getIdent())
                    + "\t" + alignment.getStartj()
                    + "\t" + alignment.getEndj()
                    + "\t" + seq.getSeqString().length()
                    + "\t" + alignment.getStarti()
                    + "\t" + alignment.getEndi()
                    + "\t" + n.dbSeq.getSeqName()
                    + "\t" + n.dbSeq.getDesc());

            out.println(">" + alignment.getAlignedSeqj());
            out.println(">" + alignment.getAlignedSeqi());
        }
    }
   
    public static void main(String[] args) throws Exception { 
        final int maxThreads;
        final int maxTasks = 1000;
        File queryFile;
        File refFile;
        AlignmentMode mode = AlignmentMode.glocal;
        int k = 1;
        int wordSize = 0 ;
        int prefilter = 10 ;  //  The top p closest protein targets
        final boolean removeBaseN;
        final PrintStream out ;

        Options options = new Options();
        options.addOption("m", "mode", true, "Alignment mode {global, glocal, local, overlap, overlap_trim} (default= glocal)");
        options.addOption("k", true, "K-nearest neighbors to return. (default = 1)");
        options.addOption("o", "out", true, "Redirect output to file instead of stdout");
        options.addOption("p", "prefilter", true, "The top p closest targets from kmer prefilter step. Set p=0 to disable the prefilter step. (default = 10) ");
        options.addOption("w", "word-size", true, "The word size used to find closest targets during prefilter. (default " + ProteinWordGenerator.WORDSIZE 
                + " for protein, " + GoodWordIterator.DEFAULT_WORDSIZE  + " for nucleotide)");
        options.addOption("n", false, "Remove Ns from the query. Default is false");
        options.addOption("t", "threads", true, "#Threads to use. This process is CPU intensive. (default 1)");

        try {
            CommandLine line = new PosixParser().parse(options, args);

            if (line.hasOption("threads")) {
                maxThreads = Integer.valueOf(line.getOptionValue("threads"));
                if ( maxThreads >= Runtime.getRuntime().availableProcessors()) {
                   System.err.println(" Runtime.getRuntime().availableProcessors() " + Runtime.getRuntime().availableProcessors()); 
                }                
            } else {
                maxThreads = 1;
            }
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
            }else {
                out = new PrintStream(System.out);
            }
            if (line.hasOption('n')) {
                removeBaseN = true;
            }else {
                removeBaseN = false;
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
        
        SequenceType querySeqType = SeqUtils.guessSequenceType(queryFile);
        SequenceType refSeqType = SeqUtils.guessSequenceType(refFile);

        if ( querySeqType !=  refSeqType) {
            throw new RuntimeException("reference seqs and query seqs must be the same type, either protein or nucleotide. " );
        }
        final PairwiseKNN theObj = new PairwiseKNN( refFile, mode, k, wordSize, prefilter);
        
        final AtomicInteger outstandingTasks = new AtomicInteger();        
        ExecutorService service = Executors.newFixedThreadPool(maxThreads);
        out.println("#query file: " + queryFile.getName() + " db file: " + refFile.getName() + " k: " + k + " mode: " + mode + " usePrefilter: " + prefilter);
        out.println("#seqname\tk\torientation\tscore\tident\tquery_start\tquery_end\tquery_length\tref_start\tref_end\tref_seqid\tref_desc");
                
        SequenceReader queryReader = new SequenceReader(queryFile);
        Sequence seq;        
        while ( (seq = queryReader.readNextSequence()) !=null){            
            final Sequence threadSeq = seq;

            Runnable r = new Runnable() {
                public void run() {
                    try {
                        List<Neighbor> alignments = theObj.findMatch(threadSeq, removeBaseN);
                        theObj.printAlignment(threadSeq, alignments, out);
                        outstandingTasks.decrementAndGet();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            };

            outstandingTasks.incrementAndGet();
            service.submit(r);

            while (outstandingTasks.get() >= maxTasks);   
        }
        
        service.shutdown();
        service.awaitTermination(1, TimeUnit.DAYS);        
        queryReader.close();
        out.close();        
    }
}
