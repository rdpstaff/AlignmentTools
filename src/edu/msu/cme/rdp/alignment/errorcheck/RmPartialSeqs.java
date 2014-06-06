/*
 * Copyright (C) 2013 wangqion
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package edu.msu.cme.rdp.alignment.errorcheck;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAligner;
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAlignment;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.alignment.pairwise.rna.DistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.IdentityDistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.OverlapCheckFailedException;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.utils.kmermatch.KmerMatchCore;
import edu.msu.cme.rdp.readseq.utils.kmermatch.NuclSeqMatch;
import edu.msu.cme.rdp.readseq.utils.kmermatch.ProteinSeqMatch;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author wangqion
 */
public class RmPartialSeqs {
    private static final char gapChar = '-';     
    private ScoringMatrix scoringMatrix ;
    private SequenceType seqType ;
    private AlignmentMode mode = AlignmentMode.overlap;
    private static final Options options = new Options();
    private static DistanceModel dist = new IdentityDistanceModel();
    private HashMap<String, Sequence> refSeqMap = new HashMap<String, Sequence>();
    private ArrayList<Sequence> seqList = new ArrayList<Sequence>();
    private KmerMatchCore sabCalculator = null;
    private int knn = 20;
    private int min_begin_gaps = 50;
    private int min_end_gaps = 50;

    
    static {          
        options.addOption("a", "alignment-mode", true, "Alignment mode: overlap, glocal, local or global. default = overlap");
        options.addOption("g", "min_gaps", true, "The minimum number of continuous gaps in the beginning or end of the query alignment. If above the cutoff, the query is marked as partial. default = 50");
        options.addOption("k", "knn", true, "The top k closest targets using a heuristic method. (default = 20)");
        options.addOption("o", "alignment-out", true, "The output file containing the pairwise alignment");
    }
    
    
    public RmPartialSeqs(String trainseqFile, String testFile, AlignmentMode mode, int knn, int min_gaps) throws IOException, OverlapCheckFailedException{
        this.mode = mode;
        this.knn = knn;
        this.min_begin_gaps = min_gaps;
        this.min_end_gaps = min_gaps;
        
        SequenceReader parser = new SequenceReader(new File(trainseqFile));
        Sequence seq;
        SequenceType seqType = null;
        while ( (seq=parser.readNextSequence()) != null) {
            if ( seqType == null){
                seqType = SeqUtils.guessSequenceType(seq);
            }
            refSeqMap.put(seq.getSeqName(), seq);
        }
        parser.close();
   
        parser = new SequenceReader(new File(testFile));        
        while ( (seq=parser.readNextSequence()) != null) {
            seqList.add( seq);
        }
        parser.close();
        
        if (seqType == SequenceType.Nucleotide) {
            scoringMatrix = ScoringMatrix.getDefaultNuclMatrix();
            sabCalculator = new NuclSeqMatch(trainseqFile);
        } else {
            scoringMatrix = ScoringMatrix.getDefaultProteinMatrix();
            sabCalculator = new ProteinSeqMatch(trainseqFile);
        }
                
    }
       
    public HashSet<Sequence> checkPartial(PrintStream seqOutstream, PrintStream alignOutstream) throws OverlapCheckFailedException, IOException{
        HashSet<Sequence> partialSeqs = new HashSet<Sequence>();
        for ( int i= 0; i < seqList.size(); i++){
            Sequence seqx = seqList.get(i);
            PairwiseAlignment bestResult = null;
            int bestScore = Integer.MIN_VALUE;
            Sequence bestSeqy = null;
            
            ArrayList<NuclSeqMatch.BestMatch> matchResults = sabCalculator.findTopKMatch(seqx, knn);
            for (NuclSeqMatch.BestMatch match: matchResults){
            
                Sequence seqy = refSeqMap.get(match.getBestMatch().getSeqName());
                PairwiseAlignment result = PairwiseAligner.align(seqx.getSeqString().replaceAll("U", "T"), seqy.getSeqString().replaceAll("U", "T"), scoringMatrix, mode);
               
                if ( bestResult == null || result.getScore() >= bestScore){               
                    bestResult = result;
                    bestScore = result.getScore() ;
                    bestSeqy = seqy;                                       
                }
                
            }
            double distance = dist.getDistance(bestResult.getAlignedSeqj().getBytes(), bestResult.getAlignedSeqi().getBytes(), 0);
            
            int beginGaps = getBeginGapLength(bestResult.getAlignedSeqi());
            int endGaps = getEndGapLength(bestResult.getAlignedSeqi()) ;
            if ( ( beginGaps >= this.min_begin_gaps) || ( endGaps >= this.min_end_gaps)){
                partialSeqs.add(seqx);                
            }else {
                seqOutstream.println(">" + seqx.getSeqName() + "\t" + seqx.getDesc() + "\n" + seqx.getSeqString());
            }
            if ( alignOutstream != null){
                alignOutstream.println(">\t" + seqx.getSeqName() + "\t" + bestSeqy.getSeqName() + "\t" + String.format("%.3f", distance) 
                        + "\tmissingBegin=" + ( beginGaps >= this.min_begin_gaps) + "\tmissingEnd=" + ( endGaps >= this.min_end_gaps)
                        + "\tbeginGaps=" + beginGaps + "\tendGaps=" + endGaps);
                alignOutstream.print(bestResult.getAlignedSeqi() + "\n");
                alignOutstream.print(bestResult.getAlignedSeqj() + "\n");
            }
        }   
        seqOutstream.close();
        if ( alignOutstream != null) alignOutstream.close();
        
        return partialSeqs;
    }
    
   
    /**
     * 
     * @param s
     * @return the number of continuous gaps in the beginning
     */
    int getBeginGapLength(String s){
        int length = 0;
        for ( int i = 0; i < s.length(); i++){
            if ( s.charAt(i) == '-'){
                length ++;
            }else {
                return length;
            }
        }
        return length;
    }
    
        /**
     * 
     * @param s
     * @return the number of continuous gaps in the end
     */
    int getEndGapLength(String s){
        int length = 0;
        for ( int i = s.length() -1; i > 0; i--){
            if ( s.charAt(i) == '-'){
                length ++;
            }else {
                return length;
            }
        }
        return length;
    }
    
    
     /**
     * This program detects partial sequences based on the best pairwise alignment for each query sequence, 
     * @param args
     * @throws Exception 
     */
    public static void main(String[] args) throws Exception {
       
        String trainseqFile = null;
        String queryFile = null;   
        PrintStream seqOutStream = null;
        PrintStream alignOutStream = null;
        AlignmentMode mode = AlignmentMode.overlap;
        int k = 10;
        int min_gaps = 50;
        
        
         try {
            CommandLine line = new PosixParser().parse(options, args);

            if (line.hasOption("alignment-mode")) {
                String m = line.getOptionValue("alignment-mode").toLowerCase();
                mode = AlignmentMode.valueOf(m);              
            }
           
            if (line.hasOption("min_gaps")) {
                min_gaps = Integer.parseInt(line.getOptionValue("min_gaps"));
            }
            if (line.hasOption("knn")) {
                k = Integer.parseInt(line.getOptionValue("knn"));
            }
            if (line.hasOption("alignment-out")) {
                alignOutStream = new PrintStream(new File(line.getOptionValue("alignment-out")));
            }
            args = line.getArgs(); 
            if ( args.length != 3){
                throw new Exception("wrong number of arguments");
            }
                       
            trainseqFile = args[0];
            queryFile = args[1];
            seqOutStream = new PrintStream(new File(args[2]));
         }catch (Exception e) {
             System.err.println("Error: " + e.getMessage());
             new HelpFormatter().printHelp(80, " [options] fulllengthSeqFile queryFile passedSeqOutFile\n  sequences can be either protein or nucleotide", "", options, "");
            return;
        }
              
        RmPartialSeqs theObj = new RmPartialSeqs(trainseqFile, queryFile, mode, k, min_gaps);
        
        theObj.checkPartial(seqOutStream, alignOutStream);
   }
}
