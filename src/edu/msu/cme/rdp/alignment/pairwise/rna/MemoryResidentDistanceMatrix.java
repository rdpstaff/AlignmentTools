package edu.msu.cme.rdp.alignment.pairwise.rna;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class MemoryResidentDistanceMatrix {
	List<DistMatrixSeq> xSeqs;
	List<DistMatrixSeq> ySeqs;
	DistanceModel distanceModel;
	double[][] distanceMatrix;

	protected MemoryResidentDistanceMatrix(List<DistMatrixSeq> sequenceList, DistanceModel distanceModel, int overlapLimit ) throws OverlapCheckFailedException {
        this.xSeqs = this.ySeqs = sequenceList;
		this.distanceModel = distanceModel;

		distanceMatrix = new double[sequenceList.size()][sequenceList.size()];
		Map<String, Set<String>> errorMap = new HashMap<String, Set<String>>();

        for(int seq1 = 0;seq1 < sequenceList.size();seq1++) {
            distanceMatrix[seq1][seq1] = 0;
            Set<String> nonOverlappingSet = new HashSet();
            for(int seq2 = seq1 + 1;seq2 < sequenceList.size();seq2++) {
				try {
                    double dist = distanceModel.getDistance( sequenceList.get(seq1).getSeqBytes(), sequenceList.get(seq2).getSeqBytes(), overlapLimit );
					distanceMatrix[seq1][seq2] = dist;
                    distanceMatrix[seq2][seq1] = dist;
				} catch (OverlapCheckFailedException ex) {
					// set the value to Double.NaN
					distanceMatrix[seq1][seq2] = Double.NaN;
                    distanceMatrix[seq2][seq1] = Double.NaN;
					// add to the error set
					nonOverlappingSet.add(sequenceList.get(seq2).getSeqid());

				}
            }
			if ( ! nonOverlappingSet.isEmpty() ) {
				errorMap.put(sequenceList.get(seq1).getSeqid(), nonOverlappingSet);
			}
        }

		if ( ! errorMap.isEmpty()) {
			throw new OverlapCheckFailedException(errorMap);
		}
	}
	
	protected MemoryResidentDistanceMatrix(List<DistMatrixSeq> xSeqList, List<DistMatrixSeq> ySeqList, DistanceModel distanceModel, int overlapLimit ) throws OverlapCheckFailedException {
		this.xSeqs = xSeqList;
		this.ySeqs = ySeqList;
		this.distanceModel = distanceModel;
		
		distanceMatrix = new double[xSeqList.size()][ySeqList.size()];
		int x = 0;
		Map<String, Set<String>> errorMap = new HashMap<String, Set<String>>();
		for(DistMatrixSeq xSeq : xSeqList) {
			int y = 0;
			Set <String> nonOverlappingSet = new HashSet<String>();
			for (DistMatrixSeq ySeq : ySeqList) {
				try {
					distanceMatrix[x][y] = distanceModel.getDistance( xSeq.getSeqBytes(), ySeq.getSeqBytes(), overlapLimit );
				} catch (OverlapCheckFailedException ex) {					
					// set the value to Double.NaN
					distanceMatrix[x][y] = Double.NaN;
					// add to the error set					
					nonOverlappingSet.add(ySeq.getSeqid());
					
				}
                y++;
			}
			if ( ! nonOverlappingSet.isEmpty() ) {
				errorMap.put(xSeq.getSeqid(), nonOverlappingSet);
			}
            x++;
		}
		if ( ! errorMap.isEmpty()) {
			throw new OverlapCheckFailedException(errorMap);
		}
	}
	
	public List<DistMatrixSeq> getXSeqs() {
		return xSeqs;
	}
	
	public List<DistMatrixSeq> getYSeqs() {
		return ySeqs;
	}
	
	public double[][] getDistanceMatrix() {
		return distanceMatrix;
	}
	
	public double getDistance( int x, int y ) {
		return distanceMatrix[x][y];
	}
}
