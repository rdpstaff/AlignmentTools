package edu.msu.cme.rdp.alignment.pairwise.rna;

import java.util.Map;
import java.util.Set;

public class OverlapCheckFailedException extends Exception {
	
	private Map <String, Set<String>> errorMap; // key is xseqid, value is set of seq ids that does not have enough overlapping 
	
	public OverlapCheckFailedException(String arg0) {
		super(arg0);
		// TODO Auto-generated constructor stub
	}
	
	public OverlapCheckFailedException(Map <String, Set<String>>arg0) {
		super();
		this.errorMap = arg0;
	}
	
	public Map<String, Set<String>> getErrorMap() {
		return this.errorMap;
	}
}
