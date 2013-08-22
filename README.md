##RDP Alignment tools package

### Intro
The Alignment tools package contains commands to do defined community analysis (Fish et al., 2013. FunGene: the functional gene pipeline and repository. Frontiers, Accepted), 
pairwise alignment and hidden markov model (HMMER3 models, no training).

### Setup
This project depends on https://github.com/rdpstaff/ReadSeq and https://github.com/rdpstaff/fungene_pipeline. See RDPTools (https://github.com/rdpstaff/RDPTools) to install.

### Usage
	
* Defined Community Analysis

	This command compares input nucleotide reads to the set of known sequences for amplification targets in the sequenced DNA. 
	It computes a global alignment between the read and each of the reference sequences. The pairwise alignment producing the highest alignment score 
	is used to identify the source organism, substitutions, indels and associated quality scores. 
	These output files can be processed by the parseErrorAnalysis.py script (https://github.com/rdpstaff/fungene_pipeline) to summarize the types of errors found in the input reads. 
	The tab-delimited summary file can be imported to Excel programs to make graphs.

		java -jar /path/to/AlignmentTools.jar compare-error-type
		usage: CompareErrorType [options] <ref_nucl> (<query_nucl> | <query_nucl.fasta> <query_nucl.qual>)
				-s,--stem <arg>   Output stem (default <query_nucl.fasta>)
 
	The sample input and output files can be downloaded from RDP tutorial http://rdp.cme.msu.edu/tutorials/defined_community/RDPtutorial_DEFINED-COM.html.
	Example commands:
		
		java -jar /path/to/AlignmentTools.jar compare-error-type -s mid01 /path/to/nifH_control_refseq_nucl_slice.fa  /path/to/mid01_trimmed.fasta /path/to/mid01_trimmed.qual
				
		/path/to/fungene_pipeline/parseErrorAnalysis.py -i chimera_contaminant.id -q mid01_qual.txt mid01_alignments.txt mid01_mismatches.txt mid01_indels.txt /path/to/nifH_control_refseq_nucl_slice.fa> mid01_erroranalysis_summary.txt

* Merge Alignment	
	This command merges multiple alignment files into one file. The input files can be in fasta or stockholm format, each with a reference sequence named "#=GC_RF" to mark model positions.
		
		java -jar /path/to/AlignmentTools.jar alignment-merger alignment merged_aligned.fasta
		
* Compute k-nearest-neighbors by pairwise alignment
	
	This command computes pairwise alignment between the read and each reference sequences. It returns the alignments with the top k highest scores. 
	It reports the following for each alignment: 	
	seqname	k	ref_seqid	ref_desc	orientation	score	ident	query_start	query_end	query_length	ref_start	ref_end

		java -jar /path/to/AlignmentTools.jar pairwise-knn         
		usage: PairwiseKNN <options> <queryFile> <dbFile>
 		-k <arg>          K-nearest neighbors to return
 		-m,--mode <arg>   Alignment mode {global, glocal, local, overlap,
                   overlap_trimmed} (default= glocal)
		-o,--out <arg>    Redirect output to file instead of stdout

		        
