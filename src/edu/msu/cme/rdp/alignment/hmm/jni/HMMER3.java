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
package edu.msu.cme.rdp.alignment.hmm.jni;

import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import java.io.File;
import java.util.Iterator;

/**
 *
 * @author fishjord
 */
public class HMMER3 {
    private final int numModels;

    public HMMER3(String hmmdb) throws Exception {
        System.loadLibrary("hmmerwrapper");
        numModels = initHmmer(hmmdb);
    }

    public HMMER3(String hmmdb, String hmmerWrapperLib) throws Exception {
        System.load(hmmerWrapperLib);
        numModels = initHmmer(hmmdb);
    }

    @Override
    public void finalize() throws Throwable {
        super.finalize();
        destroyHmmer();
    }

    public HMMER3Hit[] findHits(Sequence seq) {
	return hmmer3(seq.getSeqString());
    }

    public HMMER3Hit[] findHits(String seq) {
	return hmmer3(seq);
    }

    private native int initHmmer(String hmmdb);
    private native void destroyHmmer();
    private native HMMER3Hit[] hmmer3(String seq);

    public static void main(String[] args) throws Exception {
	if(args.length != 2) {
	    System.err.println("USAGE: test <hmm> <seqs>");
	    System.exit(1);
	}
        HMMER3 hmmer = new HMMER3(args[0]);

	SeqReader reader = new SequenceReader(new File(args[1]));
	Sequence seq;

	while((seq = reader.readNextSequence()) != null) {
	    for(HMMER3Hit hit : hmmer.findHits(seq)) {
		System.out.println(hit.modelName + "\t" + hit.bits + "\t" + hit.hmmStart + "\t" + hit.hmmEnd + "\t" + hit.seqStart + "\t" + hit.seqEnd);
                System.out.println(hit.alignedSeq);
	    }
	}

	//	HMMER3 hmmer = new HMMER3(args[0]);
	//hmmer.hmmer3("GRGVITSINFLEENGAYDDVDYVSYDVLGDVVCGGFAMPIRENKAQEIYIVMSGEMMALYAANNIAKGILKYANSGGVRLGGLICNERKTDRELELAEAL");
    }
}
