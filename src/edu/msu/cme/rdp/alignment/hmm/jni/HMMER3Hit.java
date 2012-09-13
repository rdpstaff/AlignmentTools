package edu.msu.cme.rdp.alignment.hmm.jni;

public class HMMER3Hit {
    String modelName;
    String alignedSeq;
    double bits;
    int hmmStart, hmmEnd;
    int seqStart, seqEnd;
    int hmmLength;

    public double getBits() {
        return bits;
    }

    public String getModelName() {
        return modelName;
    }

    public String getAlignedSeq() {
        return alignedSeq;
    }

    public int getHmmEnd() {
        return hmmEnd;
    }

    public int getHmmStart() {
        return hmmStart;
    }

    public int getSeqEnd() {
        return seqEnd;
    }

    public int getSeqStart() {
        return seqStart;
    }
}
