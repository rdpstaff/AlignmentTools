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

package edu.msu.cme.rdp.alignment.pairwise.rna;

import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.Serializable;

/**
 *
 * @author fishjord
 */
public class DistMatrixSeq implements Serializable {    
    private String seqid;
    private String seq;
    private String name = "";
    private String status = "";
    private byte[] seqBytes;

    public DistMatrixSeq(String seqid, String seq) {
        this.seqid = seqid;
        this.seq = seq;
    }

    public DistMatrixSeq(String seqid, String seq, byte[] seqBytes) {
        this.seqid = seqid;
        this.seq = seq;
        this.seqBytes = seqBytes;
    }
    
    public DistMatrixSeq(String seqid, String seq, String status) {
        this.seqid = seqid;
        this.seq = seq;
        this.status = status;
    }

    public void setSeqName(String name) {
        if(name != null)
            this.name = name;
    }

    public String getSeqName() {
        return name;
    }

    
    public String getStatus(){
    	return status;
    }
    
    protected void translateSeq(String seq) {
        seqBytes = SeqUtils.toBytes(seq);
    }

    public String getSeqid() {
        return seqid;
    }

    public String getSeq() {
        return seq;
    }

    public byte[] getSeqBytes() {
        if(seqBytes == null)
            translateSeq(this.seq);
        return seqBytes;
    }
}
