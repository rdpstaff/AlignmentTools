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
package edu.msu.cme.rdp.alignment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.StringReader;
import java.io.StringWriter;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class AlignmentMergerTest {

    @Test
    public void testAlignmentMerger() throws IOException {

        Sequence[] seqs = new Sequence[]{
            new Sequence("seq1", "", "AAAt--A--uatGcu"),
            new Sequence("seq2", "", "gtCAGg--cC--Cat"),
            new Sequence("seq3", "", "AUTg--G--taGcat"),
            new Sequence("seq4", "", "ttacAGTgt--T--tUt"),
            new Sequence("seq5", "", "GGAgt--G--tG")};

        String[] expectedLines = {
            ">seq1",
            "....AAAt.--.A--tatGct.",
            ">seq2",
            "gt..CAGg.--cC--...Cat.",
            ">seq3",
            "....ATTg.--.G--ta.Gcat",
            ">seq4",
            "ttacAGTgt--.T--t..Tt..",
            ">seq5",
            "....GGAgt--.G--t..G...",
            ">#=GC_RF",
            "....xxx..xx.xxx...x..."
        };

        List<File> mergeFiles = new ArrayList();
        for (Sequence seq : seqs) {

            File tmp = File.createTempFile("temp", "file");
            tmp.deleteOnExit();
            mergeFiles.add(tmp);

            FastaWriter out = new FastaWriter(tmp);
            out.writeSeq(seq);

            String refSeq = seq.getSeqString().replaceAll("[a-z]", ".").replaceAll("[A-Z\\-]", "x");
            out.writeSeq("#=GC_RF", refSeq);

            out.close();
        }

        File mergeOut = File.createTempFile("temp", "file");
        mergeOut.deleteOnExit();

        Map<File, String> dropped = AlignmentMerger.mergeAlignment(mergeFiles, mergeOut, new HashSet());

        if (!dropped.isEmpty()) {
            fail("Expected no dropped files, but some were: " + dropped);
        }

        BufferedReader reader = new BufferedReader(new FileReader(mergeOut));
        String line;
        int lineno = 0;

        while ((line = reader.readLine()) != null) {
            assertEquals(expectedLines[lineno++], line.trim());
        }
        reader.close();
    }

    @Test
    public void oldAlignmentMergerTest() throws IOException {

        String stk = "# STOCKHOLM 1.0\n"
                + "#=GF AU    Infernal 0.81\n"
                + "\n"
                + "DE02039D02             .-------------------.-AAA--gg...-----------------\n"
                + "DE02039F03             aGAACCACATGGTTTTGATAtCAAAGAttttaTCGCTTGAAGATGGACT\n"
                + "#=GC SS_cons           .<<<<<____>>>>>>>>>>.,,,,<<.....>>,)))))))--)))))\n"
                + "#=GC RF                .atccCgcAtGggataagga.GAAAGa.....tCGCTaaagGATGggCC\n"
                + "\n"
                + "DE02039D02             GGT\n"
                + "DE02039F03             CGC\n"
                + "#=GC SS_cons           )))\n"
                + "#=GC RF                cGC\n"
                + "//\n";
        /*
        String stk2 = "# STOCKHOLM 1.0\n" +
        "#=GF AU    Infernal 0.81\n" +
        "\n" +
        "DE02039H01             GUGUCGCAUGGCACUUAUGu.CAAAGAuuuaUCGCUGAAAGAUGGCCUCGC.\n" +
        "DE02039F02             -GG----------------aa------....--------------------g\n" +
        "#=GC SS_cons           <<<<<____>>>>>>>>>>..,,,,<<....>>,)))))))--)))))))).\n" +
        "#=GC RF                auccCgcAuGggauaagga..GAAAGa....uCGCUaaagGAUGggCCcGC.\n" +
        "//\n";
         */

        String fasta = ">DE02039H01\nGTGTCGCATGGCACTTATGt.CAAAGAtttaTCGCTGAAAGATGGCCTCGC.\n"
                   + ">DE02039F02  \n-GG----------------aa------....--------------------g\n"
                   + ">#=GC_SS_cons\n<<<<<____>>>>>>>>>>..,,,,<<....>>,)))))))--)))))))).\n"
                   + ">#=GC_RF     \nauccCgcAuGggauaagga..GAAAGa....uCGCUaaagGAUGggCCcGC.\n";


        List<File> mergeFiles = new ArrayList();
        File tmpFile;
        BufferedWriter out;

        tmpFile = File.createTempFile("temp", "stk");
        tmpFile.deleteOnExit();
        mergeFiles.add(tmpFile);

        out = new BufferedWriter(new FileWriter(tmpFile));
        out.write(stk);
        out.close();

        tmpFile = File.createTempFile("temp", "fasta");
        tmpFile.deleteOnExit();
        mergeFiles.add(tmpFile);

        out = new BufferedWriter(new FileWriter(tmpFile));
        out.write(fasta);
        out.close();

        File mergeOut = File.createTempFile("temp", "merge");
        mergeOut.deleteOnExit();

        Map<File, String> dropped = AlignmentMerger.mergeAlignment(mergeFiles, mergeOut, new HashSet());
        
        if (!dropped.isEmpty()) {
            fail("Expected no dropped files, but some were: " + dropped);
        }

        String[] expectedLines =
                (">DE02039D02\n"
                + ".-------------------..-AAA--gg...-----------------GGT.\n"
                + ">DE02039F03\n"
                + "aGAACCACATGGTTTTGATAt.CAAAGAttttaTCGCTTGAAGATGGACTCGC.\n"
                + ">DE02039H01\n"
                + ".GTGTCGCATGGCACTTATGt.CAAAGAttta.TCGCTGAAAGATGGCCTCGC.\n"
                + ">DE02039F02\n"
                + ".-GG----------------aa------.....--------------------g\n"
                + ">#=GC_RF\n"
                + ".xxxxxxxxxxxxxxxxxxx..xxxxxx.....xxxxxxxxxxxxxxxxxxxx."
                ).split("\n");

        BufferedReader reader = new BufferedReader(new FileReader(mergeOut));
        String line;
        int lineno = 0;

        while ((line = reader.readLine()) != null) {
            assertEquals(expectedLines[lineno++], line.trim());
        }
        reader.close();
    }
}
