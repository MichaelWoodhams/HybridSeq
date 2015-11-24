package hybridseq;
import java.io.PrintWriter;

import pal.alignment.Alignment;

/**
 * You'd have thought this capability would be in PAL somewhere
 * @author woodhams
 *
 */
public class FastaWriter {
	public static final int DEFAULT_MAX_LINE_LENGTH = 80;
	public static void writeFasta(Alignment align, PrintWriter out) {
		writeFasta(align,out,DEFAULT_MAX_LINE_LENGTH);
	}
	public static void writeFasta(Alignment align, PrintWriter out, int maxLine) {
		for (int i=0; i<align.getIdCount(); i++) {
			out.printf(">%s\n", align.getIdentifier(i));
			String seq = align.getAlignedSequenceString(i);
			int j;
			for (j=0; j<seq.length()/maxLine; j++) {
				out.println(seq.substring(j*maxLine, (j+1)*maxLine));
			}
			if (seq.length()%maxLine != 0) out.println(seq.substring(j*maxLine));
		}
	}
}
