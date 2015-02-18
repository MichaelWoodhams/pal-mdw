package pal.alignment;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import pal.misc.*;

import pal.datatype.Nucleotides;

/**
 * Align a set of sequences. 
 * 
 * This is a very preliminary class - it should be much more capable than it is.
 * 
 * For now, the only method available is calling the external program 'MUSCLE'.
 *    http://www.drive5.com/muscle/,
 *    Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Research 32(5), 1792-97. 
 * I've tried to allow for other methods to be easily added. I've used MUSCLE because it was one of a handful
 * recommended to me (Mike Charleston, personal communication) and because it can read/write sequences
 * via stdin and stdout, so I don't have to mess with temporary files.
 * 
 * ALSO: For now, it is only possible to use the default settings on MUSCLE. We should be able to set parameters
 * such as gap penalties, and whether to use slow-and-accurate or fast-and-loose alignment etc.
 * 
 * @author Michael Woodhams
 *
 */

public class Align {
	public static final int methodMuscle = 1;
	
	/**
	 * Take the sequences in an alignment and run a multiple alignment method on them.
	 * (Uses the default method - MUSCLE) 
	 * @param rawAlignment (does not get modified)
	 * @return Aligned alignment
	 */
	public static Alignment align(Alignment rawAlignment) {
		return align(rawAlignment, methodMuscle);
	}
	
	/**
	 * Take the sequences in an alignment and run a multiple alignment method on them
	 * using the specified method (e.g. Align.methodMuscle)
	 * 
	 * TO DO: allow many extra parameters - speed vs accuracy, gap penalty, ...
	 * 
	 * @param rawAlignment (does not get modified)
	 * @param method An integer indicating which alignment method to use.
	 * @return Aligned alignment
	 */
	public static Alignment align(Alignment rawAlignment, int method) {
		switch (method) {
		case methodMuscle :
			return alignByMuscle(rawAlignment);
		default:
			throw new RuntimeException("Unrecognized alignment method");
		}
	}
	
	
	/**
	 * Take the sequences in an alignment and run a multiple alignment method on them
	 * using the specified method (e.g. Align.methodMuscle)
	 * 
	 * @param rawAlignment rawAlignment (does not get modified)
	 * @param subset An array of sequence indicies from rawAlignment indicating which sequences to align
	 * @param method An integer indicating which alignment method to use.
	 * @return The multiple alignment of the chosen sequences
	 */
	public static Alignment alignSubset(Alignment rawAlignment, int[] subset, int method) {
		
		int n = subset.length;
		Identifier[] ids = new Identifier[n];
		String[] sequences = new String[n];
		for (int i=0; i<n; i++) {
			int seqNumber = subset[i];
			ids[i] = rawAlignment.getIdentifier(seqNumber);
			sequences[i] = rawAlignment.getAlignedSequenceString(seqNumber);
		}
		SimpleAlignment subalignment = new SimpleAlignment(ids, sequences, null, rawAlignment.getDataType());
		return align(subalignment, method);
	}
	
	/**
	 * As above, but use default alignment method (MUSCLE).
	 * @param rawAlignment
	 * @param subset
	 * @return
	 */
	public static Alignment alignSubset(Alignment rawAlignment, int[] subset) {
		return alignSubset(rawAlignment, subset, methodMuscle);
	}
	
	/*
	 * For now, just die if anything fails (throw RuntimeException.)
	 * TODO: allow a more graceful recovery.
	 * 
	 */
	
	private static Alignment alignByMuscle(Alignment rawAlignment) {
		
		Alignment aligned = null;
		try {
			// Create the 'muscle' process (running the external program)
			Process process = Runtime.getRuntime().exec("muscle");
			// Attach toProc and fromProc to its stdin and stdout respectively.
			// printFASTA requires a PrintWriter
			PrintWriter toProc = new PrintWriter(process.getOutputStream());
			// readFastaSequences requires a Reader
			BufferedReader fromProc = new BufferedReader(
								new InputStreamReader(process.getInputStream()));
			// Send the raw alignment in FASTA format to muscle
			AlignmentUtils.printFASTA(rawAlignment, toProc);
			toProc.close();
			// Pipe muscle's output (also FASTA format) through an alignment reader to get an alignment
			aligned = AlignmentReaders.readFastaSequences(fromProc, Nucleotides.DEFAULT_INSTANCE);
			fromProc.close();
			process.destroy();	
		}
		catch (IOException e) {
			System.out.println("Test failed with an IOException: " + e.getMessage());
			e.printStackTrace();
		}
		return aligned;
	}
}
