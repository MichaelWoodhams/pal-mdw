// NucleotideModel.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.substmodel;

import pal.misc.*;
import pal.datatype.*;

import java.io.*;


/**
 * base class for nucleotide rate matrices
 *
 * @version $Id: NucleotideModel.java,v 1.10 2003/11/30 05:29:22 matt Exp $
 *
 * @author Korbinian Strimmer
 */
abstract public class NucleotideModel extends AbstractRateMatrix implements RateMatrix, Serializable {
	//
	// Public stuff
	//

	/**
	 * Create nucleotide substitution model according to model type
	 *
	 * @param modelID model code
	 * @param params model parameters
	 * @param freq  model frequencies
	 *
	 * @return nucleotide rate matrix
	 */
	public static NucleotideModel getInstance(int modelID, double[] params, double[] freq)
	{
		if (modelID == NucleotideModelID.GTR)
		{
			return new GTR(params, freq);
		}
		else if (modelID == NucleotideModelID.TN)
		{
			return new TN(params, freq);
		}
		else if (modelID == NucleotideModelID.HKY)
		{
			return new HKY(params, freq);
		}
		else if (modelID == NucleotideModelID.F84)
		{
			return new F84(params, freq);
		}
		else if (modelID == NucleotideModelID.F81)
		{
			return new F81(freq);
		}
		else
		{
			return new F81(freq);
		}
	}

	// interface Report (inherited, remains abstract)

	// interface Parameterized (inherited, remains abstract)


	//
	// Protected stuff (for use in derived classes)
	//

	// Constructor
	protected NucleotideModel(double[] f)
	{
		// Dimension = 4
		super(4);

		setDataType(new Nucleotides());
		setFrequencies(f);
	}

	protected void printFrequencies(PrintWriter out)
	{
		out.println("Nucleotide frequencies:");
		super.printFrequencies(out);
	}

	protected void printRatios(PrintWriter out)
	{
		computeRatios();
		out.print("Expected transition/transversion ratio: ");
		format.displayDecimal(out, expectedTsTvRatio, 2);
		out.println();
		out.print("Expected pyrimidine transition/purine transition ratio: ");
		format.displayDecimal(out, expectedYRTsRatio, 2);
		out.println();
	}
	/*protected void reportFrequencies(ReportContainer rc) {
		rc.addItem("Nucleotide frequencies",null, getFrequencies());
	}
	protected void reportRatios(ReportContainer rc) {
		computeRatios();
		ReportContainer ratio = rc.addContainer("Ratios",null);
		ratio.addItem("Expected transition/transversion ratio",new Integer(expectedTsTvRatio));
		ratio.addItem(("Expected pyrimidine transition/purine transition ratio",new Integer(expectedYRTsRatio));
	}*/

	//
	// Private stuff
	//

	private void computeRatios() {
		double[] frequency = getEquilibriumFrequencies();
		double[][] rate = getRelativeRates();
		// Compute expectation ratios
		int A = 0; double piA = frequency[0];
		int C = 1; double piC = frequency[1];
		int G = 2; double piG = frequency[2];
		int T = 3; double piT = frequency[3];

		double numYTs = piC*rate[C][T] + piT*rate[T][C];
		double numRTs = piA*rate[A][G] + piG*rate[G][A];
		double numTv =
			piA*rate[A][C] + piC*rate[C][A] +
			piA*rate[A][T] + piT*rate[T][A] +
			piC*rate[C][G] + piG*rate[G][C] +
			piG*rate[G][T] + piT*rate[T][G];

		expectedTsTvRatio = (numYTs + numRTs)/numTv;
		expectedYRTsRatio = numYTs/numRTs;
	}

	// Expected transition-transversion ratio
	private double expectedTsTvRatio;

	// Expected ratio of pyrimidine and purine transitions
	private double expectedYRTsRatio;
}
