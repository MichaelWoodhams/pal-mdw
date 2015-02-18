package pal.distance;

import pal.math.IntMatrix;
import Jama.Matrix;

/**
 * Static class.
 * Calculates logdet distance (and in some cases the standard error in the distance) for
 * given divergence matrix, which can be either integer or real.
 * 
 * The definition of logdet distance used here includes a 1/N factor (where N is the number of states).
 * This factor is in the LogDet definition in Lockhart et al. MBE 11 605 (1994)
 * and in Massingham and Goldman, MBE 24 2277 (2007) ("M&G") but not in
 * the paralinear distance of Lake, PNAS 91 1455 (1994). 
 * 
 * @author Michael Woodhams, School of IT, University of Sydney
 * @date 2009-05-05.
 *
 */

public class LogDetDistance {
	public static double MAX_DISTANCE = 1000;
	/**
	 * LogDet distance from a real valued divergence matrix - e.g. theoretical expected values.
	 * Does not require the matrix to sum to one - calculation is insensitive to arbitrary scaling of matrix.
	 * Returns a maximum value of 1000 (which is effectively indistinguishible from infinity.)
	 * @param divergence
	 * @return distance
	 */
	public static double distance(Matrix divergence) {
		
		int numStates = divergence.getRowDimension();

		double[] colSum = new double[numStates];  // "C(D)" in M&G
		double[] rowSum = new double[numStates];  // "R(D)" in M&G
		for (int row=0; row<numStates; row++) {
			for (int col=0; col<numStates; col++) {
				colSum[col] += divergence.get(row, col);
				rowSum[row] += divergence.get(row, col);
			}
		}
		double det = (double) divergence.det();
		double rowColumnSumLogs = 0;
		for (int state = 0; state < numStates; state++) {
			rowColumnSumLogs += Math.log(colSum[state]*rowSum[state]);
		}
		double dist = -(Math.log(det)-0.5*rowColumnSumLogs)/numStates;
		// Sampling noise can cause 'det' to be negative or very small.
		// Cap the output at 1000:
		if (Double.isNaN(dist) || Double.isInfinite(dist) || dist > MAX_DISTANCE ) { dist = MAX_DISTANCE; }
		return dist;
	}

	
	public static double distance(double[][] divergence) {
		return distance(new Matrix(divergence));
	}
	
	
	/**
	 *  LogDet distance from an integer divergence matrix (actual count of sites)
	 * @param divergence - divergence matrix
	 * @return distance
	 */
	public static double distance(IntMatrix divergence) {
		return distancePossiblyWithError(divergence, false)[0];
	}
	public static double distance(int[][] divergence) {
		return distance(new IntMatrix(divergence));
	}
	
	/**
	 *  LogDet distance from an integer divergence matrix (actual count of sites)
	 *  with standard error. 
	 *  Std Error of logdet distance taken from Massingham and Goldman, MBE 24 2277 (2007) ("M&G")
	 *  but **I HAVE NOT TESTED THIS!!**. Convince yourself it is correct before using for anything that matters.
	 * @param divergence
	 * @return array of {distance, std error in distance}.
	 */
	public static double[] distanceWithError(IntMatrix divergence) {
		return distancePossiblyWithError(divergence, true);
	}
	
	private static double[] distancePossiblyWithError(IntMatrix divergence, boolean withError) {
		double[] answer = new double[2];
		
		int numStates = divergence.columns;

		int[] colSum = new int[numStates];  // "C(D)" in M&G
		int[] rowSum = new int[numStates];  // "R(D)" in M&G
		int numSites = 0;
		for (int row=0; row<numStates; row++) {
			for (int col=0; col<numStates; col++) {
				colSum[col] += divergence.m[row][col];
				rowSum[row] += divergence.m[row][col];
				numSites += divergence.m[row][col];
			}
		}
		double det = (double) divergence.det();
		double rowColumnSumLogs = 0;
		for (int state = 0; state < numStates; state++) {
			rowColumnSumLogs += Math.log(colSum[state]*rowSum[state]);
		}
		
		answer[0] = -(Math.log(det)-0.5*rowColumnSumLogs)/numStates;
		
		if (!withError) {
			return answer;
		}

		// divAdj = det *(div)^-1 
		IntMatrix divAdj = divergence.adjugate();
		// Matrix 'u' as defined in M&G (convert to real numbers at this point - although we could stick with ints for a little longer.)
		Matrix u = new Matrix(numStates,numStates);
		for (int row = 0; row < numStates; row++) {
			for (int col = 0; col < numStates; col++) {
				u.set(row,col,  1./(2*numStates)*(1./colSum[col]+1./rowSum[row]-2*divAdj.m[col][row]/det));
			}
		}
		Matrix uSquared = u.times(u);
		double sum = 0;
		for (int row = 0; row < numStates; row++) {
			for (int col = 0; col < numStates; col++) {
				sum += uSquared.get(row, col)*divergence.m[row][col];
			}
		}
		answer[1] = Math.sqrt(sum/numSites); // standard error in distance
		return answer;
	}
	

}
