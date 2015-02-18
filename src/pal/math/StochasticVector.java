package pal.math;

import pal.misc.Parameterized;

/**
 * A stochastic vector is one where all entries are non-negative and which sums to one.
 * Often we need to find an an optimal stochastic vector, but our optimisation routines
 * aren't set up to directly handle the 'sums to one' constraint. This class works
 * around that restriction. 
 * 
 * This class provides a function to map (n-1) arguments in the range 0-1 to a length n
 * stochastic vector, with the following properties:
 * * The mapping is continuous and (almost) monotonic 
 * * All possible arguments map to a stochastic vector
 * * All possible stochastic vectors are within the domain of the mapping.
 * 
 * Each argument specifies what proportion of the remaining probability its corresponding
 * entry in the stochastic vector gets. I.e, if 'a' = arguments, 'v' = vector then
 * 
 * v[i] = a[i]*(1- sum(v[0..i-1]))
 * v[n] = (1- sum(v[0..n-1])
 * 
 * Also provided: A similar encoding for normalized rate matrices. 
 * Constraints are: off diagonals are non-negative, rows sum to zero, diagonal sums to -1.
 * For n x n matrix, that gives n(n-1)-1 free parameters. (For n=4, 11 parameters.)
 * To allow easy upper and lower bounds for the parameters, I construct them as follows.
 * All rate matrix parameters are in the range 0-1.
 * Parameters are (for n=4) (D1, D2, A, B, C, a1, a2, b1, b2, c1, c2, d1, d2))
 * D1, D2 = branch lengths (0 to infinity) 
 * R[0][0] = -A
 * The three off-diagonals on row 0 are -R[0][0]*a1, -R[0][0]*(1-a1)*a2, -R[0][0]*(1-a1)*(1-a2)
 * R[1][1] = -B*(1-A)  
 * The three off-diagonals on row 1 are are -R[1][1]*b1, -R[1][1]*(1-b1)*b2, -R[1][1]*(1-b1)*(1-b2)
 * R[2][2] = -C*(1-B)*(1-A)
 * (off-diagonals from c1, c2, as above)
 * R[3][3] = -(1-C)*(1-B)*(1-A)
 * (off-diagonals from d1, d2, as above)
 * Thus any set of parameters in the range 0-1 will give a legitimate normalized rate matrix, and
 * all legitimate rate matrices have a parameterization.
 * 
 * Inverses of both mappings are provided, and various interfaces.
 * 
 * This is a static class, containing methods only.
 * 
 * Possible TODO: allow it to be a real object, implementing Parameterized. (But still keep the static methods.)
 * 
 * TODO: Rate matrix methods which read/write Jama Matrices instead of double[][]'s. 
 * (And then see if any of the invocations of the old routines should use the new ones instead.)
 * 
 * TODO: Tests for these routines are in RateMatrixLikelihood. Move them here.
 * 
 * @author woodhams
 *
 */

public class StochasticVector implements Parameterized {
	double[] vector_;
	double[] parameters_;
	int n_; // = vector_.length;
	
	/**
	 * This constructor starts with a 'flat' vector: all values equal. 
	 * @param len
	 */
	public StochasticVector(int len) {
		allocateSpace(len);
		this.setFlat();
	}
	/**
	 * 
	 * @param vector: a vector with non-negative elements summing to one
	 */
	public StochasticVector(double[] vector) {
		allocateSpace(vector.length);
		setVector(vector);
		makeParameters(vector_, parameters_, 0);
	}
	
	public StochasticVector clone() {
		return new StochasticVector(this.vector_);
	}
	
	/**
	 * @return The stochastic vector
	 */
	public double[] getVector() {
		return vector_;
	}
	/**
	 * @return The stochastic vector
	 */
	public void getVector(double[] store) {
		makeVector(parameters_,store);
	}
	/**
	 * Set the stochastic vector (validity of vector not checked.)
	 * @param vector
	 */
	public void setVector(double[] vector) {
		allocateSpace(vector.length);
		vector_ = vector.clone();
		makeParameters(vector_, parameters_, 0);
	}
	/**
	 * Set vector to be flat: i.e. vector = [1/n, 1/n ... 1/n].
	 */
	public void setFlat() {
		for (int i=0; i<n_; i++) {
			vector_[i] = 1./n_;
		}
		makeParameters(vector_, parameters_, 0);
	}
	/**
	 * @return The length of the stochastic vector (one more than the number of parameters)
	 */
	// c.f. getNumParameters()
	public int length() {
		return n_;
	}
	/**
	 * If vector length has changed, updates n_ and allocates a new parameters_ array.
	 * Does NOT change vector_ array, as whatever is calling allocateSpace takes care of this.
	 * @param newLen
	 */
	private void allocateSpace(int newLen) {
		if (newLen != n_) {
			n_ = newLen;
			parameters_ = new double[n_-1];
		}
	}
	
	/**
	 * Methods for Parameterized interface:
	 */
	public double getDefaultValue(int n) { return 1./(n_-n); } // (used to be zero). This default makes a flat vector.
	public double getLowerLimit(int n)   { return 0; }
	public double getUpperLimit(int n)   { return 1; }
	public int getNumParameters()        { return n_-1; }
	public double getParameter(int n)    { return parameters_[n]; }
	public void setParameter(double param, int n) {
		parameters_[n] = param;
		makeVector(parameters_,vector_);
	}
	public void setParameterSE(double paramSE, int n) { throw new RuntimeException("Not implemented");	}
	
	
	/* **************************************************************************************************************
	 * The static methods:
	 */
	
	/**
	 * The most general method. 
	 * @param n The length of the stochastic vector to be generated
	 * @param parameters A vector containing within it the n-1 parameters to be mapped to the stochastic vector 
	 * @param parametersOffset The starting point within the parameters vector 
	 * @param vector A vector to store the results within
	 * @param vectorOffset Where within vector to start storing results.
	 */
	public static void makeVector(int n, double[] parameters, int parametersOffset, double[] vector, int vectorOffset) {
		double remainder = 1;
		for (int i=0; i<n-1; i++) {
			double value = remainder * parameters[i+parametersOffset];
			vector[i+vectorOffset] = value;
			remainder -= value;
		}
		vector[n-1+vectorOffset] = remainder;
	}
	
	/**
	 * Write stochastic vector into the start of provided storage space,
	 * using all of the provided parameters vector. 
	 * @param parameters
	 * @param vector
	 */
	public static void makeVector(double[] parameters, double[] vector) {
		makeVector(parameters.length+1, parameters, 0, vector, 0);
	}
	
	public static double[] makeVector(double[] parameters) {
		double[] vector = new double[parameters.length+1];
		makeVector(parameters, vector);
		return vector;
	}
	
	/**
	 * Inverse mapping: given a stochastic vector, find the corresponding parameters.
	 * @param n The length of the stochastic vector supplied
	 * @param vector A vector containing within it the stochastic vector 
	 * @param vectorOffset Where to find the stochastic vector within 'vector'
	 * @param parameters A vector to store the results within
	 * @param parametersOffset Where within vector to start storing results. (Will store (n-1) numbers.)
	 */
	public static void makeParameters(int n, double[] vector, int vectorOffset, double[] parameters, int parametersOffset) {
		double cumulativeSum = 0;
		// Extra complexity ('if' statements) ensures all values strictly in 0-1 range even
		// in the presence of rounding errors or a sequence of zeros at the end of 'vector'.
		for (int i=0; i<n-1; i++) {		
			if (cumulativeSum <1) {
				parameters[i+parametersOffset] = vector[i+vectorOffset]/(1-cumulativeSum);
				cumulativeSum += vector[i+vectorOffset];
				if (cumulativeSum > 1) {
					parameters[i+parametersOffset] = 1;
				}
			} else {
				// any value would do here 
				parameters[i+parametersOffset] = 0;
			}
		}
	}
	/**
	 * Inverse mapping: given a stochastic vector, find the corresponding parameters.
	 * @param vector The stochastic vector (of length n)
	 * @param parameters A vector to store the results within
	 * @param parametersOffset Where within vector to start storing results. (Will store (n-1) numbers.)
	 */
	public static void makeParameters(double[] vector, double[] parameters, int parametersOffset) {
			makeParameters(vector.length, vector, 0, parameters, parametersOffset);
	}
	
	/**
	 * Set the parameters so that they produce a flat stochastic vector.
	 * I.e. the vector produced by these parameters = [1/n, 1/n, ... 1/n].
	 * This is achieved by parameters = [1/n, 1/(n-1), 1/(n-2), ... 1/2].
	 * @param vectorLength
	 * @param parameters
	 * @param parametersOffset
	 */
	public static void getFlatParameters(int vectorLength, double[] parameters, int parametersOffset) {
		for (int i=0; i<vectorLength-1; i++) parameters[parametersOffset+i] = 1./(vectorLength-i);
	}
	public static void makeParametersFlat(double[] parameters) {
		getFlatParameters(parameters.length+1,parameters,0);
	}
	

	
	/**
	 * Map set of parameters in range 0-1 into a normalized rate matrix. 
	 * For a size n rate matrix, we require n(n-1)-1 parameters.
	 * 'n' is determined from the provided 'rateMatrix' array.
	 * @param parameters Vector containing the parameters (all in the range 0 to 1) somewhere in it
	 * @param parametersOffset Starting point in 'parameters'
	 * @param rateMatrix Pre-allocated storage space for the rate matrix 
	 */
	public static void parametersToNormRateMatrix(double[] parameters, int parametersOffset, double[][] rateMatrix) {
		int n = rateMatrix.length;
		double[] workspace = new double[n];
		
		// First n-1 parameters determine the row sums
		makeVector(n, parameters, parametersOffset, workspace, 0);
		for (int i=0; i<n; i++) {
			rateMatrix[i][i] = -workspace[i];
		}
		
		// Now set off-diagonals with the rest of the parameters, n-2 parameters per row
		// generating an n-1 length stochastic vector, which gets scaled to match the
		// diagonal element (already set.) 
		for (int row=0; row<n; row++) {
			makeVector(n-1, parameters, parametersOffset + (n-1) + row*(n-2), workspace, 0);
			double sum = -rateMatrix[row][row];
			for (int j=0; j<n-1; j++) {
				int col = j + (j>=row ? 1 : 0);
				rateMatrix[row][col] = sum*workspace[j];
			}
		}
	}
	
	
	public static void normRateMatrixToParameters(double[][] rateMatrix, double[] parameters, int parametersOffset) {
		int n = rateMatrix.length;
		// 'rowSum' is the sum of the off-diagonal elements in a row. As each row sums to zero, this
		// is also the negative of the diagonal element in that row.
		double[] rowSum = new double[n];
		for (int row=0; row<n; row++) {
			rowSum[row] = -rateMatrix[row][row];
		}
		// rowSum is a stochastic vector, store the parameterization (n-1 numbers) 
		makeParameters(n, rowSum, 0, parameters, parametersOffset);
		
		double[] offDiagonal = new double[n-1];
		for (int row=0; row<n; row++) {
			if (rowSum[row] > 0) {
				for (int j=0; j<n-2; j++) {
					int col = j + (j>=row ? 1 : 0);
					offDiagonal[j] = rateMatrix[row][col]/rowSum[row];
				};
				makeParameters(n-1, offDiagonal, 0, parameters, parametersOffset+(n-1)+row*(n-2));
			} else {
				// rare case - null row. Any values would do here
				for (int j=0; j<n-2; j++) {
					parameters[parametersOffset+(n-1)+row*(n-2)+j] = 0;
				}
			}
		}
	}
	
	/**
	 * Return parameterization of the given rate matrix.
	 * This method creates a new double[] to hold the parameterization.
	 * @param rateMatrix
	 * @return
	 */
	public static double[] normRateMatrixToParameters(double[][] rateMatrix) {
		int n = rateMatrix.length;
		double[] parameters = new double[n*(n-1)-1];
		normRateMatrixToParameters(rateMatrix, parameters, 0);
		return parameters;
	}
}
