// RateMatrixHandler.java
//
// (c) 1999-2004 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package pal.substmodel;

/**
 * <p>Title: RateMatrixHandler </p>
 * <p>Description: A utility class to manage the new style of rate matrices </p>
 * @author Matthew Goode
 * @version 1.0
 * <ul>
 *  <li> 11 May 2004 - Created, still a work in progress </li>
 *  <li> 2009-08-03 modified by M Woodhams (MDW). The parameterization now
 *  (optionally) allows the base frequencies to be optimised  also.
 * </ul>
 * 
 * 
 */
import pal.math.EigenSystem;
import pal.math.StochasticVector;
import pal.misc.*;

import java.io.PrintStream;
import java.io.PrintWriter;
public class RateMatrixHandler implements Parameterized, java.io.Serializable {
	private final NeoRateMatrix rateMatrix_;
	private final double[] parameters_;
	private final double[] parametersSE_;
	private final double[] defaultParameters_;
	private boolean updateMatrix_ = true;
	private final double[] equilibriumFrequencies_;
	private final int dimension_;
	private final boolean fixedFrequencies_;
	private final int nParameters_;
	private final int firstFrequencyParameter_;

	private final double[][] relativeRateStore_;
	private final double[][] qMatrixStore_;

	//private final MatrixExponential matrixExp_;
	private final EigenSystem eigenSystem_;

	private final boolean reversible_;
	
	private RateMatrixHandler(RateMatrixHandler toCopy) {
		this(toCopy.rateMatrix_, toCopy.equilibriumFrequencies_, toCopy.fixedFrequencies_);
		System.arraycopy(toCopy.parameters_, 0, parameters_,0,nParameters_);
		System.arraycopy(toCopy.parametersSE_, 0, parametersSE_,0,nParameters_);
		System.arraycopy(toCopy.defaultParameters_, 0, defaultParameters_,0,nParameters_);
	}

	// Simplified constructor for the common case: if we specified the equilibrium frequencies,
	// chances are they are variable.
	public RateMatrixHandler(NeoRateMatrix rateMatrix, double[] equilibriumFrequencies) {
		this(rateMatrix,equilibriumFrequencies,false);
	}
	/*
	 *  Simplified constructor for the common case of equal base frequencies
	 *  OR a NeoRateMatrix which is not independent of equilibrium frequencies
	 *  (e.g. GeneralLinearRateMatrix.)
	 */
	public RateMatrixHandler(NeoRateMatrix rateMatrix) {
		this(rateMatrix,rateMatrix.isIndependentOfEqbmFreq() ? equalBaseFrequencies(rateMatrix) : null, true);
	}
	// Helper method for RateMatrixHandler(NeoRateMatrix)
	private static double[] equalBaseFrequencies(NeoRateMatrix rateMatrix) {
		int dimension = rateMatrix.getDimension();
		double[] baseFrequencies = new double[dimension];
		for (int i=0; i<dimension; i++) {
			baseFrequencies[i] = 1./dimension;
		}
		return baseFrequencies;
	}
	
	// The full constructor
	public RateMatrixHandler(NeoRateMatrix rateMatrix, double[] equilibriumFrequencies, boolean fixedFrequencies) {
			
		if (rateMatrix.isIndependentOfEqbmFreq()) {
			if (equilibriumFrequencies == null) 
				throw new IllegalArgumentException("Can only pass null equilibrium frequences for a non-independent-equilibrium-frequencies rate matrix");
		} else {
			if (!fixedFrequencies)
				throw new IllegalArgumentException("Rate matrix not independent of equilibrium frequencies must have 'fixed frequencies'");
			// It wouldn't hurt to omit this check, but if you know what you're doing you'll never trip it.
			if (equilibriumFrequencies != null)
				throw new IllegalArgumentException("Equilibrium frequences should be null for rate matrix not independent of equilibrium frequencies");
		}
		this.rateMatrix_ = rateMatrix;
		this.dimension_ = rateMatrix_.getDimension();
		this.fixedFrequencies_ = fixedFrequencies;
		this.firstFrequencyParameter_ = rateMatrix.getNumberOfRateParameters();
		if (fixedFrequencies) {
			this.nParameters_ = firstFrequencyParameter_;
		} else {
			this.nParameters_ = firstFrequencyParameter_ + dimension_ - 1;
		}
		this.parameters_ = new double[nParameters_];
		this.parametersSE_ = new double[parameters_.length];
		this.defaultParameters_ = new double[parameters_.length];
		rateMatrix.getDefaultRateParameters(defaultParameters_,0);
		// Set the supplied frequencies as the default
		if (!fixedFrequencies) {
			StochasticVector.makeParameters(equilibriumFrequencies.length, equilibriumFrequencies,0,defaultParameters_,firstFrequencyParameter_);
		} 
		System.arraycopy(defaultParameters_,0,parameters_,0,parameters_.length);
		if (rateMatrix_.isIndependentOfEqbmFreq()) {
			this.equilibriumFrequencies_ = pal.misc.Utils.getCopy(equilibriumFrequencies);
		} else {
			this.equilibriumFrequencies_ = new double[dimension_];
		}
		this.relativeRateStore_ = new double[dimension_][dimension_];
		this.qMatrixStore_ = new double[dimension_][dimension_];
		eigenSystem_ = new EigenSystem();
		this.reversible_ = rateMatrix_.isReversible();
  }

	public final RateMatrixHandler getCopy() { return new RateMatrixHandler(this); }

	public final double[] getEquilibriumFrequencies() { 
		if (!fixedFrequencies_ || rateMatrix_.isIndependentOfEqbmFreq()) {
			checkMatrix();
		}
		return equilibriumFrequencies_; 
	}
	
	/**
	 * Does this rate matrix have fixed equilibrium frequencies (e.g. Jukes Cantor)
	 * or variable (e.g. HKY, GTR)?
	 * @return
	 */
	public final boolean getFixedFrequencies() {
		return fixedFrequencies_;
	}
	
	// MDW: added new method, required for MosaicModelEdgeParameters, MosaicSubstitutionModel
	public void getRateMatrix(double[][] matrixStore) {
		checkMatrix();
		pal.misc.Utils.copy(qMatrixStore_, matrixStore);
	}

	private final void checkMatrix() {
		if(updateMatrix_) {
			if (rateMatrix_.isIndependentOfEqbmFreq()) {
				rateMatrix_.createRelativeRates(relativeRateStore_,parameters_,0);
				if (!fixedFrequencies_) {
					StochasticVector.makeVector(dimension_, parameters_, firstFrequencyParameter_, equilibriumFrequencies_, 0);
				}
				fromQToR(relativeRateStore_,equilibriumFrequencies_,qMatrixStore_,dimension_,reversible_);
				double scale = makeValid(qMatrixStore_,equilibriumFrequencies_,dimension_);
				scale(qMatrixStore_,dimension_,scale);
				eigenSystem_.setNewLogMatrix(qMatrixStore_);
			} else {
				rateMatrix_.createRelativeRates(qMatrixStore_,parameters_,0);
				double sum=0;
				for (int i=0; i<dimension_; i++) sum -= qMatrixStore_[i][i];
				scale(qMatrixStore_,dimension_,sum/dimension_);
				eigenSystem_.setNewLogMatrix(qMatrixStore_);
				/*
				// Find the equilibrium frequencies as the eigenvector corresponding to the eigenvalue of 1.
				double[] eVals = eigenSystem_.getRealEigenvalues();
				double bestDiff = 100;
				int bestIndex = -1;
				for (int i=0; i<dimension_; i++) {
					double diff = Math.abs(eVals[i]-1);
					if (diff < bestDiff) {
						bestDiff = diff;
						bestIndex = i;
					}
				}
				double[][] eVecs = eigenSystem_.getEigenVectorsA();
				for (int i=0; i<dimension_; i++) equilibriumFrequencies_[i]=eVecs[i][bestIndex];
				*/
				equilibriumFromQ();
			}
			updateMatrix_ = false;
		}
	}
	
	/**
	 * Use the Q matrix to determine the equilibrium frequencies.
	 */
	/* 
	 * If pi is the equilibrium frequency vector, then we have
	 * pi Q = 0
	 * (This is just the detailed balance requirement.)
	 * So pi (Q+I) = pi, so pi is an eigenvector of Q+I 
	 */
	// TODO: make qplusI, eigen into data members so that we don't have to keep creating new objects
	private void equilibriumFromQ() {
		double[][] qplusI = new double[dimension_][dimension_];
		for (int i=0; i<dimension_; i++) {
			for (int j=0; j<dimension_; j++) {
				qplusI[j][i] = qMatrixStore_[i][j] + ((i==j) ? 1 : 0);
			}
		}
		EigenSystem eigen = new EigenSystem(qplusI);
		double[] eVals = eigen.getRealEigenvalues();
		double bestDiff = 100;
		int bestIndex = -1;
		for (int i=0; i<dimension_; i++) {
			double diff = Math.abs(eVals[i]-1);
			if (diff < bestDiff) {
				bestDiff = diff;
				bestIndex = i;
			}
		}
		/*
		 * To find the equilibrium vector pi: 
		 * Markov matrix P(t) = exp(Qt) = U^-1 exp(D t) U
		 * where U = matrix with rows = eigenvectors, D = diagonal matrix of eigenvalues.
		 * Given first e-value is 1, rest are 0<x<1, for t-> infinity we end up with
		 * P(t) -> matrix where every column is the first column of U^-1, and this
		 * column is the vector pi.
		 * Alternatively: pi Q = 0 (by detailed balance)
		 * so pi (Q+I) = pi, so can find eigenvectors of Q+I. (I didn't get this to work)
		 * However, given we have the eigensystem for Q already, easiest thing is to 
		 * throw a long branch length at it. 
		 */
		eigenSystem_.getMatrixPowerA(50, relativeRateStore_);
		for (int i=0; i<dimension_; i++) equilibriumFrequencies_[i]=relativeRateStore_[0][i];
	}
	
	
	public void getTransitionProbabilities(double distance, double[][] store ) {
		checkMatrix();
		eigenSystem_.getMatrixPowerA(distance, store);
	}
	public void getTransitionProbabilitiesTranspose(double distance, double[][] store ) {
		checkMatrix();
		eigenSystem_.getMatrixPowerTransposeA(distance, store);
	}


	// interface Report (remains abstract)

	// interface Parameterized (remains abstract)

	/** Computes normalized rate matrix from Q matrix (general reversible model)
	 * - Q_ii = 0
	 * - Q_ij = Q_ji
	 * - Q_ij is stored in R_ij (rate)
	 * - only upper triangular is used
	 * Also updates related MatrixExponential
	 */
	private static final void fromQToR(double[][] relativeRates, double[] equilibriumFrequencies, double[][] qMatrix, int dimension, boolean reversible) {
		if(reversible) {
			for( int i = 0; i<dimension; i++ ) {
				for( int j = i+1; j<dimension; j++ ) {
					qMatrix[i][j] = relativeRates[i][j]*equilibriumFrequencies[j];
					qMatrix[j][i] = relativeRates[i][j]*equilibriumFrequencies[i];
				}
			}
		} else {
		  for( int i = 0; i<dimension; i++ ) {
				for( int j = i+1; j<dimension; j++ ) {
					qMatrix[i][j] = relativeRates[i][j]*equilibriumFrequencies[j];
					//This is the only difference
					qMatrix[j][i] = relativeRates[j][i]*equilibriumFrequencies[i];
				}
			}
		}
	}

	//
	// Private stuff
	//

	/** Make it a valid rate matrix (make sum of rows = 0)
		* @return current rate scale
		*/
	private static final  double makeValid(double[][] relativeRates, double[] equilibriumFrequencies, int dimension) {
		double total = 0;
		for (int i = 0; i < dimension ; i++){
			double sum = 0.0;
			for (int j = 0; j < dimension ; j++)	{
				if (i != j)	{
					sum += relativeRates[i][j];
				}
			}
			relativeRates[i][i] = -sum;
			total+=equilibriumFrequencies[i]*sum;
		 }
		 return total;
	}
	private final static double calculateNormalScale(double[][] relativeRates, double[] equilibriumFrequencies, int dimension) {
	  double scale = 0.0;

		for (int i = 0; i < dimension; i++)	{
			scale += -relativeRates[i][i]*equilibriumFrequencies[i];
		}
		return scale;
	}

	// Normalize rate matrix to one expected substitution per unit time
	private static final void normalize(final double[][] relativeRates, double[] equilibriumFrequencies, int dimension) {
		scale(relativeRates,dimension, calculateNormalScale(relativeRates,equilibriumFrequencies,dimension));
	}
	 // Normalize rate matrix by a certain scale to achieve an overall scale (used with a complex site class model)
	private static final void scale(final double[][] relativeRates,int dimension, double scale)  {
		for (int i = 0; i < dimension; i++)  {
			for (int j = 0; j < dimension; j++)  {
				relativeRates[i][j] = relativeRates[i][j]/scale;
			}
		}
	}

	/**
	 * Reporting stuff
	 * @param out where to report too
	 */
	public void report(PrintWriter out) {
	  out.println("RateMatrixHandler: parameters = [" + pal.misc.Utils.toString(parameters_) + "], rate matrix:");
	  out.println("Reporting Not functioning yet...");
	}

	private final void parametersChanged() { this.updateMatrix_ = true; }

	public int getNumParameters() { 	return parameters_.length;	}
	public void setParameter(double param, int n) { this.parameters_[n] = param; parametersChanged(); }
	public double getParameter(int n) { return parameters_[n]; }
	public void setParameterSE(double paramSE, int n) { this.parametersSE_[n] = paramSE; }
	/*
	 *  Frequency parameters are always in the range 0-1. (See StochasticVector class.)
	 * However, checkMatrix() method dies if the first frequency parameter = 1 
	 * (gives scale = 0) so avoid this. The part of parameter space excluded is
	 * never going to matter in real world applications.
	 * 
	 * Also: If any of the frequency parameters are equal to zero, the corresponding base frequency
	 * is zero, which is liable to cause non-numeric likelihood evaluations and mess up optimizations.
	 * So set the lower limit to a small number rather than zero.
	 */
	public double getLowerLimit(int n) { 
		return (n >= firstFrequencyParameter_) ? 1e-6 : rateMatrix_.getRateParameterLowerBound(n); 
	}
	public double getUpperLimit(int n) { 
		return (n >= firstFrequencyParameter_) ? 0.9999 : rateMatrix_.getRateParameterUpperBound(n); 
	}
	public double getDefaultValue(int n) { return defaultParameters_[n]; }
	
	// MDW added for easier debugging. Possibly this is what 'report()' is supposed to do.
	public void printObject() { this.printObject(System.out, ""); }
	public void printObject(String ps) {this.printObject(System.out, ps); }
	public void printObject(PrintStream out, String ps) {
		out.println(ps+"Object: RateMatrixHandler");
		rateMatrix_.printObject(out, ps+"  ",parameters_);
		System.out.println(ps+"unfinished...");
	}
	
	/**
	 * test:
	 * Generate a Markov matrix, check that it gives correct result (i.e. the same as it used to.)
	 * Run this if you make substantive changes to the class. Add new test cases as needed.
	 */
	public static void test() {
		PrintWriter output = new PrintWriter(System.out, true);

		
		double[] equalibriumFreq = new double[]{0.1,0.2,0.3,0.4};
		RateMatrixHandler rmh1t = new RateMatrixHandler(GeneralREVRateMatrix.createGTR(),equalibriumFreq,true);
		RateMatrixHandler rmh1f = new RateMatrixHandler(GeneralREVRateMatrix.createGTR(),equalibriumFreq,false);
		RateMatrixHandler rmh2t = new RateMatrixHandler(new GeneralPoissonRateMatrix(4),equalibriumFreq,true);
		RateMatrixHandler rmh2f = new RateMatrixHandler(new GeneralPoissonRateMatrix(4),equalibriumFreq,false);
		 
		double[] GTRparameters = new double[]{1.5,0.7,1.2,0.5,1.8};
		for (int i=0; i<GTRparameters.length; i++) {
			rmh1t.setParameter(GTRparameters[i],i);
			rmh1f.setParameter(GTRparameters[i],i);
		}
		 
		double[][] gtrAnswers = new double[][]
		        {{0.4393, 0.1533, 0.1408, 0.2666},
				 {0.0767, 0.4653, 0.1228, 0.3353},
				 {0.0469, 0.0818, 0.6407, 0.2306},
				 {0.0667, 0.1676, 0.1729, 0.5928}};
		double[][] poissonAnswers = new double[][]
		        {{0.4311, 0.1264, 0.1896, 0.2528},
				 {0.0632, 0.4943, 0.1896, 0.2528},
				 {0.0632, 0.1264, 0.5575, 0.2528},
				 {0.0632, 0.1264, 0.1896, 0.6207}};
		                                          
		 
		double[][] markov = new double[4][4];
		double branchLength = 0.7;
		 
		 
		rmh1t.getTransitionProbabilities(branchLength, markov);
		/*
		for (int i=0; i<4; i++) output.printf("[%6.4f %6.4f %6.4f %6.4f]\n", markov[i][0],markov[i][1],markov[i][2],markov[i][3]);
		output.println();
		*/
		output.println(matrixInfinityNorm(markov,gtrAnswers)<0.0008 ? "Passed" : "Failed");

		rmh1f.getTransitionProbabilities(branchLength, markov);
		output.println(matrixInfinityNorm(markov,gtrAnswers)<0.0008 ? "Passed" : "Failed");
		 
		rmh2t.getTransitionProbabilities(branchLength, markov);
		output.println(matrixInfinityNorm(markov,poissonAnswers)<0.0008 ? "Passed" : "Failed");	
		 
		rmh2f.getTransitionProbabilities(branchLength, markov);
		output.println(matrixInfinityNorm(markov,poissonAnswers)<0.0008 ? "Passed" : "Failed");
		 
		RateMatrixHandler rmh3 = new RateMatrixHandler(GeneralLinearRateMatrix.getK3ST_F81());
		/*
		double[] parameters = new double[]{1.5,0.7,1.2,0.5,1.8};
		for (int i=0; i<parameters.length; i++) rmh3.setParameter(parameters[i],i);
		rmh3.getTransitionProbabilities(branchLength, markov);
		for (int i=0; i<4; i++) output.printf("[%6.4f %6.4f %6.4f %6.4f]\n", markov[i][0],markov[i][1],markov[i][2],markov[i][3]);
		output.println();
		output.println("Eigenvectors:");
		rmh3.eigenSystem_.getEigenVectorsM().print(output, 8, 4);
		output.printf("Real part of Eigenvalues: %s\n", pal.misc.Utils.toString(rmh3.eigenSystem_.getRealEigenvalues()));
		output.printf("Imag part of Eigenvalues: %s\n", pal.misc.Utils.toString(rmh3.eigenSystem_.getImagEigenvalues()));
		*/
		/*
		 *  This example was generated in Mathematica. 
		 *
		 * Need multiple rescalings to get agreement.
		 * First, GeneralLinearRateMatrix.getK3ST_F81() takes 6th parameter as =1, and 
		 * treats the rest as being relative to this - so parameters need scaling by 6th one.
		 * Second, RMH scales the Q matrix to have trace = - dimension.
		 * The trace of the raw Q matrix is 10*(sum of parameters + 1) so this will
		 * multiply Q matrix by .4/(sum of parameters + 1).
		 * To counteract these, need to scale branch length inversely:
		 * PAL branch length = Mathematica branch length * rawParameter[6] * (sum of parameters+1)*2.5
		 * have 6th parameter = 1.
		 */
		
		double[][] K3ST_F81Answers = new double[][]
		          {{0.200099, 0.244414, 0.266634, 0.288853},
		           {0.199975, 0.244538, 0.266634, 0.288853},
		           {0.199975, 0.244414, 0.266757, 0.288853},
		           {0.199975, 0.244414, 0.266634, 0.288977}};
		testGLRM(1.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, markov);
		output.println(matrixInfinityNorm(markov,K3ST_F81Answers)<0.000008 ? "Passed" : "Failed");
		
		K3ST_F81Answers = new double[][]
		          {{0.239751, 0.293333, 0.228351, 0.238565},
		           {0.239733, 0.293350, 0.228347, 0.238570},
		           {0.239723, 0.293318, 0.228379, 0.238580},
		           {0.239718, 0.293322, 0.228362, 0.238598}};
		testGLRM(1, 0.856775, 0.353035, 0.266428, 0.510916, 0.695405, 0.524003, markov);
		output.println(matrixInfinityNorm(markov,K3ST_F81Answers)<0.000008 ? "Passed" : "Failed");		
		output.println("Finished RateMatrixHandler test");
	}
	/*
	 * This gets around all the annoying scaling that GeneralLinearRateMatrix does.
	 * Minor point: order of parameters in my test set differs from my test set. 
	 * Scaling 1: the final parameter (c') is set to be 1 by K3ST_F81, but is free parameter in the test data
	 * Scaling 2: RMH scales the Q matrix to have trace of -dimension (i.e. -4).
	 * This routine compensates for those scalings, and returns results in 'markov'.
	 */
	private static void testGLRM(double t, double a, double b, double c, double cp, double bp, double ap, double[][] markov) {
		RateMatrixHandler rmh = new RateMatrixHandler(GeneralLinearRateMatrix.getK3ST_F81());
		double[] unscaledParams = new double[]{a,ap,b,bp,c,cp};
		double lastParam = unscaledParams[unscaledParams.length-1];
		double sum = 0;
		for (int i=0; i<unscaledParams.length-1; i++) {
			double scaledParam = unscaledParams[i]/lastParam;
			rmh.setParameter(scaledParam,i);
			sum += scaledParam;
		}
		double qScale = 1/(lastParam * (sum+1) * 2.5); // to convert back to raw Q matrix, from original parameters		
		double scaledT = t/qScale;
		rmh.getTransitionProbabilities(scaledT, markov);
	}
	private static double matrixInfinityNorm(double[][]A, double[][]B) {
		double sum = 0;
		for (int i=0; i<A.length; i++) {
			for (int j=0; j<A[0].length; j++) {
				sum += Math.abs(A[i][j]-B[i][j]);
			}
		}
		return sum;
	}
}