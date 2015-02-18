package pal.substmodel;

import java.io.PrintStream;

import pal.datatype.DataType;
import pal.math.StochasticVector;

/**
 * This models any rate matrix where the entries are linear in the parameters.
 * (We do allow a constant to be added to each rate matrix entry, violating
 * the strict mathematical definition of linearity. This is so that we can
 * in effect constraint the scaling of the matrix, by requiring that one of
 * the parameters equal 1.)
 */

@SuppressWarnings("serial")
public class GeneralLinearRateMatrix implements NeoRateMatrix {
	private static final double DEFAULT_UPPER_BOUND = 100;
	private static final double DEFAULT_LOWER_BOUND = 0;
	private static final double STOCHASTIC_UPPER_BOUND = 0.9999; // just short of 1. (Actually 1 might cause problems.)
	
	private final double[][][] basis_;
	private final int nParameters_;
	private final double[][] constant_;
	private final double[] defaultParameters_;
	private final int dimension_;
	private final double[] upperBounds_;
	private final double[] lowerBounds_; // allowing non-zero here is experimental and may get reversed.
	private final boolean useStochasticVectorParameters_;
	private final double[] weights_; // preallocated mem space when useStochasticVectorParameters_ = true.
	/*
	 * If useStochasticVectorParameters_, then the parameters are treated as a StochasticVector,
	 * so n parameters in range 0-1 get expanded to n+1 nonnegative weights summing to one. This
	 * constrains the overall scale of the Q matrix, avoiding degenerate overparameterization, where
	 * Q matrix scale and overall scale of tree branch lengths are degenerate.
	 * Using the 'constant_' matrix is another way of constraining the scale of the Q matrix: 
	 * essentially this constrains one of the parameters to be equal to 1. 
	 */
	

	/**
	 * 
	 * @param dimension
	 * @param constant may be null, else double[dimension][dimension]
	 * @param basis double[nParameters][dimension][dimension]
	 */
	public GeneralLinearRateMatrix(int dimension, double[][] constant, double[][][] basis, boolean weightsSumToOne) {
		
		dimension_ = dimension;
		constant_ = constant;
		basis_ = basis;
		if (weightsSumToOne) {
			useStochasticVectorParameters_ = true;
			nParameters_ = basis.length-1; 
			weights_ = new double[basis.length];
		} else {
			useStochasticVectorParameters_ = false;
			nParameters_ = basis.length; 
			weights_ = null;
		}
		defaultParameters_ = new double[nParameters_];
		upperBounds_ = new double[nParameters_];
		lowerBounds_ = new double[nParameters_];
		if (weightsSumToOne) {
			setUpperBounds(STOCHASTIC_UPPER_BOUND);
			setLowerBounds(0.0);
			for (int i=0; i<nParameters_; i++) {
				defaultParameters_[i] = 1.0/(basis.length-i);
			}
		} else {
			setDefaultParameters(1);
			setUpperBounds(DEFAULT_UPPER_BOUND);
			setLowerBounds(DEFAULT_LOWER_BOUND);
		}
	}
	
	/**
	 * Return a new GeneralLinearRateMatrix created by permuting the rows/columns of an existing one
	 * @param permutation
	 * @return
	 */
	public GeneralLinearRateMatrix(GeneralLinearRateMatrix base, int[] permutation) {
		if (permutation.length != base.dimension_) throw new IllegalArgumentException("Permutation wrong length");
		nParameters_ = base.nParameters_;
		defaultParameters_ = base.defaultParameters_.clone();
		upperBounds_ = base.upperBounds_.clone();
		lowerBounds_ = base.lowerBounds_.clone();
		useStochasticVectorParameters_ = base.useStochasticVectorParameters_;
		if (useStochasticVectorParameters_) {
			weights_ = new double[nParameters_+1];
		} else {
			weights_ = null;
		}
		
		dimension_ = base.dimension_;
		constant_ = (base.constant_==null) ? null : new double[dimension_][dimension_];
		basis_ = new double[base.basis_.length][dimension_][dimension_];
		for (int row=0; row<dimension_; row++) {
			int pRow = permutation[row];
			for (int col=0; col<dimension_; col++) {
				int pCol = permutation[col];
				if (constant_!=null) {
					constant_[row][col] = base.constant_[pRow][pCol];
				}
				for (int i=0; i<basis_.length; i++) {
					basis_[i][row][col] = base.basis_[i][pRow][pCol];
				}
			}
		}
	}
	
	/*
	 * This is the heart of the class:
	 */
	public void createRelativeRates(double[][] rateStore, double[] rateParameters, int startIndex) {
		int n = nParameters_;
		if (useStochasticVectorParameters_) {
			n++;
			StochasticVector.makeVector(n, rateParameters, startIndex, weights_, 0);
			rateParameters = weights_;
			startIndex = 0;
		}

		// Indices: r = row, c = column, p = parameter
		if (constant_ == null) {
			for (int r=0; r<dimension_; r++) 
				for (int c=0; c<dimension_; c++)
					rateStore[r][c] = 0;
		} else {
			for (int r=0; r<dimension_; r++) 
				for (int c=0; c<dimension_; c++)
					rateStore[r][c] = constant_[r][c];
		}
		for (int p=0; p<n; p++) {
			double param = rateParameters[p+startIndex];
			for (int r=0; r<dimension_; r++) 
				for (int c=0; c<dimension_; c++)
					rateStore[r][c] += param*basis_[p][r][c];
		}
	}

	/*
	 * I'm not sure that setUpperBounds will ever get used
	 * bounds should have length nParameters_.
	 */
	public void setUpperBounds(double[] bounds) {
		if (useStochasticVectorParameters_) throw new RuntimeException("Can't alter parameter bounds when using weights sum to one");
		System.arraycopy(bounds,0,upperBounds_,0,nParameters_);
	};
	public void setUpperBounds(double bound) {
		if (useStochasticVectorParameters_ && bound != STOCHASTIC_UPPER_BOUND) 
			throw new RuntimeException("Can't alter parameter bounds when using weights sum to one");
		for (int i=0; i<nParameters_; i++) upperBounds_[i] = bound;
	}
	public void setLowerBounds(double bound) {
		if (useStochasticVectorParameters_ && bound != 0) throw new RuntimeException("Can't alter parameter bounds when using weights sum to one");
		for (int i=0; i<nParameters_; i++) lowerBounds_[i] = bound;
	}

	
	public void setDefaultParameters(double[] defaultParams) {
		System.arraycopy(defaultParams,0,defaultParameters_,0,nParameters_);
	};
	public void setDefaultParameters(double defaultParam) {
		for (int i=0; i<nParameters_; i++) defaultParameters_[i] = defaultParam;
	}

	
	public String getUniqueName() { return "General Linear Rate Matrix (dimension "+dimension_+", basis size "+nParameters_+")"; }
	// TODO: (possibly) can be reversible in special circumstances. Account for this?
	public boolean isReversible() { return false; }
	// Parameters determine equilibrium freq through complex formula - can't set eqbm freqs. independently.
	public boolean isIndependentOfEqbmFreq() { return false; }
	public int getDimension() { return dimension_; }
	public boolean isDataTypeCompatible(DataType dt) { return dt.getNumStates()==dimension_; }

	

	public int getNumberOfRateParameters() { return nParameters_; }
	public double getRateParameterLowerBound(int parameter) { return lowerBounds_[parameter]; }
	public double getRateParameterUpperBound(int parameter) { return upperBounds_[parameter]; }
	public void getDefaultRateParameters(double[] store, int startIndex) {
		System.arraycopy(defaultParameters_,0,store,startIndex,defaultParameters_.length);
	}

	public void printObject(PrintStream out) {this.printObject(out, "");}
	public void printObject(PrintStream out, String ps) {
		out.printf(ps+"Object: GeneralLinearRateMatrix, dimension = %d, number of parameters = %d",dimension_,nParameters_);
		if (useStochasticVectorParameters_) {
			out.println(" (stochastic vector weights)");
		} else {
			out.println();
		}
		if (constant_ == null) {
			out.println(ps+"No constant matrix");
		} else {
			out.println(ps+"Constant matrix:");
			for (int r = 0; r < dimension_; r++) {
				for (int c = 0; c < dimension_; c++) {
					out.printf("% 8f%s", constant_[r][c], (c+1 == dimension_) ? "\n"+ps : ", ");
				}
			}
		}
		for (int p=0; p<basis_.length; p++) {
			out.printf("%sBasis matrix %d:\n",ps,p);
			for (int r = 0; r < dimension_; r++) {
				for (int c = 0; c < dimension_; c++) {
					out.printf("% 8f%s", basis_[p][r][c], (c+1 == dimension_) ? "\n"+ps : ", ");
				}
			}	
		}
	}
	
	public void printObject(PrintStream out, String ps, double[] parameters) {
		this.printObject(out, ps);
		out.println("Unfinished...");
	}
	
	
	/*
	 * 
	 * All the following deprecated - still used by some test routines. 
	 * Use LieMarkovModel class objects, or construct them in a similar way.
	 * 
	 */
	
	/**
	 * Returns a GeneralLinearRateMatrix for the Lie algebra conforming
	 * K3ST+F81 model (6 free parameters, 5 after fixing scaling.)
	 * 
	 * See "Lie Markov Models" Fern ́andez-S ́anchez, Peter Jarvis and Jeremy Sumner
	 * (unpublished paper, as of the writing of this code.)
	 * 
	 * @deprecated use LieMarkovModel.
	 * @return
	 */
	private static final double[][] W12 = new double[][]
	    {{-2, 2, 0, 0},
		 { 2,-2, 0, 0},
		 { 1, 1,-3, 1},
		 { 1, 1, 1,-3}};
	private static final double[][] W34 = new double[][]
	    {{-3, 1, 1, 1},
	     { 1,-3, 1, 1},
	     { 0, 0,-2, 2},
	     { 0, 0, 2,-2}};
	private static final double[][] W13 = new double[][]
	    {{-2, 0, 2, 0},
	     { 1,-3, 1, 1},
	     { 2, 0,-2, 0},
	     { 1, 1, 1,-3}};
	private static final double[][] W24 = new double[][]
	    {{-3, 1, 1, 1},
	     { 0,-2, 0, 2},
	     { 1, 1,-3, 1},
	     { 0, 2, 0,-2}};
	private static final double[][] W14 = new double[][]
	    {{-2, 0, 0, 2},
	     { 1,-3, 1, 1},
	     { 1, 1,-3, 1},
	     { 2, 0, 0,-2}};
	private static final double[][] W23 = new double[][]
	    {{-3, 1, 1, 1},
	     { 0,-2, 2, 0},
	     { 0, 2,-2, 0},
	     { 1, 1, 1,-3}};
	@Deprecated
	public static GeneralLinearRateMatrix getK3ST_F81() {
		double[][] constant = W23;
		double[][][] basis = new double[][][] {W12,W34,W13,W24,W14};
		return new GeneralLinearRateMatrix(4, constant, basis, false); 
	}
	
	/*
	 * Return a 7 parameter K3ST+F81 model. This will be doubly redundant:
	 * overall scaling is redundant with a scaling of the branch lengths,
	 * and the 7 matrices span a space of only 6 dimensions - however,
	 * it allows us to limit parameters to be non-negative without omitting
	 * any part of the space, which the 6 parameter formulation cannot.
	 * 
	 * I've chosen scaling such that all these matrices have the same trace,
	 * hence the weights-sum-to-one will become a trace-sums-to -12 constraint.
	 */
	private static final double[][] R1 = new double[][]
	     {{ 0, 0, 0, 0},
	      { 4,-4, 0, 0},
	      { 4, 0,-4, 0},
	      { 4, 0, 0,-4}};	
	private static final double[][] R2= new double[][]
	     {{-4, 4, 0, 0},
	      { 0, 0, 0, 0},
	      { 0, 4,-4, 0},
	      { 0, 4, 0,-4}};	
	private static final double[][] R3= new double[][]
	     {{-4, 0, 4, 0},
	      { 0,-4, 4, 0},
	      { 0, 0, 0, 0},
	      { 0, 0, 4,-4}};	
	private static final double[][] R4= new double[][]
	     {{-4, 0, 0, 4},
	      { 0,-4, 0, 4},
	      { 0, 0,-4, 4},
	      { 0, 0, 0, 0}};
	private static final double[][] La= new double[][]
	     {{-3, 3, 0, 0},
	      { 3,-3, 0, 0},
	      { 0, 0,-3, 3},
	      { 0, 0, 3,-3}};	
	private static final double[][] Lb= new double[][]
	     {{-3, 0, 3, 0},
	      { 0,-3, 0, 3},
	      { 3, 0,-3, 0},
	      { 0, 3, 0,-3}};	
	private static final double[][] Lg= new double[][]
	     {{-3, 0, 0, 3},
	      { 0,-3, 3, 0},
	      { 0, 3,-3, 0},
	      { 3, 0, 0,-3}};
	// These two should be equivalent, (can produce exactly the same set of Q matrices, up to scaling.)
	// 'Stochastic' uses one fewer parameter and fixes overall scale of the Q matrix.
	@Deprecated
	public static GeneralLinearRateMatrix getRedundantK3ST_F81() {
		double[][][] basis = new double[][][] {R1,R2,R3,R4,La,Lb,Lg};
		return new GeneralLinearRateMatrix(4, null, basis, false); 
	}
	@Deprecated
	public static GeneralLinearRateMatrix getStochasticRedundantK3ST_F81() {
		double[][][] basis = new double[][][] {R1,R2,R3,R4,La,Lb,Lg};
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	// 10 parameter Lie Markov model which allows transitions/transversions to be distinguishable.
	private static final double[][] Q_CT_A = new double[][]
	     {{-0, 0, 0, 0},
	      { 1,-1, 0, 0},
	      { 0, 0,-0, 0},
	      { 1, 0, 0,-1}};
	private static final double[][] Q_CT_G= new double[][]
	     {{-0, 0, 0, 0},
	      { 0,-1, 1, 0},
	      { 0, 0,-0, 0},
	      { 0, 0, 1,-1}};
	private static final double[][] Q_AG_C = new double[][]
	     {{-1, 1, 0, 0},
	      { 0,-0, 0, 0},
	      { 0, 1,-1, 0},
	      { 0, 0, 0,-0}};	   
	private static final double[][] Q_AG_T = new double[][]
	     {{-1, 0, 0, 1},
	      { 0,-0, 0, 0},
	      { 0, 0,-1, 1},
	      { 0, 0, 0,-0}};	   
	private static final double[][] Q_C_A__T_G = new double[][]
	     {{-0, 0, 0, 0},
	      { 1,-1, 0, 0},
	      { 0, 0,-0, 0},
	      { 0, 0, 1,-1}};	                                      
	private static final double[][] Q_C_G__T_A = new double[][]
	     {{-0, 0, 0, 0},
	      { 0,-1, 1, 0},
	      { 0, 0,-0, 0},
	      { 1, 0, 0,-1}};	                                      
	private static final double[][] Q_A_C__G_T = new double[][]
	     {{-1, 1, 0, 0},
	      { 0,-0, 0, 0},
	      { 0, 0,-1, 1},
	      { 0, 0, 0,-0}};	                                      
	private static final double[][] Q_A_T__G_C = new double[][]
	     {{-1, 0, 0, 1},
	      { 0,-0, 0, 0},
	      { 0, 1,-1, 0},
	      { 0, 0, 0,-0}};	                                      
	private static final double[][] Q_C_T = new double[][]
	     {{-0, 0, 0, 0},
	      { 0,-2, 0, 2},
	      { 0, 0,-0, 0},
	      { 0, 0, 0,-0}};	                                      
	private static final double[][] Q_T_C = new double[][]
	     {{-0, 0, 0, 0},
		  { 0,-0, 0, 0},
	      { 0, 0,-0, 0},
	      { 0, 2, 0,-2}};	                                      
	private static final double[][] Q_G_A = new double[][]
	     {{-0, 0, 0, 0},
	      { 0,-0, 0, 0},
	      { 2, 0,-2, 0},
	      { 0, 0, 0,-0}};	                                      
	private static final double[][] Q_A_G = new double[][]
	     {{-2, 0, 2, 0},
	      { 0,-0, 0, 0},
	      { 0, 0,-0, 0},
	      { 0, 0, 0,-0}};
	@Deprecated
	public static GeneralLinearRateMatrix getStochasticRedundant10a() {
		double[][][] basis = new double[][][] {
				Q_CT_A, Q_CT_G, Q_AG_C, Q_AG_T, 
				Q_C_A__T_G, Q_A_C__G_T, Q_C_G__T_A, Q_A_T__G_C,
				Q_C_T, Q_T_C, Q_G_A, Q_A_G};
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
}
