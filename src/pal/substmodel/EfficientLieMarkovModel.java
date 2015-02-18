package pal.substmodel;

import java.io.PrintStream;

import pal.datatype.DataType;
import pal.math.MathUtils;
import pal.math.StochasticVector;

/**
 * The models from this class are wrappers around a LieMarkovModel, but
 * with two redundant parameters eliminated.
 *
 * @author Michael
 *
 */

/*
 * TODO:
 * I'm not happy with the structure here. I extend LieMarkovModel so as to have access
 * to the private ray families etc, but that is a static class and this one is not.
 * Perhaps change LieMarkovModel?
 */

// TODO: 

@SuppressWarnings("serial")
public class EfficientLieMarkovModel extends LieMarkovModel implements NeoRateMatrix {
	private GeneralLinearRateMatrix lmm_;
	private int firstFamilySize_;
	private int stochasticVectorStart_; // = firstFamilySize_ - 1
	private int stochasticVectorLen_; // = basis length - firstFamilySize
	private double[] newParams_;
	private double[] lowerBounds_;
	private double[] upperBounds_;
	private String uniqueName_;
	private double[] lmmParameters_; // preallocated workspace, need not be persistent.
	
	// Same bound used in GeneralLinearRateMatrix
	private static final double STOCHASTIC_UPPER_BOUND = 0.999999; // just short of 1. (Actually =1 might cause problems.)
	private static final double DEFAULT_BOUND = 100;
	
	/*
	 *  Can't call it 'getModel6c' as return type conflicts with LieMarkovModel.getModel6c.
	 *  This part of the 'I haven't got the class structure right' problem mentioned above.
	 */
	public static EfficientLieMarkovModel getEfficient6_7a() {
		// must put rayFam412a first, for first family weighting to work.
		// (Alternative is to have a way to specify firstFamilySize explicitly to = 3,
		// so (rayFam140, rayFam204a) together get first family weighting scheme.)
		return new EfficientLieMarkovModel(new double[][][][] {orbit412a,orbit110,orbit201a}, "6.7a");
	}
	public static EfficientLieMarkovModel getEfficient10_12() {
		// must put one of families 15 and 21 first, as 14 does not have the
		// 'add a constant to all weights' redundancy with either of the other two families.
		return new EfficientLieMarkovModel(new double[][][][] {orbit401c,orbit401e,orbit410a}, "10.12");
	}
	
	private EfficientLieMarkovModel(double[][][][] families, String name) {
		firstFamilySize_ = families[0].length;

		uniqueName_ = name;
		if (firstFamilySize_ <= 1) throw new IllegalArgumentException("First ray family must be size 2 or bigger");
		int nBasis = 0;
		for (double[][][] family : families) nBasis += family.length;
		stochasticVectorStart_ = firstFamilySize_ - 1;
		stochasticVectorLen_ = nBasis - firstFamilySize_;
		double[][][] basis = new double[nBasis][][];
		int i=0;
		for (double[][][] family : families)
			for (double[][] ray : family)
				basis[i++] = ray;
		int dimension = basis[0][0].length;
		lmm_ = new GeneralLinearRateMatrix(dimension, null, basis, false);
		int nNewParams = nBasis - 2;
		newParams_ = new double[nNewParams];
		lowerBounds_ = new double[nNewParams];
		upperBounds_ = new double[nNewParams];
		this.getDefaultRateParameters(newParams_, 0);
		// Set the bounds for the stochastic vector parameters:
		for (i=stochasticVectorStart_; i<lowerBounds_.length; i++) {
			lowerBounds_[i] = 0;
			upperBounds_[i] = STOCHASTIC_UPPER_BOUND;
		}
		// and set bounds for first family parameters:
		this.setBounds(DEFAULT_BOUND);
		lmmParameters_ = new double[nBasis]; 
	}
	
	// NeoRateMatrix methods:
	
	// This is where we translate the parameters. This is the only place we use preallocated array lmmParameters_.
	public void createRelativeRates(double[][] rateStore, double[] rateParameters, int startIndex) {
		// Holds the minimum of the first family parameters and zero
		double min = MathUtils.min(0,MathUtils.getMinimum(rateParameters, 0, firstFamilySize_-2));
		for (int i=0; i<firstFamilySize_-1; i++) lmmParameters_[i] = rateParameters[i]-min;
		lmmParameters_[firstFamilySize_-1] = -min;
		// And StochasticVector treatment for the rest
		StochasticVector.makeVector(stochasticVectorLen_, rateParameters, stochasticVectorStart_, lmmParameters_, stochasticVectorStart_+1);
	
		lmm_.createRelativeRates(rateStore, lmmParameters_, startIndex);
}
	
	public String getUniqueName() { return uniqueName_; }
	public boolean isReversible() {	return lmm_.isReversible(); } // false
	public boolean isIndependentOfEqbmFreq() {return lmm_.isIndependentOfEqbmFreq(); } // false
	public int getDimension() { return lmm_.getDimension(); } // 4
	public boolean isDataTypeCompatible(DataType dt) { return lmm_.isDataTypeCompatible(dt); }

	public int getNumberOfRateParameters() { return newParams_.length; }
	public double getRateParameterLowerBound(int parameter) { return lowerBounds_[parameter]; }
	public double getRateParameterUpperBound(int parameter) { return upperBounds_[parameter]; }

	public void getDefaultRateParameters(double[] parameterStore, int startIndex) {
		// For the 'first family', we have firstFamilySize-1 parameters. (A savings of 1,
		// due to the constraint that one of the first family rays will have weight 0.)
		// Default is all equal zero. (Could reconsider this?)
		for (int i=0; i<firstFamilySize_-1; i++) parameterStore[startIndex+i]=0;
		// For the remainder, we use a StochasticVector to impose the constraint
		// that all of them are non-negative, and they sum to 1.
		StochasticVector.getFlatParameters(stochasticVectorLen_, parameterStore, stochasticVectorStart_);
	}

	/**
	 * Note: StochasticVector parameters are NOT changed by this,
	 * only the 'first family' parameters.
	 * 
	 * Best to use setBounds(double) instead.
	 */
	public void setUpperBounds(double bound) {
		if (bound <=0) throw new IllegalArgumentException("Bad upper bound");
		for (int i=0; i<firstFamilySize_-1; i++) upperBounds_[i]=bound;
	}

	/**
	 * Note: StochasticVector parameters are NOT changed by this.
	 * Any elements in bounds beyond those needed for the 'first family'
	 * parameters are ignored.
	 * 
	 * Best to use setBounds(double) instead.
	 */
	public void setUpperBounds(double[] bounds) {
		// 
		for (int i=0; i<firstFamilySize_-1; i++) {
			if (bounds[i] <=0) throw new IllegalArgumentException("Bad upper bound");
			upperBounds_[i]=bounds[i];
		}
	}

	/**
	 * Best to use setBounds(double) instead
	 * @param bound
	 */
	public void setLowerBounds(double bound) {
		if (bound ==0) throw new IllegalArgumentException("Bad lower bound");
		bound = -Math.abs(bound); // make it negative if it isn't already
		for (int i=0; i<firstFamilySize_-1; i++) lowerBounds_[i]=bound;
	}
	
	/**
	 * For first family, set bounds to +/- bound.
	 * @param bound
	 */
	public void setBounds(double bound) {
		setLowerBounds(-bound);
		setUpperBounds(bound);
	}

	public void printObject(PrintStream out, String ps, double[] parameters) {
		this.printObject(out, ps);
		out.println("Unfinished...");
	}
	public void printObject(PrintStream out) {this.printObject(out, "");}
	public void printObject(PrintStream out, String prefixString) {
		out.printf(prefixString+"Object: EfficientLieMarkovModel, firstFamilySize = %d, number of parameters = %d, wrapper for the following:\n",firstFamilySize_,newParams_.length);
		lmm_.printObject(out, prefixString + ">");
	}
}
