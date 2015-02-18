package pal.substmodel;

import java.io.PrintWriter;

import pal.math.StochasticVector;

/**
 * Class for rate distribution where we can set any rates and probabilities we wish, 
 * restricted only by mathematical validity.
 * We can optionally demand that the rates are normalized (average rate = 1).
 * 
 * Parameterization:
 * First n-1 parameters are stochastic vector parameterisation of probability vector.
 * Next n parameters are the rates. NOTE: If rates are normalized, we really only have
 * n-1 free parameters here. I should figure out an appropriate (n-1) parameterisation.
 * For now, best not to let these be free parameters in an optimisation.
 * 
 * @author woodhams
 *
 */

@SuppressWarnings("serial")
public class ArbitraryRates extends RateDistribution {
	public static final double MAX_RATE=1000;
	
	private boolean normalizedRates_; // if true, inner product of 'rate' and 'probability' = 1
	private double[] rawRate_;        // only needed if normalizedRates_ == true, else is synonymous with 'rate'.
	private double[] probabilityParameterisation_; // StochasticVector parameterisation of 'probability'
	// Inherited:
	// public int numRates;
	// public double[] rate;
	// public double[] probability;
	// !! Why are those 'public'?
	
	public ArbitraryRates(double[] rate, double[] probability, boolean normalized) {
		super(rate.length);
		normalizedRates_ = normalized;
		this.probability = probability;
		rawRate_ = rate;
		probabilityParameterisation_ = new double[numRates-1];
		StochasticVector.makeParameters(probability,probabilityParameterisation_, 0);
		if (!normalized) {
			this.rate = rawRate_;
		} else {
			renormalize();
		}
	}
	
	/**
	 * This constructor creates a rate distribution with n rate classes, initially
	 * all equally probable and with rate = 1.
	 */
	public ArbitraryRates(int n, boolean normalized) {
		super(n);
		normalizedRates_ = normalized;
		probabilityParameterisation_ = new double[n-1];
		for (int i=0; i<n; i++) {
			probability[i] = 1/(double)n;
			rate[i] = 1;
		}
		StochasticVector.makeParameters(probability,probabilityParameterisation_, 0);
		if (!normalized) {
			rawRate_ = rate;
		} else {
			rawRate_ = new double[n];
			for (int i=0; i<n; i++) { rawRate_[i] = 1; }
		}
	}
	

	private void renormalize() {
		if (!normalizedRates_) return;
		double sum = 0;
		for (int i=0; i<numRates; i++) {
			sum += rawRate_[i]*probability[i];
		}
		for (int i=0; i<numRates; i++) {
			rate[i] = rawRate_[i]/sum;
		}
	}
	

	public double getDefaultValue(int n) {
		if (n<numRates-1) {
			return 0.5;
		} else {
			return 1;
		}
	}

	public double getLowerLimit(int n) { return 0; }
	public int getNumParameters()      { return numRates*2-1; }


	public double getParameter(int n) {
		if (n<numRates-1) {
			return probabilityParameterisation_[n];
		} else {
			return rate[n-numRates-1];
		}
	}

	public double getUpperLimit(int n) {
		if (n<numRates-1) {
			return 1;
		} else {
			return MAX_RATE;
		}
	}

	public void setParameter(double param, int n) {
		if (n<numRates-1) {
			probabilityParameterisation_[n] = param;
			StochasticVector.makeVector(probabilityParameterisation_, probability);
			renormalize();
		} else {
			rawRate_[n-numRates-1] = param;
			renormalize();
		}
	}

	public void setParameterSE(double paramSE, int n) {
		return;
	}

	public void report(PrintWriter out) {
		out.printf("Model of rate heterogeneity: Arbitrary, normalized = %b, probability vector = %s, rate vector = %s\n",
				normalizedRates_, pal.misc.Utils.toString(probability),pal.misc.Utils.toString(rate));
	}
}
