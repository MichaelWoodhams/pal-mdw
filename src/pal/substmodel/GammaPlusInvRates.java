package pal.substmodel;

import java.io.PrintWriter;

import pal.statistics.GammaDistribution;

/**
 * Discrete Gamma plus invariant sites distribution.
 * 
 * @author Michael Woodhams
 *
 */

// Created with the aid of cut-and-paste from GammaRates and InvariableSites
@SuppressWarnings("serial")
public class GammaPlusInvRates extends RateDistribution {
	private boolean showSE;

	// p[0] = alpha (shape parameter), p[1] = frac invariable.
	// se[] = standard errors.
	private double[] p;
	private double[] se;
	
	private static final double[] LOWER_BOUNDS = new double[]{0.001,0};
	private static final double[] UPPER_BOUNDS = new double[]{100,0.99};
	private static final double[] DEFAULT_VALS = new double[]{0.5,0.5};
	
	/**
	 * construct discrete Gamma distribution (mean = 1.0)
	 *
	 * @param n number of rate classes in the Gamma distribution
	 * @param a shape parameter (alpha)
	 * @param f fraction of invariable sites
	 */
	public GammaPlusInvRates(int n, double a, double f)
	{
		super(n+1);
		p = new double[]{a,f};
		se = new double[]{0,0};
		showSE = false;
		makeDistribution();
	}

	public String toString() {
		return "Discretized Gamma plus Invariant Sites distribution, fraction invariant="+Double.toString(p[1])+", " 
		+Integer.toString(numRates-1)+" Gamma rate classes, alpha = "+Double.toString(p[0]);
	}
	
	// interface Report

	public void report(PrintWriter out)
	{
		out.println("Model of rate heterogeneity: Discrete Gamma plus Invariant sites");
		out.printf("Fraction of invariant sites: %.3f",p[1]);
		if (showSE) {
			out.printf("  (S.E. %.3f)\n", se[1]);
		} else {
			out.println();
		}
		out.println("Number of rate categories: " + numRates);
		out.print("Gamma distribution parameter alpha: ");
		format.displayDecimal(out, p[0], 2);
		if (showSE)
		{
			out.print("  (S.E. ");
			format.displayDecimal(out, se[0], 2);
			out.println(")");
		}
		else
		{
			out.println();
		}
		out.println();
		printRates(out);
	}
	
	// interface Parameterized

	public int getNumParameters()
	{
		return 2;
	}

	public void setParameter(double param, int n)
	{
		p[n]=param;
		makeDistribution();
	}

	public double getParameter(int n)
	{
		return p[n];
	}

	public void setParameterSE(double paramSE, int n)
	{
		se[n] = paramSE;
		showSE = true;
	}

	public double getLowerLimit(int n)   { return LOWER_BOUNDS[n]; }
	public double getUpperLimit(int n)   { return UPPER_BOUNDS[n]; }
	public double getDefaultValue(int n) { return DEFAULT_VALS[n]; }


	//
	// Private stuff
	//


	private void makeDistribution()
	{
		double alpha = p[0]; // aliases used for ease of understanding 
		double invFrac = p[1];
		int numGammaRates = numRates-1;
		
		double mean = 0.0;
		double gammaClassWeight = (1-invFrac)/numGammaRates; 
		for (int i = 0; i < numGammaRates; i++)
		{
			rate[i] = GammaDistribution.quantile((2.0*i+1.0)/(2.0*numRates), alpha, 1.0/alpha);
			probability[i] = gammaClassWeight;
			mean += rate[i]*gammaClassWeight;
		}
		rate[numRates-1] = 0;
		probability[numRates-1] = invFrac;
		// Normalize so mean rate = 1:
		for (int i = 0; i < numRates; i++)
		{
			rate[i] /= mean;
		}
		fireParametersChangedEvent();
	}


}
