// InvariableSites.java
//
// (c) 2000--2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.substmodel;

import pal.misc.*;
import pal.math.*;

import java.io.*;


/**
 * invariable sites model (two-rate model with mean rate = 1.0)
 *
 * @version $Id: InvariableSites.java,v 1.4 2003/06/11 05:26:46 matt Exp $
 *
 * @author Korbinian Strimmer
 */
/*
 * MDW: Modified to allow upper limit to be changed.
 */
public class InvariableSites extends RateDistribution implements Serializable
{
	private double upperLimit_=1.0;
	//
	// Public stuff
	//

	/**
	 * construct discrete rate distribution with two rates
	 * (one invariable and one variable)
	 *
	 * @param f fraction of invariable sites
	 */
	public InvariableSites(double f)
	{
		super(2);
		frac = f;
		showSE = false;

		makeDistrib(frac);
	}

	public String toString() {
		return "Invariable sites distribution, "+Double.toString(frac)+" of sites invariable";
	}
	
	public void setUpperLimit(double limit) { upperLimit_=limit; };
	
	// interface Report

	public void report(PrintWriter out)
	{
		out.println("Model of rate heterogeneity: Invariable sites model");
		out.println("Number of rate categories: " + numRates);
		out.print("Fraction of invariable sites: ");
		format.displayDecimal(out, frac, 2);
		if (showSE)
		{
			out.print("  (S.E. ");
			format.displayDecimal(out, fracSE, 2);
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
		return 1;
	}

	public void setParameter(double param, int n)
	{
		frac = param;
		makeDistrib(frac);
	}

	public double getParameter(int n)
	{
		return frac;
	}

	public void setParameterSE(double paramSE, int n)
	{
		fracSE = paramSE;
		showSE = true;
	}

	public double getLowerLimit(int n)
	{
		return 0.0;
	}

	/**
	 * When evaluating likelihoods on the upper bound of 1.0, we'll get likelihood zero
	 * and hence non-numeric log likelihood, which causes some optimizers (notably
	 * ConjugateDirectionSearch) to fail. 
	 */
	public double getUpperLimit(int n)
	{
		return upperLimit_;
	}

	public double getDefaultValue(int n)
	{
		return 0.0;
	}


	//
	// Private stuff
	//

	private boolean showSE;

	// fraction of invariable sites
	private double frac;
	private double fracSE;

	private void makeDistrib(double f)
	{
		rate[0] = 0.0;
		rate[1] = 1.0/(1.0-f);  // ensures that mean rate = 1.0

		probability[0] = f;
		probability[1] = 1.0-f;
		fireParametersChangedEvent();
	}
	public boolean isDistributionIndependentlyMutable() { return false; }
}
