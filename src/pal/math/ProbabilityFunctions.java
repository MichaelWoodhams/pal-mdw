package pal.math;

/**
 * A class to collect any probability density function (PDF) or cumulative probability function (CDF) evaluations.
 * As I'm only adding to it as I find I need something, the selection will be sparse initially.
 * 
 * @author Michael Woodhams, School of IT, University of Sydney
 * @date 2009-04-29
 */

public class ProbabilityFunctions {
	
	public static double chiSquaredCDF(double x, int n) {
		return GammaFunction.incompleteGammaP(0.5*n,0.5*x);
	}
	
}
