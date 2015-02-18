package pal.math;
import pal.math.*;

/**
 * Given a sample drawn from a multinomial distribution, generate new samples drawn from the 'same' distribution.
 * (I.e. we assume the initial sample followed the distribution perfectly. Effectively we are bootstrapping
 * from the initial sample.)
 * 
 * @author woodhams
 *
 */

public class Multinomial {
	private int nCategories; // the number of categories
	private int cumulative[]; // the cumulative counts
	private int sampleSize;
	private MersenneTwisterFast rng;
	
	// TODO: check no negative frequencies
	public Multinomial(int[] frequencies) {
		setFrequencies(frequencies);
		rng = new MersenneTwisterFast();
	}
	public Multinomial(int[] frequencies, long seed) {
		this(frequencies);
		this.setSeed(seed);
	}
	public Multinomial(double[] frequencies) {
		setFrequencies(frequencies);
		rng = new MersenneTwisterFast();
	}
	public Multinomial(double[] frequencies, long seed) {
		this(frequencies);
		this.setSeed(seed);
	}
	
	public void setFrequencies(double[] frequencies) {
		int n = frequencies.length;
		int[] intFreq = new int[n];
		double sum=0;
		for (double f : frequencies) { sum+=f; }
		// scale the doubles to *large* integers, to reduce roundoff errors.
		// '-n' prevents any danger of round-off making the sum of intFreq overflow. 
		double scaleFactor = (Integer.MAX_VALUE-n)/sum;
		for (int i=0; i<n; i++) {
			intFreq[i] = (int)(frequencies[i]*scaleFactor);
		}
		setFrequencies(intFreq);
	}
	
	public void setFrequencies(int[] frequencies) {
		nCategories = frequencies.length;
		if (cumulative == null || cumulative.length != nCategories) {
			cumulative = new int[nCategories];
		}
		cumulative[0]=frequencies[0];
		for (int i=1; i<nCategories; i++) {
			cumulative[i] = frequencies[i]+cumulative[i-1];
		}
		sampleSize = cumulative[nCategories-1];
	}
	
	public void setSeed(long seed) {
		rng.setSeed(seed);
	}
	
	
	// TODO: versions of newSample where you pass the array the data is to be stored in. (More efficient)
	
	/**
	 * Draw a new sample of size n
	 * @return Vector of frequencies of the new sample
	 */
	public int[] newSample(int n) {
		int[] result = new int[nCategories];
		for (int i=0; i<n; i++) {
			result[draw()]++;
		}
		return result;
	}
	
	/**
	 * Draw a new sample of the same size as the original sample
	 * @return Vector of frequencies of the new sample
	 */
	public int[] newSample() {
		return newSample(sampleSize);
	}
	
	/**
	 * Draw a single item from the distribution
	 * @return the category number of the item drawn
	 */
	public int draw() {
		int random = rng.nextInt(sampleSize);  // range 0..sampleSize-1
		int i;
		for (i=0; random >= cumulative[i]; i++) {}
		return i;
	}
	
	// A cursory test. Will give false failure result three times in 2^100.
	public static boolean test() {
		// test 1: new Multinomial
		Multinomial multi = new Multinomial(new int[] {0,10,0,10,0});
		int[] result = multi.newSample(100);
		boolean firstTest = (result[0]==0 && result[1]>0 && result[2]==0 && result[3]>0 && result[4]==0);
		// test 2: change frequencies but not nCategories
		multi.setFrequencies(new int[] {10,0,0,0,10});
		result = multi.newSample(100);
		boolean secondTest = (result[0]>0 && result[1]==0 && result[2]==0 && result[3]==0 && result[4]>0);
		// test 3: change frequencies and nCategories
		multi.setFrequencies(new int[] {10,10,0});
		result = multi.newSample(100);
		boolean thirdTest = (result[0]>0 && result[1]>0 && result[2]==0);
		return firstTest && secondTest && thirdTest;
	}
}
