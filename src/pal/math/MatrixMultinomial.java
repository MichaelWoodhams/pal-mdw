package pal.math;
//import pal.math.*;
import Jama.Matrix;

/**
 * Produces random multinomial distribution, where the class of objects being chosen by the multinomial is
 * the elements of an array. E.g. given a 'true' divergence matrix, it can return a sampled divergence matrix.
 * Internally, we store frequencies as integers. If given frequencies as reals, we multiply by a big int
 * and round to get an approximation.
 * 
 * The real work is done by the Multinomial class.
 * 
 * @author woodhams
 *
 */

public class MatrixMultinomial {
	private Multinomial flatMultinomial = null;
	private int n, m; // matrix rows/columns
	private int defaultSampleSize = 0;
	
	
	/**
	 * @param freq - integer array of frequencies. Must be non-negative (not checked for.)
	 */
	public MatrixMultinomial(int[][] freq) {
		this.setFrequencies(freq);
	}
	/**
	 * Constructor with random number seed
	 * @param freq
	 * @param seed
	 */
	public MatrixMultinomial(int[][] freq, long seed) {
		this.setFrequencies(freq);
		this.setSeed(seed);
	}
	
	
	/**
	 * Constructor for real valued frequencies (i.e. non-integer.)
	 * Frequencies must be non-negative (not checked.)
	 * @param freq
	 */
	public MatrixMultinomial(double[][] freq) {
		this.setFrequencies(freq);
	}
	public MatrixMultinomial(double[][] freq, long seed) {
		this(freq);
		this.setSeed(seed);
	}

	/**
	 * Constructor taking a Jama library Matrix
	 * @param freq
	 */
	public MatrixMultinomial(Matrix freq) {
		this(freq.getArray());
	}
	public MatrixMultinomial(Matrix freq, long seed) {
		this(freq);
		this.setSeed(seed);
	}

	
	public void setFrequencies(int[][] freq) {
		n = freq.length;
		m = freq[0].length;
		int[] flattenedFreq = new int[n*m];
		for (int i=0; i<n; i++) {
			for (int j=0; j<m; j++) {
				flattenedFreq[i*m+j]=freq[i][j];
				defaultSampleSize += freq[i][j];
			}
		}
		if (flatMultinomial == null) {
			flatMultinomial = new Multinomial(flattenedFreq);
		} else {
			flatMultinomial.setFrequencies(flattenedFreq);
		}
	}
	public void setFrequencies(double[][] freq) {
		n = freq.length;
		m = freq[0].length;
		double sum=0;
		for (int i=0; i<n; i++) {
			for (int j=0; j<m; j++) {
				sum+=freq[i][j];
			}
		}
		// scale frequencies to range nearly up to max int value. (Leave buffer against rounding.)
		double mult = (Integer.MAX_VALUE - n*m)/sum; 
		int[] flattenedFreq = new int[n*m];
		for (int i=0; i<n; i++) {
			for (int j=0; j<m; j++) {
				flattenedFreq[i*m+j]=(int)(freq[i][j]*mult);
			}
		}
		if (flatMultinomial == null) {
			flatMultinomial = new Multinomial(flattenedFreq);
		} else {
			flatMultinomial.setFrequencies(flattenedFreq);
		}
	}
	
	public void setSeed(long seed) {
		flatMultinomial.setSeed(seed);
	}

	/**
	 * If you used an integer matrix constructor, this method will return an array with the same
	 * sum (i.e. number of samples taken) as the one supplied to the constructor. 
	 * @return
	 */
	public int[][] newSample() {
		return newSample(defaultSampleSize);
	}
	
	/**
	 * Return array of integers which sum to nSamples, sampled from the given distribution.
	 * @param nSamples
	 * @return
	 */
	public int[][] newSample(int nSamples) {
		int[][] result = new int[n][m];
		for (int i=0; i<nSamples; i++) {
			int flat = flatMultinomial.draw(); // number from 0 to n*m-1
			int x = flat / m;
			int y = flat % m;
			result[x][y]++;
		}
		return result;
	}
	
	/**
	 * Draw a single time from the distribution. Returns resulting {x,y} in store.
	 * @param store
	 */
	public void draw(int[] store) {
		int flat = flatMultinomial.draw(); // number from 0 to n*m-1
		store[0] = flat / m;
		store[1] = flat % m;
	}
	
	// A cursory test. Will give false failure result very very rarely
	public static boolean integerTest() {
		MatrixMultinomial multi = new MatrixMultinomial(new int[][] {{0,10,0},{10,10,0}});
		int[][] result = multi.newSample(300);
		return (result[0][0]==0 && result[0][1]>0 && result[0][2]==0 && result[1][0]>0 && result[1][1]>0 && result[1][2]==0);
	}
	public static boolean doubleTest() {
		MatrixMultinomial multi = new MatrixMultinomial(new double[][] {{0,0.1,0},{0.1,0.1,0}});
		int[][] result = multi.newSample(300);
		return (result[0][0]==0 && result[0][1]>0 && result[0][2]==0 && result[1][0]>0 && result[1][1]>0 && result[1][2]==0);
	}
	public static boolean test() {
		return integerTest() && doubleTest();
	}

	
}
