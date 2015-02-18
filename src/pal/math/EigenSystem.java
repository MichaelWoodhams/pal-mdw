package pal.math;

//import pal.jama.Matrix;
//import pal.jama.EigenvalueDecomposition;
import Jama.Matrix;  // overrides pal.math.Matrix
import Jama.EigenvalueDecomposition;

/**
 * Finds and stores e-vectors, e-values, inverse-evectors-matrix for a real
 * possibly non-symmetric matrix. Can then calculate matrix log, exponential etc
 * etc.
 * 
 * IMPORTANT!!!!
 * This class uses the Jama Matrix library. Class 'Matrix' is a Jama Matrix, not a PAL matrix.
 * For methods which return a matrix, I have provided two versions: an "M" version which returns
 * a Jama matrix, and an A version which returns a Java 2D array.
 * 
 * @author woodhams
 */
/*
 *
 * 
 * Typical usages:
 * We have a rate matrix and a branch length and desire a Markov matrix
 * 	 EigenSystem es = new EigenSystem(rateMatrix, true);
 *   Matrix markov = es.getMatrixPowerM(branchLength);
 * (can repeat with different branch lengths.)
 * 
 * We have a Markov matrix, and desire a rate matrix. (Not always possible, see
 * for example "Why PAM matrices can't be written as an exponential of a matrix"
 * by Carolin Kosiol and Nick Goldman.) 
 *   EigenSystem es = new EigenSystem(markovMatrix,false); // 'false' parameter is optional
 *   double[][] rateMatrix = es.getMatrixLogA();
 * This returns rate matrix Q such that initial Markov matrix P = exp(Q). 
 * 
 *  If you want to reuse an existing EigenSystem object for a new matrix, 
 *  use the setNew[Log]Matrix methods. 
 */

/* 
 * TODO: Routines like getMatrixPower{A|M} should have an option where instead of returning a
 * matrix/array, you pass them a matrix/array for the data to be stored in. This can avoid
 * unnecessary object creation and deletion. 
 * 
 * TODO: deal with case when log of a matrix is not real. Currently throws a RuntimeException.
 * Perhaps change to throw a non-runtime exception, and then deal with the fallout
 * in routines which call it.
 */

// TODO: getMatrixExp function - optionally take a constant multiplier as an argument
// (i.e. return exp(k*A)). The routine will look very much like the getMatrixLog routine.
// In the mean time, you can create an EigenSystem by specifying the log of the matrix,
// and then ask for the original matrix to power 1 (or to power k.)
// TODO: static 'matrixExp' and 'matrixLog' methods which create an EigenSystem and then
// return the required matrix, and then forget the EigenSystem object.

public class EigenSystem {
	// eigenvectors, eigenvalues contained 'eigen'. In the case of complex
	// eigenvalues, the 'eigenvector' matrix is
	// not really the eigenvectors: it satisfies AV = VD (A = original matrix, V
	// = 'eigenvector' matrix, D = diagonal
	// eigenvalue matrix.) If all eigenvalues are real, D is diagonal and V is
	// the eigenvector matrix. If eigenvalues
	// are complex (in which case they come in conjugate pairs) D is
	// block-diagonal with blocks of the form
	// {{x, y},{-y,x}} where x is the real part of the eigenvalue and y the
	// imaginary part. The two 'eigenvectors'
	// (corresponding columns in V) are real vectors which can be linearly
	// combined (with complex weights)
	// to create the actual complex eigenvectors (which are also a conjugate
	// pair.)
	private int n;              // matrix size, cached for convenience.
	// Following three taken from EigenvalueDecomposition. I can't just store the EigenvalueDecomposition because we need
	// access to these for scaling (rescaleByPower, normaliseToDeterminant).
	private Matrix eVec;        // Matrix of eigenvectors (except see above for what happens if there are complex eigenvalues.)
	private double[] reEVal;    // Real part of eigenvalues
	private double[] imEVal;    // Imaginary part of eigenvalues
	private double[] reLogEVal; // Real part of log of eigenvalues. Will be null
								// if there are real non-positive e-values.
	private double[] imLogEVal; // Imaginary part. Will null if not required
	private Matrix invEVec;     // inverse of the eigenvectors matrix
	
	// TODO: methods to allow reusing object with a new matrix (of same size) to avoid the need for much object creation and destruction

	/**
	 * Null constructor. Call setNewMatrix or setNewLogMatrix after using this.
	 */
	public EigenSystem() {		
	}
	
	
	/**
	 * Constructor taking a Jama library Matrix as input
	 * If 'fromLogMatrix' is true, the supplied matrix is the (matrix) logarithm of our
	 * target matrix, otherwise it is the target matrix itself.
	 * For example, an instantaneous rate matrix is the log of a Markov matrix.
	 */
	public EigenSystem (Matrix a, boolean fromLogMatrix) {
		this(); 
		if (fromLogMatrix) {
			this.setNewLogMatrix(a);
		} else {
			this.setNewMatrix(a);
		}
	}
	
	/**
	 * Constructor taking a Jama library Matrix as input
	 */
	public EigenSystem(Matrix a) {
		this();
		this.setNewMatrix(a);
	}
	
	/**
	 * Constructor taking array as input
	 */
	public EigenSystem(double[][] a) {
		this(new Matrix(a));
	}
	public EigenSystem(double[][] a, boolean fromLogMatrix) {
		this(new Matrix(a), fromLogMatrix);
	}

	
	
	/**
	 * Set the matrix to be processed. 
	 */
	public void setNewMatrix(double[][] a) {
		setNewMatrix(new Matrix(a));
	}
	public void setNewMatrix(Matrix a) {
		checkMatrix(a);
		n = a.getColumnDimension();
		EigenvalueDecomposition eigen = a.eig();
		eVec = eigen.getV();
		reEVal = eigen.getRealEigenvalues();
		imEVal = eigen.getImagEigenvalues();
		boolean complexLogEVals = false;
		boolean nonpositiveEVals = false;
		reLogEVal = new double[n];
		imLogEVal = new double[n];
		for (int i = 0; i < n; i++) {
			if (reEVal[i] > 0 && imEVal[i] == 0) {
				reLogEVal[i] = Math.log(reEVal[i]);
			} else {
				if (imEVal[i] != 0) {
					complexLogEVals = true;
					reLogEVal[i] = Complex.realPartComplexLog(reEVal[i],
							imEVal[i]);
					imLogEVal[i] = Complex.imagPartComplexLog(reEVal[i],
							imEVal[i]);
				} else {
					// real, non-positive eigenvalue. I haven't figured out
					// how/if this can be handled.
					nonpositiveEVals = true;
				}
			}
		}
		if (nonpositiveEVals) {
			reLogEVal = null;
			imLogEVal = null;
		}
		if (!complexLogEVals) {
			imLogEVal = null;
		}
		// TODO: try to handle the case where V is poorly conditioned or
		// singular.
		invEVec = eigen.getV().inverse();
	}
	
	/**
	 * Counterpart to setNewMatrix. We are given the log matrix, so get the log eigenvalues directly
	 * from the eigendecomposition, and then have to find the 'proper' eigenvalues.
	 */
	public void setNewLogMatrix(double[][] a) {
		setNewLogMatrix(new Matrix(a));
	}
	public void setNewLogMatrix(Matrix a) {
		checkMatrix(a);
		n = a.getColumnDimension();
		EigenvalueDecomposition eigen = a.eig();
		eVec = eigen.getV();
		reLogEVal = eigen.getRealEigenvalues();
		imLogEVal = eigen.getImagEigenvalues();
		reEVal = new double[n];
		imEVal = new double[n];
		boolean complexLogEVals = false;
		for (int i = 0; i < n; i++) {
			if (imLogEVal[i] == 0) {
				reEVal[i] = Math.exp(reLogEVal[i]);
				imEVal[i] = 0;
			} else {
				complexLogEVals = true;
				reEVal[i] = Complex.realPartComplexExp(reLogEVal[i], imLogEVal[i]);
				imEVal[i] = Complex.imagPartComplexLog(reLogEVal[i], imLogEVal[i]);
			}
		}
		if (!complexLogEVals) {
			imLogEVal = null;
		}
		// TODO: try to handle the case where V is poorly conditioned or
		// singular.
		invEVec = eigen.getV().inverse();
	}
	
	// Throw error if matrix has any non-finite entries
	private void checkMatrix(Matrix a) {
		for (int i=0; i<a.getRowDimension(); i++) {
			for (int j=0; j<a.getColumnDimension(); j++) {
				if (Double.isInfinite(a.get(i, j)) || Double.isNaN(a.get(i,j))) {
					throw new IllegalArgumentException("Non-finite entry in EigenSystem input matrix");
				}
			}
		}
	}

	// 'get' methods which return a Jama Matrix:
	/**
	 * @return Jama Matrix of eigenvectors as columns.
	 * NOTE: If there are complex eigenvalues, special terms apply.
	 * See the Jama EigenDecomposition class for details
	 */
	public Matrix getEigenVectorsM() {
		return eVec;
	}
	/**
	 * @return 2D array containing eigenvectors as columns.
	 * NOTE: If there are complex eigenvalues, special terms apply.
	 * See the Jama EigenDecomposition class for details
	 */
	public double[][] getEigenVectorsA() {
		return eVec.getArray();
	}
	/**
	 * @return Inverse of the eigenvector matrix, as a Jama Matrix.
	 * See getEigenVectorsM() for warnings about complex eigenvalues.
	 */
	public Matrix getInverseEigenMatrixM() {
		return invEVec;
	}
	/**
	 * @return Inverse of the eigenvector matrix, as a 2D array.
	 * See getEigenVectorsM() for warnings about complex eigenvalues.
	 */
	public double[][] getInverseEigenMatrixA() {
		return invEVec.getArray();
	}

	/**
	 * @return The real parts of the eigenvalues. 
	 * Order corresponds to the order of eigenvectors returned by getEigenVectors{A|M}.
	 */
	public double[] getRealEigenvalues() {
		return reEVal;
	}
	/**
	 * @return The imaginary parts of the eigenvalues. 
	 * Order corresponds to the order of eigenvectors returned by getEigenVectors{A|M}.
	 */
	public double[] getImagEigenvalues() {
		return imEVal;
	}

	/**
	 * Return the matrix logarithm of the specified matrix.
	 * An important application of this is that if the specified matrix was a Markov matrix,
	 * the matrix logarithm supplies the instantaneous rate matrix which would generate
	 * that Markov matrix over a time of 1. (The rate matrix can be linearly scaled
	 * to change the time scale.)
	 * @return Matrix logarithm as a 2D array.
	 */
	public double[][] getMatrixLogA() {
		return getMatrixLogM().getArray();
	}
	/**
	 * Return the matrix logarithm of the specified matrix.
	 * An important application of this is that if the specified matrix was a Markov matrix,
	 * the matrix logarithm supplies the instantaneous rate matrix which would generate
	 * that Markov matrix over a time of 1. (The rate matrix can be linearly scaled
	 * to change the time scale.)
	 * @return Matrix logarithm as a Jama Matrix
	 */

	public Matrix getMatrixLogM() {
		if (reLogEVal == null) {
			// The simplest case - we can't do anything (or maybe we can and I
			// haven't figured out what)
			throw new RuntimeException("Log of this matrix is not real");
		}

		if (imLogEVal == null) {
			// The simple case - all eigenvalues were real and positive. We
			// return V*log(D)*V^-1
			// where D is the diagonal matrix of eigenvalues. Knowing that D is
			// diagonal we can
			// avoid the need for a full matrix multiplication of 3 matrices.
			double[][] logA = new double[n][n];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < n; k++) {
						logA[i][j] += eVec.get(i, k) * reLogEVal[k]
								* invEVec.get(k, j);
					}
				}
			}
			return new Matrix(logA);
		}
		// Else the complex (literally) case. We return V*logD*V^-1 where logD
		// is a block-diagonal matrix
		// Where eigenvalue lambda[i] real, positive, we get a diagonal
		// entry log(lambda[i]). Where
		// it is complex, the conjugate pair of eigenvalues get a block with
		// Re(log(lambda[i])) on the
		// diagonal and +/- Im(log(lambda[i])) on the other corners of the
		// block.
		Matrix logD = new Matrix(n, n);
		for (int i = 0; i < n; i++) {
			if (imLogEVal[i] == 0) {
				logD.set(i, i, reLogEVal[i]);
			} else {
				// Complex case. Deal with two entries at once.
				// This code relies on complex eigenvalues coming in
				// adjacent conjugate pairs.
				logD.set(i, i, reLogEVal[i]);
				logD.set(i, i + 1, imLogEVal[i]);
				i++;
				logD.set(i, i, reLogEVal[i]);
				logD.set(i, i - 1, imLogEVal[i]);
			}
		}
		// This could be more efficient if I took into account the tridiagonal
		// nature of logD.
		// As mostly this is used for 4x4 matrices, the penalty is fairly low,
		// so I sacrifice speed for program clarity.
		return eVec.times(logD.times(invEVec));
	}
	

	/**
	 * Return original matrix raised to a power. (NOT restricted to integer powers)
	 * @param power
	 * @return matrix^power as a 2D array
	 */
	public double[][] getMatrixPowerA(double power) {
		return getMatrixPowerM(power).getArray();
	}
	/**
	 * Puts original matrix raised to a power into supplied storage space. (NOT restricted to integer powers)
	 * @param power
	 * @param store
	 * NOTE: Doesn't actually avoid any new array creation, and adds more work. Implemented to ease 
	 * conversion from obsolete PAL's MatrixExponential class. Doing this properly (avoiding
	 * array creation) would require changing JAMA to also allow avoiding array creation.
	 */
	public void getMatrixPowerA(double power, double[][] store) {
		double[][] a = getMatrixPowerM(power).getArray();
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				store[i][j] = a[i][j];
			}
		}
	}
	/**
	 * Also for ease of conversion from obsolete PAL's MatrixExponential class. See comments on 
	 * getMatrixPowerA(double power, double[][] store).
	 * @param power
	 * @param store
	 */
	public void getMatrixPowerTransposeA(double power, double[][] store) {
		double[][] a = getMatrixPowerM(power).getArray();
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				store[i][j] = a[j][i];
			}
		}
	}

	/**
	 * Return original matrix raised to a power. (NOT restricted to integer powers)
	 * @param power
	 * @return matrix^power as a Jama Matrix
	 */
	public Matrix getMatrixPowerM(double power) {
		if (reLogEVal == null) {
			throw new RuntimeException("Exp of this matrix may not be real");
		}
		// The real case:
		if (imLogEVal == null) {
			double[][] expScaleA = new double[n][n];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < n; k++) {
						expScaleA[i][j] += eVec.get(i, k)
								* Math.exp(power * reLogEVal[k])
								* invEVec.get(k, j);
					}
				}
			}
			return new Matrix(expScaleA);
		}
		// Else the complex case.
		Matrix D = new Matrix(n, n); // D = exp(scale * logD)
		for (int i = 0; i < n; i++) {
			if (imLogEVal[i] == 0) {
				D.set(i, i, Math.exp(power * reLogEVal[i]));
			} else {
				// Complex case. Deal with two entries at once. Relies on
				double re = Complex.realPartComplexExp(power * reLogEVal[i],
						power * imLogEVal[i]);
				double im = Complex.imagPartComplexExp(power * reLogEVal[i],
						power * imLogEVal[i]);
				D.set(i, i, re);
				D.set(i, i + 1, im);
				i++;
				D.set(i, i, re);
				D.set(i, i - 1, -im);
			}
		}
		return eVec.times(D.times(invEVec));
	}

	/**
	 * Effectively replace original matrix A by A^power.
	 * I.e. from now on, getMatrixPower and getMatrixLog will return values as
	 * if the constructor had been called on (A^power) instead of on A.
	 */
	public void rescaleByPower(double power) {
		if (reLogEVal != null) {
			for (int i=0; i<n; i++) {
				reLogEVal[i] *= power;
			}
		}
		if (imLogEVal != null) {
			for (int i=0; i<n; i++) {
				imLogEVal[i] *= power;
			}
		}
		
		for (int i=0; i<n; i++) {
			if (imEVal[i] == 0) {
				reEVal[i] = Math.pow(reEVal[i], power);
			} else {
				reEVal[i] = Complex.realPartComplexExp(reLogEVal[i], imLogEVal[i]);
				imEVal[i] = Complex.imagPartComplexExp(reLogEVal[i], imLogEVal[i]);
			}
		}
	}
	
	/**
	 * Find the power to raise the matrix to to make it have the desired determinant,
	 * then raise it to that power. After this operation, getMatrixPower(1) will
	 * return a matrix with the desired determinant.
	 * Current determinant and new determinant must have the same sign. 	
	 * @param det
	 */
	public void normalizeToDeterminant(double det) {
		double currentDet = 1;
		for (int i=0; i<n; i++) {
			if (imEVal[i] ==0) {
				currentDet *= reEVal[i];
			} else {
				// take care of the conjugate pair of eigenvalues both at once:
				currentDet *= (reEVal[i]*reEVal[i]+imEVal[i]*imEVal[i]);
				i++;
			}
		}
		if (currentDet * det <= 0) {
			throw new RuntimeException("Determinant to normalize to has wrong sign (or original matrix is singular)");
		}
		double power = Math.log(det)/Math.log(currentDet);
		this.rescaleByPower(power);
	}

	/**
	 * For testing and debugging this class.
	 * A 4x4 matrix with complex eigenvalues
	 */
	public static double[][] testMatrix1 = new double[][] {
		{ 0.93, 0.01, 0.05, 0.01 }, { 0.04, 0.92, 0.01, 0.03 },
		{ 0.03, 0.02, 0.94, 0.01 }, { 0.01, 0.05, 0.01, 0.97 } };
	/**
	 * For testing and debugging this class.
	 * A 4x4 matrix with all real eigenvalues
	 */
	public static double[][] testMatrix2 = new double[][] {
		{ 0.70, 0.20, 0.05, 0.05 }, { 0.10, 0.60, 0.10, 0.20 },
		{ 0.05, 0.05, 0.80, 0.10 }, { 0.10, 0.05, 0.05, 0.80 } };
	
	/**
	 * Run 'test' on 'testMatrix1'.
	 */
	public static void test() {
		test(testMatrix1);
	}

	/**
	 * Runs tests on the supplied array. Prints results to stdout.
	 * @param array
	 */
	public static void test(double[][] array) {
		Matrix a = new Matrix(array);
		// EigenvalueDecomposition ed = new EigenvalueDecomposition(a);
		EigenSystem es = new EigenSystem(a);
		int n = array.length;

		System.out.print("=========== Test of EigenSystem class ============\nTest matrix:");
		a.print(18, 13);
		System.out.printf("Determinant = %f\n", a.det());

		Matrix eVec = es.getEigenVectorsM();
		System.out.print("Eigenvectors: ");
		eVec.print(18, 13);

		System.out.print("Eigenvalues:\n");
		for (int i = 0; i < n; i++) {
			System.out.printf("%12f + %12fi\n", es.getRealEigenvalues()[i], es.getImagEigenvalues()[i]);
		}

		System.out.print("Matrix log:");
		Matrix logA = es.getMatrixLogM();
		logA.print(18, 13);

		System.out.print("Recovering the original matrix by matrix exp with scale of 1:");
		Matrix test = es.getMatrixPowerM(1);
		test.print(18, 13);
		System.out.printf("Max deviation from match = %e\n", maxdiff(a, test));

		Matrix sqrtA = es.getMatrixPowerM(0.5);
		System.out.print("Recovering the original matrix by squaring the square root:");
		test = sqrtA.times(sqrtA);
		test.print(18, 13);
		System.out.printf("Max deviation from match = %e\n", maxdiff(a, test));
		
		System.out.print("Recovering orignal matrix by scaling to power 1/3, then cubing: ");
		es.rescaleByPower(1./3.);
		test = es.getMatrixPowerM(3);
		System.out.printf("Max deviation from match = %e\n", maxdiff(a, test));
		
		System.out.print("Testing creating eigensystem from given log matrix. Should recover original matrix.");
		EigenSystem esfl = new EigenSystem(logA,true);
		test = esfl.getMatrixPowerM(1);
		test.print(18, 13);
		System.out.printf("Max deviation from match = %e\n", maxdiff(a, test));
		
		System.out.print("Scaling to determinant of 0.01");
		es.normalizeToDeterminant(0.01);
		test = es.getMatrixPowerM(1);
		test.print(18, 13);
		System.out.printf("Determinant = %g\n",test.det());
	}

	// For use by the test routine
	private static double maxdiff(Matrix a, Matrix b) {
		double maxerr = 0;
		int n = a.getColumnDimension();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				maxerr = Math.max(maxerr, a.get(i, j) - b.get(i, j));
			}
		}
		return maxerr;
	}
}