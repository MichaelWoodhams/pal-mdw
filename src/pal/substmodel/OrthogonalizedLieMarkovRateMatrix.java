package pal.substmodel;

import java.io.PrintStream;
//import java.util.concurrent.Callable;
import static java.lang.Math.abs;

import pal.datatype.DataType;
import pal.datatype.Nucleotides;

public abstract class OrthogonalizedLieMarkovRateMatrix implements NeoRateMatrix {
	private static final long serialVersionUID = -2905869052890989373L;
	
	private static final int A = Nucleotides.A_STATE;
	private static final int C = Nucleotides.C_STATE;
	private static final int G = Nucleotides.G_STATE;
	private static final int T = Nucleotides.UT_STATE;
	
	private enum Permutation {
		RY (new int[]{A,G,C,T}),
		WS (new int[]{A,T,C,G}),
		MK (new int[]{A,C,G,T});
		
		protected final int[] permutation;
		Permutation(int[] perm) {
			this.permutation = perm;
		}
		public static Permutation fromString(String perm) {
			if (perm == null) return RY; // default to transitions vs transversions
			switch (perm) {
			case "RY":
				return RY;
			case "WS":
				return WS;
			case "MK":
				return MK;
			default:
				throw new IllegalArgumentException("Unrecognized permutation name "+perm);
			}
		}
		
	}

	private String name_;
	private int nParam_;
	private SymmetryMatrix[] matrices_; // nParam_ of them.
	private double[] nativeParameters_; // (a,b,c etc.) the multipliers for the symmetry matrices. Preallocated storage space. Length = nParam_.
	private Permutation permutation_; // RY, WS or MK
	
	// Only for calling by subclass constructors.
	protected OrthogonalizedLieMarkovRateMatrix(String permName) {
		permutation_ = Permutation.fromString(permName);
	}
	
	/**
	 * Inner class SymmetryMatrix:
	 * 
	 * 
	 * @author woodhams
	 *
	 */
	private static class SymmetryMatrix {
		private int[][] positive; // indicies of cells to add to. n x 2.
		private int[][] negative; // indices of cells to subtract from. Will subtract 1 from each cell unless subtract3 is true
		private boolean subtract3; // only used for smB
		
		private SymmetryMatrix(int[][] pos, int[][] neg) {
			this(pos,neg,false);
		}
		
		private SymmetryMatrix(int[][] pos, int[][] neg, boolean sub3) {
			positive = pos;
			negative = neg;
			subtract3 = sub3;
		}
		/*
		 * IMPORTANT NOTE!
		 * All the 'pictures' below are in columns-sum-to-zero order.
		 * pal works in rows-sum-to-zero, so these will get transposed in process of 
		 * building the rate matrix.
		 */
		
		protected void applyToRates(double[][] rateStore, double weight, Permutation perm) {
			// Note transposing coordinates in building rateStore 
			for (int[] coord : positive) {
				rateStore[perm.permutation[coord[1]]][perm.permutation[coord[0]]] += weight;
			}
			if (subtract3) weight *= 3;
			for (int[] coord : negative) {
				rateStore[perm.permutation[coord[1]]][perm.permutation[coord[0]]] -= weight;
			}
		}
		protected void applyToStrings(StringBuffer[][] strings, char variable) {
			// Leave in column-sum order.
			for (int[] coord : positive) {
				strings[coord[0]][coord[1]].append('+').append(variable);
			}
			String factor = subtract3 ? "3" : "";
			for (int[] coord : negative) {
				strings[coord[0]][coord[1]].append('-').append(factor).append(variable);
			}
		}

		/* 
		 * Weight matrix for variable 'b' in all models.
		 * =+++
		 * +=++
		 * ++=+
		 * +++=
		 * where '=' means -3.
		 */
		public static final SymmetryMatrix smB = new SymmetryMatrix(
				new int[][]{{0,1},{0,2},{0,3},{1,0},{1,2},{1,3},{2,0},{2,1},{2,3},{3,0},{3,1},{3,2}},
				new int[][]{{0,0},{1,1},{2,2},{3,3}}, true);
		/*
		 * Weight matrix for variable 'a' in all models.
		 * -+00
		 * +-00
		 * 00-+
		 * 00+-
		 */
		public static final SymmetryMatrix smA = new SymmetryMatrix(
				new int[][]{{0,1},{1,0},{2,3},{3,2}},
				new int[][]{{0,0},{1,1},{2,2},{3,3}});
		/*
		 * 'Twisted' version of smA
		 * -+00
		 * +-00
		 * 00+-
		 * 00-+
		 */
		public static final SymmetryMatrix smAt = new SymmetryMatrix(
				new int[][]{{0,1},{1,0},{2,2},{3,3}},
				new int[][]{{0,0},{1,1},{2,3},{3,2}});
		/*
		 * 'Reversed' version of smA
		 * 00+-
		 * 00-+
		 * +-00
		 * -+00
		 */
		public static final SymmetryMatrix smAr = new SymmetryMatrix(
				new int[][]{{0,2},{1,3},{2,0},{3,1}},
				new int[][]{{0,3},{1,2},{2,1},{3,0}});
		/*
		 * Reversed and twisted smA
		 * 00+-
		 * 00-+
		 * -+00
		 * +-00
		 */
		public static final SymmetryMatrix smArt = new SymmetryMatrix(
				new int[][]{{0,2},{1,3},{2,1},{3,0}},
				new int[][]{{0,3},{1,2},{2,0},{3,1}});
		/*
		 * Horizontal matrix
		 * ++++
		 * ++++
		 * ----
		 * ----
		 */
		public static final SymmetryMatrix smH = new SymmetryMatrix(
				new int[][]{{0,0},{0,1},{0,2},{0,3},{1,0},{1,1},{1,2},{1,3}},
				new int[][]{{2,0},{2,1},{2,2},{2,3},{3,0},{3,1},{3,2},{3,3}});
		/*
		 * 'Up' matrix
		 * ++++
		 * ----
		 * 0000
		 * 0000
		 */
		public static final SymmetryMatrix smU = new SymmetryMatrix(
				new int[][]{{0,0},{0,1},{0,2},{0,3}},
				new int[][]{{1,0},{1,1},{1,2},{1,3}});
		/*
		 * Twisted up matrix
		 * ++--
		 * --++
		 * 0000
		 * 0000
		 */
		public static final SymmetryMatrix smUt = new SymmetryMatrix(
				new int[][]{{0,0},{0,1},{1,2},{1,3}},
				new int[][]{{1,0},{1,1},{0,2},{0,3}});
		/*
		 * Negative Twisted up matrix
		 * --++
		 * ++--
		 * 0000
		 * 0000
		 */
		public static final SymmetryMatrix smUtm = new SymmetryMatrix(
				new int[][]{{1,0},{1,1},{0,2},{0,3}},
				new int[][]{{0,0},{0,1},{1,2},{1,3}});
		/*
		 * 'down' matrix (avoid 'lower' to leave L free for 'left')
		 * 0000
		 * 0000
		 * ++++
		 * ----
		 */
		
		public static final SymmetryMatrix smD = new SymmetryMatrix(
				new int[][]{{2,0},{2,1},{2,2},{2,3}},
				new int[][]{{3,0},{3,1},{3,2},{3,3}});
		/*
		 * Twisted down matrix
		 * 0000
		 * 0000
		 * ++--
		 * --++
		 */
		public static final SymmetryMatrix smDt = new SymmetryMatrix(
				new int[][]{{2,0},{2,1},{3,2},{3,3}},
				new int[][]{{3,0},{3,1},{2,2},{2,3}});
		/*
		 * Twisted left matrix
		 * -+00
		 * -+00
		 * +-00
		 * +-00
		 */
		public static final SymmetryMatrix smLt = new SymmetryMatrix(
				new int[][]{{0,1},{1,1},{2,0},{3,0}},
				new int[][]{{0,0},{1,0},{2,1},{3,1}});
		/*
		 * Twisted right matrix
		 * 00+-
		 * 00+-
		 * 00-+
		 * 00-+
		 */
		public static final SymmetryMatrix smRt = new SymmetryMatrix(
				new int[][]{{0,2},{1,2},{2,3},{3,3}},
				new int[][]{{0,3},{1,3},{2,2},{3,2}});
	}

	/**
	 * Implicitly performs a normalization, so for 8 dimensional model, we'll have 7 parameters.
	 * (The normalization is a^2+b^2=1.)
	 * All parameters are in the range -1 to 1. 
	 * 
	 * Writes into nativeParameters_.
	 */
	protected abstract void translate(double[] parameters, int startIndex);
	

	/*
	 * 'min' and 'max' functions for various numbers of parameters.
	 * Some of this *should* have been avoidable by 'import static java.lang.Math.min',
	 * but either Java or Eclipse insisted on seeing min(a,b) as an error because min(a,b,c) existed.
	 */
	private static double min(double a, double b)           { return Math.min(a,b); }
	private static double max(double a, double b)           { return Math.max(a,b); }
	private static double min(double a, double b, double c) { return min(min(a, b),c);}
	private static double max(double a, double b, double c) { return max(max(a, b),c);}
	private static double min(double a, double b, double c, double d, double e) {
		return min(min(a,b,c),min(d,e));
	}
	
	/*
	 * NeoRateMatrix methods
	 */
	public void createRelativeRates(double[][] rateStore,
			double[] rateParameters, int startIndex) {
		// set nativeParameters
		translate(rateParameters, startIndex);
		// Zero rateStore
		int dim = getDimension();
		for (int i=0; i<dim; i++) {
			for (int j=0; j<dim; j++) {
				rateStore[i][j] = 0;
			}
		}
		int n = nativeParameters_.length;
		// Add weighted symmetry matrices
		for (int i=0; i<n; i++) {
			matrices_[i].applyToRates(rateStore, nativeParameters_[i], permutation_);
		}
		// Check for illegal values (should never happen) and 
		// remove negative-by-rounding-error values while we're at it.
		for (int row=0; row<4; row++) {
			for (int col=0; col<4; col++) {
				if (col != row && rateStore[row][col]<0) {
					if (rateStore[row][col]<-1e-12) 
						throw new RuntimeException("Non-zero off-diagonal");
					rateStore[row][col]=0;
				}
			}
		}
	}
	
	public String getUniqueName() { return name_; 	}
	public boolean isReversible() {	return false; 	}
	public boolean isIndependentOfEqbmFreq() { return false;	}
	public int getDimension() {		return 4;	}
	public int getNumberOfRateParameters() { return nParam_; }
	public boolean isDataTypeCompatible(DataType dt) { return dt.getNumStates()==4; }
	/*
	 * Technically, the allowable range for all parameters is -1 to 1. However, values on the
	 * boundaries lead to zeros in the rate matrix, which often leads to zero likelihood (hence
	 * infinite log likelihood) and ConjugateDirectionSearch in particular is both prone to
	 * evaluate on the boundary, and is unable to handle infinities. So I hack it a bit by
	 * setting the range a little smaller than -1 to 1.
	 * 
	 * (In particular, letting parameter 0 be -1 gives b=0 and transversions are forbidden.)
	 */
	public double getRateParameterLowerBound(int parameter) { return -0.99999; }
	public double getRateParameterUpperBound(int parameter) { return  0.99999; }
	public void getDefaultRateParameters(double[] parameterStore, int startIndex) {
		int n = getNumberOfRateParameters();
		for (int i=startIndex; i<startIndex+n; i++) parameterStore[i] = 0;
	}

	// Not yet coded output that processes parameters.
	public void printObject(PrintStream out, String prefixString, double[] parameters) {
		this.printObject(out, prefixString);
		out.printf("%s...unfinished,  sorry...%n",prefixString); 
	}

	// I've been lazy about superfluous '+' at start of expression, ',' and end of row.
	private static final char[] LITTLE_VARIABLE_NAMES = new char[]{'a','b','c','d','e','f','g','h','i','j'};
	public void printObject(PrintStream out, String prefixString) {
		out.printf("%sOrthogLieMarkRateMat %s:\n",prefixString, name_);
		StringBuffer[][] matrix = new StringBuffer[4][4];
		for (int i=0; i<4; i++) for (int j=0; j<4; j++) matrix[i][j] = new StringBuffer();
		for (int i=0; i<=nParam_; i++) {
			matrices_[i].applyToStrings(matrix, LITTLE_VARIABLE_NAMES[i]);
		}
		int maxlen=0;
		for (int i=0; i<4; i++) for (int j=0; j<4; j++) maxlen = Math.max(maxlen, matrix[i][j].length());
		String format = "%-"+Integer.toString(maxlen)+"s, ";
		for (int row=0; row<4; row++) {
			out.print(prefixString+'{');
			for (int col=0; col<4; col++) {
				out.printf(format,matrix[row][col].toString());
			}
			out.printf("}%n");
		}	
	}

	public void printObject(PrintStream out) {
		this.printObject(out,"");
	}

	public void setUpperBounds(double[] bounds) {
		throw new RuntimeException("Can't alter parameter bounds");
	}


	public void setUpperBounds(double bound) {
		throw new RuntimeException("Can't alter parameter bounds");
	}
	
	/*
	 * The different models need different 'translate' functions. I've achieved this through 
	 * making them inner classes. However, there are other possibilities, such as a
	 * Callable data member. To preserve the option to change implementation, I
	 * provide a set of functions to return an instance of each model
	 */
	
	/*
	 * For models 2.2b, 3.3a, 3.3b, 3.3c, 4.4a, 4.4b, 6.6, 8.8 the 'weighted rays' implementation
	 * is as good as this one. I include these models for consistency. (Also avoids the need
	 * to mess with the weighted rays implementations to make them play nice with ConjugateDirectionSearch.)
	 */
	// 2.2b: Also known as Kimura 2 parameter model
	public static OrthogonalizedLieMarkovRateMatrix getM2r2bRY() {return new M2r2b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM2r2bWS() {return new M2r2b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM2r2bMK() {return new M2r2b("MK");}
	private static class M2r2b extends OrthogonalizedLieMarkovRateMatrix {
		protected M2r2b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "2.2b";
			super.nParam_ = 1;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,SymmetryMatrix.smB};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		protected void translate(double[] parameters, int startIndex) {
			double B, a, b;
			
			B = parameters[startIndex  ]; double maxp=abs(B);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
		}
	}
	
	// 3.3a is Kimura 3 parameter model. As it is fully symmetric, there is no RY/WS/MK distinction.
	public static OrthogonalizedLieMarkovRateMatrix getM3r3a() {return new M3r3a("RY");}
	private static class M3r3a extends OrthogonalizedLieMarkovRateMatrix {
		protected M3r3a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "3.3a";
			super.nParam_ = 2;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,SymmetryMatrix.smB,SymmetryMatrix.smAr};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// same as 3.3b
		protected void translate(double[] parameters, int startIndex) {
			double B, C, a, b, c;
			
			B = parameters[startIndex  ]; double maxp=abs(B);
			C = parameters[startIndex+1]; maxp = max(maxp, C);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;	
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM3r3bRY() {return new M3r3b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM3r3bWS() {return new M3r3b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM3r3bMK() {return new M3r3b("MK");}
	private static class M3r3b extends OrthogonalizedLieMarkovRateMatrix {
		protected M3r3b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "3.3b";
			super.nParam_ = 2;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,SymmetryMatrix.smB,SymmetryMatrix.smArt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// same as 3.3a
		protected void translate(double[] parameters, int startIndex) {
			double B, C, a, b, c;
			
			B = parameters[startIndex  ]; double maxp=abs(B);
			C = parameters[startIndex+1]; maxp = max(maxp, C);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;	
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM3r3cRY() {return new M3r3c("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM3r3cWS() {return new M3r3c("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM3r3cMK() {return new M3r3c("MK");}
	private static class M3r3c extends OrthogonalizedLieMarkovRateMatrix {
		protected M3r3c(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "3.3c";
			super.nParam_ = 2;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,SymmetryMatrix.smB,SymmetryMatrix.smAt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// same as 3.3a
		protected void translate(double[] parameters, int startIndex) {
			double B, C, a, b, c;
			
			B = parameters[startIndex  ]; double maxp=abs(B);
			C = parameters[startIndex+1]; maxp = max(maxp, C);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*(a+b);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;	
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM3r4RY() {return new M3r4("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM3r4WS() {return new M3r4("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM3r4MK() {return new M3r4("MK");}
	private static class M3r4 extends OrthogonalizedLieMarkovRateMatrix {
		protected M3r4(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "3.4";
			super.nParam_ = 2;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,SymmetryMatrix.smB,SymmetryMatrix.smH};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		protected void translate(double[] parameters, int startIndex) {
			double B, C, a, b, c;
			
			B = parameters[startIndex  ]; double maxp=abs(B);
			C = parameters[startIndex+1]; maxp = max(maxp, C);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*min(b, a+b);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;	
		}
	}
	
	// The Felsenstein 81 model. As it is fully symmetric, there is no RY/WS/MK distinction.
	public static OrthogonalizedLieMarkovRateMatrix getM4r4a() {return new M4r4a("RY");}
	private static class M4r4a extends OrthogonalizedLieMarkovRateMatrix {
		public M4r4a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "4.4a";
			super.nParam_ = 3;
			// This is the only model not to have smA in the first place
			super.matrices_ = new SymmetryMatrix[]
					{SymmetryMatrix.smH,SymmetryMatrix.smB,SymmetryMatrix.smU,SymmetryMatrix.smD};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, a, b, c, d;
			
			B = parameters[startIndex  ]; double maxp= abs(B);
			C = parameters[startIndex+1]; maxp = max(maxp, C);
			D = parameters[startIndex+2]; maxp = max(maxp, D);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// Not having smA as the first SymmetryMatrix, 'a' and 'b' are derived differently.
			// Sum of off-diagonals is 12b, so for trace=-4 we have 12b=4, or b=1/3. 
			b=1/3.;
			a=B*b;
			c=C*(b+a);
			d=D*(b-a);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;	
			super.nativeParameters_[3]=d;	
		}
		/*
		 *  To ensure that no row of off-diagonals is equal to zero (which will give likelihood=0)
		 *  we don't let any of the parameters actually reach 1 or -1
		 */
		public double getRateParameterLowerBound(int parameter) { return -0.99999; }
		public double getRateParameterUpperBound(int parameter) { return  0.99999; }
	}

	public static OrthogonalizedLieMarkovRateMatrix getM4r4bRY() {return new M4r4b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM4r4bWS() {return new M4r4b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM4r4bMK() {return new M4r4b("MK");}
	private static class M4r4b extends OrthogonalizedLieMarkovRateMatrix {
		public M4r4b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "4.4b";
			super.nParam_ = 3;
			super.matrices_ = new SymmetryMatrix[]
					{SymmetryMatrix.smA,SymmetryMatrix.smB,SymmetryMatrix.smAt,SymmetryMatrix.smH};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, a, b, c, d;
			
			B = parameters[startIndex  ]; double maxp= abs(B);
			C = parameters[startIndex+1]; maxp = max(maxp, C);
			D = parameters[startIndex+2]; maxp = max(maxp, D);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			d=D*b;
			c=C*(a+b)-d;
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;	
			super.nativeParameters_[3]=d;	
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM4r5aRY() {return new M4r5a("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM4r5aWS() {return new M4r5a("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM4r5aMK() {return new M4r5a("MK");}
	private static class M4r5a extends OrthogonalizedLieMarkovRateMatrix {
		public M4r5a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "4.5a";
			super.nParam_ = 3;
			super.matrices_ = new SymmetryMatrix[]
					{SymmetryMatrix.smA,SymmetryMatrix.smB,SymmetryMatrix.smAr,SymmetryMatrix.smH};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// Same as 4.5b
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, a, b, c, d;
			
			B = parameters[startIndex  ]; double maxp= abs(B);
			C = parameters[startIndex+1]; maxp = max(maxp, C);
			D = parameters[startIndex+2]; maxp = max(maxp, D);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			d=D*min(b,a+b);
			c=C*(b-abs(d));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;	
			super.nativeParameters_[3]=d;	
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM4r5bRY() {return new M4r5b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM4r5bWS() {return new M4r5b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM4r5bMK() {return new M4r5b("MK");}
	private static class M4r5b extends OrthogonalizedLieMarkovRateMatrix {
		public M4r5b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "4.5b";
			super.nParam_ = 3;
			super.matrices_ = new SymmetryMatrix[]
					{SymmetryMatrix.smA,SymmetryMatrix.smB,SymmetryMatrix.smArt,SymmetryMatrix.smH};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// Same as 4.5a
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, a, b, c, d;
			
			B = parameters[startIndex  ]; double maxp= abs(B);
			C = parameters[startIndex+1]; maxp = max(maxp, C);
			D = parameters[startIndex+2]; maxp = max(maxp, D);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			d=D*min(b,a+b);
			c=C*(b-abs(d));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;	
			super.nativeParameters_[3]=d;	
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM5r6aRY() {return new M5r6a("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r6aWS() {return new M5r6a("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r6aMK() {return new M5r6a("MK");}
	private static class M5r6a extends OrthogonalizedLieMarkovRateMatrix {
		public M5r6a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.6a";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
					                               SymmetryMatrix.smArt,
					                               SymmetryMatrix.smAr,
					                               SymmetryMatrix.smAt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*(b-abs(c));
			e=E*(a+b);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM5r6bRY() {return new M5r6b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r6bWS() {return new M5r6b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r6bMK() {return new M5r6b("MK");}
	private static class M5r6b extends OrthogonalizedLieMarkovRateMatrix {
		public M5r6b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.6b";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
					                               SymmetryMatrix.smH,
					                               SymmetryMatrix.smU,
					                               SymmetryMatrix.smD};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*min(b,a+b);
			d=D*min(a+b+c, b+c);
			e=E*min(a+b-c, b-c);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM5r11aRY() {return new M5r11a("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r11aWS() {return new M5r11a("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r11aMK() {return new M5r11a("MK");}
	private static class M5r11a extends OrthogonalizedLieMarkovRateMatrix {
		public M5r11a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.11a";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
					                               SymmetryMatrix.smAt,
					                               SymmetryMatrix.smU,
					                               SymmetryMatrix.smD};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// Same as 5.11b, 5.11c
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*(a+b);
			d=D*min(a+b+c, b);
			e=E*min(a+b-c, b);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM5r11bRY() {return new M5r11b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r11bWS() {return new M5r11b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r11bMK() {return new M5r11b("MK");}
	private static class M5r11b extends OrthogonalizedLieMarkovRateMatrix {
		public M5r11b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.11b";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
					                               SymmetryMatrix.smAt,
					                               SymmetryMatrix.smUt,
					                               SymmetryMatrix.smDt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// same as 5.11a, 5.11c
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*(a+b);
			d=D*min(a+b+c, b);
			e=E*min(a+b-c, b);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM5r11cRY() {return new M5r11c("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r11cWS() {return new M5r11c("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r11cMK() {return new M5r11c("MK");}
	private static class M5r11c extends OrthogonalizedLieMarkovRateMatrix {
		public M5r11c(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.11c";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
					                               SymmetryMatrix.smAt,
					                               SymmetryMatrix.smLt,
					                               SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// same as 5.11a, 5.11c
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*(a+b);
			d=D*min(a+b+c, b);
			e=E*min(a+b-c, b);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM5r16RY() {return new M5r16("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r16WS() {return new M5r16("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r16MK() {return new M5r16("MK");}
	private static class M5r16 extends OrthogonalizedLieMarkovRateMatrix {
		public M5r16(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.16";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
					                               SymmetryMatrix.smH,
					                               SymmetryMatrix.smLt,
					                               SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*min(b,a+b);
			d=D*min(a+b+c, b-c);
			e=E*min(a+b-c, b+c);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM5r7aRY() {return new M5r7a("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r7aWS() {return new M5r7a("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r7aMK() {return new M5r7a("MK");}
	private static class M5r7a extends OrthogonalizedLieMarkovRateMatrix {
		public M5r7a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.7a";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
					                               SymmetryMatrix.smAr,
					                               SymmetryMatrix.smU,
					                               SymmetryMatrix.smD};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// Same as 5.7b, 5.7c
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(a+b, b-abs(c));
			e=E*min(a+b, b-abs(c));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM5r7bRY() {return new M5r7b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r7bWS() {return new M5r7b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r7bMK() {return new M5r7b("MK");}
	private static class M5r7b extends OrthogonalizedLieMarkovRateMatrix {
		public M5r7b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.7b";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
					                               SymmetryMatrix.smAr,
					                               SymmetryMatrix.smUt,
					                               SymmetryMatrix.smDt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// Same as 5.7a, 5.7c
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(a+b, b-abs(c));
			e=E*min(a+b, b-abs(c));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM5r7cRY() {return new M5r7c("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r7cWS() {return new M5r7c("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM5r7cMK() {return new M5r7c("MK");}
	private static class M5r7c extends OrthogonalizedLieMarkovRateMatrix {
		public M5r7c(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "5.7c";
			super.nParam_ = 4;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
                                                   SymmetryMatrix.smB,
                                                   SymmetryMatrix.smAr,
                                                   SymmetryMatrix.smLt,
                                                   SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}

		// Same as 5.7a, 5.7b
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, a, b, c, d, e;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(a+b, b-abs(c));
			e=E*min(a+b, b-abs(c));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
		}
	}
	
	// 6.6WS is the Strand Symmetric Model
	public static OrthogonalizedLieMarkovRateMatrix getM6r6RY() {return new M6r6("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r6WS() {return new M6r6("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r6MK() {return new M6r6("MK");}
	public static OrthogonalizedLieMarkovRateMatrix getStrandSymmetricModel() {return new M6r6("WS");}
	private static class M6r6 extends OrthogonalizedLieMarkovRateMatrix {
		public M6r6(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "6.6";
			super.nParam_ = 5;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, a, b, c, d, e, f;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*b;
			f=(1-F)/2*(-b+abs(c+d))+(1+F)/2*(b-abs(c-d));
			e=E*(a+b)-f;
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
		}
	}

	// 6.7a is a fully symmetric model, so there is no RY/WS/MK distinction.
	public static OrthogonalizedLieMarkovRateMatrix getM6r7a() {return new M6r7a("RY");}
	private static class M6r7a extends OrthogonalizedLieMarkovRateMatrix {
		public M6r7a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "6.7a";
			super.nParam_ = 5;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		// same as 6.7b
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, a, b, c, d, e, f;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(a+b, b-abs(c));
			e=E*min(a+b+d, b+d-abs(c));
			f=F*min(a+b-d, b-d-abs(c));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM6r7bRY() {return new M6r7b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r7bWS() {return new M6r7b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r7bMK() {return new M6r7b("MK");}
	private static class M6r7b extends OrthogonalizedLieMarkovRateMatrix {
		public M6r7b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "6.7b";
			super.nParam_ = 5;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		// same as 6.7a
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, a, b, c, d, e, f;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(a+b, b-abs(c));
			e=E*min(a+b+d, b+d-abs(c));
			f=F*min(a+b-d, b-d-abs(c));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM6r17aRY() {return new M6r17a("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r17aWS() {return new M6r17a("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r17aMK() {return new M6r17a("MK");}
	private static class M6r17a extends OrthogonalizedLieMarkovRateMatrix {
		public M6r17a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "6.17a";
			super.nParam_ = 5;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smLt,
							                       SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		// same as 6.17b
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, a, b, c, d, e, f;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(a+b, b-abs(c));
			e=E*min(a+b+d, b-d-abs(c));
			f=F*min(a+b-d, b+d-abs(c));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM6r17bRY() {return new M6r17b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r17bWS() {return new M6r17b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r17bMK() {return new M6r17b("MK");}
	private static class M6r17b extends OrthogonalizedLieMarkovRateMatrix {
		public M6r17b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "6.17b";
			super.nParam_ = 5;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smLt,
							                       SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		// same as 6.17a
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, a, b, c, d, e, f;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(a+b, b-abs(c));
			e=E*min(a+b+d, b-d-abs(c));
			f=F*min(a+b-d, b+d-abs(c));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM6r8aRY() {return new M6r8a("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r8aWS() {return new M6r8a("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r8aMK() {return new M6r8a("MK");}
	private static class M6r8a extends OrthogonalizedLieMarkovRateMatrix {
		public M6r8a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "6.8a";
			super.nParam_ = 5;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, a, b, c, d, e, f;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			d=D*b;
			c=C*(a+b)-d;
			e=E*min(b+d,a+b+c+d);
			f=F*min(b-d,a+b-c-d);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM6r8bRY() {return new M6r8b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r8bWS() {return new M6r8b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM6r8bMK() {return new M6r8b("MK");}
	private static class M6r8b extends OrthogonalizedLieMarkovRateMatrix {
		public M6r8b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "6.8b";
			super.nParam_ = 5;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smLt,
							                       SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, a, b, c, d, e, f;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			d=D*b;
			c=C*(a+b)-d;
			e=E*min(b-d,a+b+c+d);
			f=F*min(b+d,a+b-c-d);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM8r8RY() {return new M8r8("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r8WS() {return new M8r8("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r8MK() {return new M8r8("MK");}
	private static class M8r8 extends OrthogonalizedLieMarkovRateMatrix {
		public M8r8(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "8.8";
			super.nParam_ = 7;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD,
							                       SymmetryMatrix.smUtm,
							                       SymmetryMatrix.smDt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, a, b, c, d, e, f, g, h;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			d=D*b;
			c=C*(a+b)-d;
			e=E*(a/2+b+c/2+d);
			g=(1-G)/2*max(-b-d-e,-a-b-c-d+e)+(1+G)/2*min(b+d-e,a+b+c+d+e);
			f=F*(a/2+b-c/2-d);
			h=(1-H)/2*max(-b+d-f,-a-b+c+d+f)+(1+H)/2*min(b-d-f,a+b-c-d+f);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM8r10aRY() {return new M8r10a("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r10aWS() {return new M8r10a("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r10aMK() {return new M8r10a("MK");}
	private static class M8r10a extends OrthogonalizedLieMarkovRateMatrix {
		public M8r10a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "8.10a";
			super.nParam_ = 7;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, a, b, c, d, e, f, g, h;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			e=E*(a+2*b);
			f=F*(a+b)-e;
			f=(1-F)/2*max(-b,-a-b-e)+(1+F)/2*min(b,a+b-e);
			g=G*min(b+f,a+b+e+f);
			h=H*min(b-f,a+b-e-f);
			c=C*(b-abs(g)/2-abs(h)/2);
			d=(1-D)/2*max(-b+c+f+abs(h),-b-c-f+abs(g))
		     +(1+D)/2*min( b-c+f-abs(g), b+c-f-abs(h));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM8r10bRY() {return new M8r10b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r10bWS() {return new M8r10b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r10bMK() {return new M8r10b("MK");}
	private static class M8r10b extends OrthogonalizedLieMarkovRateMatrix {
		public M8r10b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "8.10b";
			super.nParam_ = 7;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smLt,
							                       SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, a, b, c, d, e, f, g, h;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			e=E*(a+2*b);
			f=(1-F)/2*max(-b,-a-b-e)+(1+F)/2*min(b,a+b-e);
			g=G*min(b-f,a+b+e+f);
			h=H*min(b+f,a+b-e-f);
			c=C*(b-abs(g)/2-abs(h)/2);
			d=(1-D)/2*max(-b+c+f+abs(g),-b-c-f+abs(h))
		     +(1+D)/2*min( b-c+f-abs(h), b+c-f-abs(g));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM8r16RY() {return new M8r16("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r16WS() {return new M8r16("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r16MK() {return new M8r16("MK");}
	private static class M8r16 extends OrthogonalizedLieMarkovRateMatrix {
		public M8r16(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "8.16";
			super.nParam_ = 7;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD,
							                       SymmetryMatrix.smLt,
							                       SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, a, b, c, d, e, f, g, h;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*(a+2*b);
			d=(1-D)/2*max(-a-b-c,-b)+(1+D)/2*min(a+b-c,b);
			e=E*min(b+d,a+2*b,a+2*b+c);
			f=F*min(b-d,a+2*b-abs(c)-abs(e));
			g=(1-G)/2*max(-a-b-c-d-e,-b+d+abs(f))
		     +(1+G)/2*min( a+b+c+d-e, b-d-abs(f));
			h=(1-H)/2*max(-a-b+c+d-f,-b-d+abs(e))
			 +(1+H)/2*min( a+b-c-d-f, b+d-abs(e));
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM8r17RY() {return new M8r17("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r17WS() {return new M8r17("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r17MK() {return new M8r17("MK");}
	private static class M8r17 extends OrthogonalizedLieMarkovRateMatrix {
		public M8r17(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "8.17";
			super.nParam_ = 7;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD,
							                       SymmetryMatrix.smLt,
							                       SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, a, b, c, d, e, f, g, h;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(b-abs(c),a+b);
			e=E*(b+d);
			f=(1-F)/2*max(-a-2*b+abs(c-e),-b+d)
			 +(1+F)/2*min( a+2*b-abs(c+e), b-d);
			g=(1-G)/2*max(-b+d+abs(c+f),-a-b-d-e)
			 +(1+G)/2*min( b-d-abs(c-f), a+b+d-e);
			h=(1-H)/2*max(-b-d+abs(c+e),-a-b+d-f)
			 +(1+H)/2*min( b+d-abs(c-e), a+b-d-f);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM8r18RY() {return new M8r18("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r18WS() {return new M8r18("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM8r18MK() {return new M8r18("MK");}
	private static class M8r18 extends OrthogonalizedLieMarkovRateMatrix {
		public M8r18(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "8.18";
			super.nParam_ = 7;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD,
							                       SymmetryMatrix.smUtm,
							                       SymmetryMatrix.smDt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, a, b, c, d, e, f, g, h;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*min(b-abs(c),a+b);
			e=E/2*(a+2*b+2*d-abs(c));
			f=F/2*(a+2*b-2*d-abs(c));
			g=(1-G)/2*max(-b-d-e+abs(c),-a-b-d+e)
		     +(1+G)/2*min( b+d-e-abs(c), a+b+d+e);
			h=(1-H)/2*max(-b+d-f+abs(c),-a-b+d+f)
			 +(1+H)/2*min( b-d-f-abs(c), a+b-d+f);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM9r20aRY() {return new M9r20a("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM9r20aWS() {return new M9r20a("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM9r20aMK() {return new M9r20a("MK");}
	private static class M9r20a extends OrthogonalizedLieMarkovRateMatrix {
		public M9r20a(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "9.20a";
			super.nParam_ = 8;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD,
							                       SymmetryMatrix.smUtm,
							                       SymmetryMatrix.smDt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, I, a, b, c, d, e, f, g, h, i;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			I = parameters[startIndex+7]; maxp = Math.max(maxp, I);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*(b-abs(c));
			e=E*(a+b);
			f=F/2*(a+2*b+e-abs(c+d));
			g=G/2*(a+2*b-e-abs(c-d));
			h=(1-H)/2*max(-b-f+abs(c+d),-a-b-e+f)
			 +(1+H)/2*min( b-f-abs(c+d), a+b+e+f);
			i=(1-I)/2*max(-b-g+abs(c-d),-a-b+e+g)
			 +(1+I)/2*min( b-g-abs(c-d), a+b-e+g);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
			super.nativeParameters_[8]=i;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM9r20bRY() {return new M9r20b("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM9r20bWS() {return new M9r20b("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM9r20bMK() {return new M9r20b("MK");}
	private static class M9r20b extends OrthogonalizedLieMarkovRateMatrix {
		public M9r20b(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "9.20b";
			super.nParam_ = 8;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smLt,
							                       SymmetryMatrix.smRt,
							                       SymmetryMatrix.smUtm,
							                       SymmetryMatrix.smDt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, I, a, b, c, d, e, f, g, h, i;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			I = parameters[startIndex+7]; maxp = Math.max(maxp, I);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*min(b,a+2*b);
			e=E*(a+b);
			d=D*min(a+2*b,b-abs(c));
			f=F*min(b,a+2*b-abs(c));
			g=(1-G)/2*max(-b,-a-2*b+c+abs(d-e+f),-a-2*b-c+abs(d+e+f))
			 +(1+G)/2*min( b, a+2*b+c-abs(d+e-f), a+2*b-c-abs(d-e-f)); 
			h=(1-H)/2*max(-b+abs(c+d+g),-a-b-e+f)
			 +(1+H)/2*min( b-abs(c+d-g), a+b+e+f);
			i=(1-I)/2*max(-b+abs(c-d-f),-a-b+e+g)
			 +(1+I)/2*min( b-abs(c-d+f), a+b-e+g);
			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
			super.nativeParameters_[8]=i;
		}
	}

	public static OrthogonalizedLieMarkovRateMatrix getM10r12RY() {return new M10r12("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM10r12WS() {return new M10r12("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM10r12MK() {return new M10r12("MK");}
	private static class M10r12 extends OrthogonalizedLieMarkovRateMatrix {
		public M10r12(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "10.12";
			super.nParam_ = 9;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD,
							                       SymmetryMatrix.smUtm,
							                       SymmetryMatrix.smDt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, I, J, a, b, c, d, e, f, g, h, i, j;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			I = parameters[startIndex+7]; maxp = Math.max(maxp, I);
			J = parameters[startIndex+8]; maxp = Math.max(maxp, J);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			f=F*b;
			d=(1-D)/2*(-b+abs(c+f))+(1+D)/2*(b-abs(c-f));
			e=E*(a+b)-f;
			g=G*(a+2*b+e+2*f-abs(c+d))/2;
			h=H*(a+2*b-e-2*f-abs(c-d))/2;
			i=(1-I)/2*(max(-f-g+abs(c+d),-a-e-f+g)-b)
			 +(1+I)/2*(min( f-g-abs(c+d), a+e+f+g)+b);
			j=(1-J)/2*(max( f-h+abs(c-d),-a+e+f+h)-b)
			 +(1+J)/2*(min(-f-h-abs(c-d), a-e-f+h)+b);

			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
			super.nativeParameters_[8]=i;
			super.nativeParameters_[9]=j;
		}
	}
	
	public static OrthogonalizedLieMarkovRateMatrix getM10r34RY() {return new M10r34("RY");}
	public static OrthogonalizedLieMarkovRateMatrix getM10r34WS() {return new M10r34("WS");}
	public static OrthogonalizedLieMarkovRateMatrix getM10r34MK() {return new M10r34("MK");}
	private static class M10r34 extends OrthogonalizedLieMarkovRateMatrix {
		public M10r34(String permName) {
			super(permName);
			super.name_ = super.permutation_ + "10.34";
			super.nParam_ = 9;
			super.matrices_ = new SymmetryMatrix[]{SymmetryMatrix.smA,
					                               SymmetryMatrix.smB,
							                       SymmetryMatrix.smArt,
							                       SymmetryMatrix.smAr,
							                       SymmetryMatrix.smAt,
							                       SymmetryMatrix.smH,
							                       SymmetryMatrix.smU,
							                       SymmetryMatrix.smD,
							                       SymmetryMatrix.smLt,
							                       SymmetryMatrix.smRt};
			super.nativeParameters_ = new double[super.nParam_+1];
		}
		
		protected void translate(double[] parameters, int startIndex) {
			double B, C, D, E, F, G, H, I, J, a, b, c, d, e, f, g, h, i, j;
			
			B = parameters[startIndex  ]; double maxp= Math.abs(B);
			C = parameters[startIndex+1]; maxp = Math.max(maxp, C);
			D = parameters[startIndex+2]; maxp = Math.max(maxp, D);
			E = parameters[startIndex+3]; maxp = Math.max(maxp, E);
			F = parameters[startIndex+4]; maxp = Math.max(maxp, F);
			G = parameters[startIndex+5]; maxp = Math.max(maxp, G);
			H = parameters[startIndex+6]; maxp = Math.max(maxp, H);
			I = parameters[startIndex+7]; maxp = Math.max(maxp, I);
			J = parameters[startIndex+8]; maxp = Math.max(maxp, J);
			if (maxp > 1) throw new IllegalArgumentException("Parameter out of range [-1,1]");
			
			// 3b+a=1 gives a rate matrix with trace=-4. Combined with b>=0, a+b>=0 => 0<=b<=0.5 
			b=(B+1)/4;
			// Protection against rounding error making a+b<0:
			a=max(1-3*b,-b);
			c=C*b;
			d=D*b;
			f=(1-F)/2*(abs(c+d)-b)+(1+F)/2*(b-abs(c-d));
			e=E*(a+b)-f;
			g=G*min(b+f,a+2*b-abs(c),a+2*b+e,a+3*b-c-d-e-f,a+3*b+c+d-e-f);
			h=(1-H)/2*max(-b+f,-a-2*b+c+abs(d+e-g),-a-2*b-c+abs(d-e-g))
			 +(1+H)/2*min( b-f, a+2*b-c-abs(d+e+g), a+2*b+c-abs(d-e+g));
			i=(1-I)/2*max(-b+f+abs(c-d-h),-a-b-e-f-g)
			 +(1+I)/2*min( b-f-abs(c-d+h), a+b+e+f-g);
			j=(1-J)/2*max(-b-f+abs(c+d+g),-a-b+e+f-h)
			 +(1+J)/2*min( b+f-abs(c+d-g), a+b-e-f-h);

			super.nativeParameters_[0]=a;
			super.nativeParameters_[1]=b;
			super.nativeParameters_[2]=c;
			super.nativeParameters_[3]=d;
			super.nativeParameters_[4]=e;
			super.nativeParameters_[5]=f;
			super.nativeParameters_[6]=g;
			super.nativeParameters_[7]=h;
			super.nativeParameters_[8]=i;
			super.nativeParameters_[9]=j;
		}
	}


	public static void printModels() {
		new M3r4("RY").printObject(System.out);
	}
}
