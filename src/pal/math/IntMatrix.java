package pal.math;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

/**
 * Quick and dirty integer matrix functions (used by Paralinear distances.)
 * Optimised versions for size up to n=4.
 * Horribly inefficient for n>>4. 
 * Beware of overflow: if you use 'det' or 'adjugate', you will get answers on the order of 
 * (max absolute cell value)^(matrix size). E.g. a 4x4 matrix with elements of size up to 
 * 1000 could give a determinant on the order of 10^12. Data values are long ints 
 * (max value nearly 10^19). Still, for example, a 8x8 with element size ~1000 could easily overflow this.
 * If overflow is a risk, use real matrices instead.
 * 
 * For the 4x4 case, if all elements are less than 55000, we'll be OK.
 * 
 * @author Michael Woodhams
 *
 */

// TODO: See if any preexisting PAL code could be made more efficient by using this class.

public class IntMatrix {
	// Allow public access to data members for efficiency and laziness. 
	// It isn't likely we'll ever want to change the internal representation.
	public long[][] m;  
	public int rows;
	public int columns;
	
	public IntMatrix(int r, int c) {
		m = new long[r][c];
		rows = r;
		columns = c;
	}
	public IntMatrix(int n) {
		m = new long[n][n];
		rows = n;
		columns = n;
	}
	public IntMatrix(long[][] a) {
		rows = a.length;
		columns = a[0].length;
		m = a.clone();
	}
	public IntMatrix(int[][] a) {
		rows = a.length;
		columns = a[0].length;
		m = new long[rows][columns];
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {
				m[i][j] = a[i][j]; // casting int to long
			}
		}
	}

	/** Print the matrix to stdout.
     * @param w    Column width (per integer. Plus two spaces between columns.)
     */

   public void print (int w) {
      print(new PrintWriter(System.out,true),w); 
   }

   /** Print the matrix to the output stream.  Line the elements up in columns.
   @param output the output stream.
   @param width  Column width.
   */

   public void print (PrintWriter output, int width) {
	  String format = "%"+width+"d  ";
      output.println();  // start on new line.
      for (int i = 0; i < rows; i++) {
         for (int j = 0; j < columns; j++) {
        	 output.printf(format, m[i][j]);
         }
         output.println();
      }
      output.println();   // end with blank line.
   }

	public String toString() {
		StringBuffer output = new StringBuffer("[[");
	    for (int i = 0; i < rows; i++) {
	    	for (int j = 0; j < columns; j++) {
	    		output.append(m[i][j]);
	    		if (j < columns-1) { output.append(", ");}
	    	}
	    	if (i < rows-1) { output.append("], [");}
	    }
	    output.append("]]");
	    return output.toString();
	}
	
	// Explict formulae for determinants up to dimension 4. No size checking is done, for efficiency.
	public long det1() {
		return m[0][0] * m[0][0]; // Not sure this is mathematically correct, but who'd need to know anyway?
	}
	public long det2() {
		return m[0][0] * m[1][1] - m[0][1] * m[1][0];
	}
	
	public long det3() {
		return m[0][0]*m[1][1]*m[2][2] 
		     + m[0][1]*m[1][2]*m[2][0] 
		     + m[0][2]*m[1][0]*m[2][1] 
		     - m[0][0]*m[1][2]*m[2][1] 
		     - m[0][1]*m[1][0]*m[2][2] 
		     - m[0][2]*m[1][1]*m[2][0];
	}
	public long det4() {
		return
		+m[0][3]*m[1][2]*m[2][1]*m[3][0]
		-m[0][2]*m[1][3]*m[2][1]*m[3][0]
		-m[0][3]*m[1][1]*m[2][2]*m[3][0]
        +m[0][1]*m[1][3]*m[2][2]*m[3][0]
        +m[0][2]*m[1][1]*m[2][3]*m[3][0]
        -m[0][1]*m[1][2]*m[2][3]*m[3][0]
        -m[0][3]*m[1][2]*m[2][0]*m[3][1]
        +m[0][2]*m[1][3]*m[2][0]*m[3][1]
        +m[0][3]*m[1][0]*m[2][2]*m[3][1]
        -m[0][0]*m[1][3]*m[2][2]*m[3][1]
        -m[0][2]*m[1][0]*m[2][3]*m[3][1]
        +m[0][0]*m[1][2]*m[2][3]*m[3][1]
        +m[0][3]*m[1][1]*m[2][0]*m[3][2]
        -m[0][1]*m[1][3]*m[2][0]*m[3][2]
        -m[0][3]*m[1][0]*m[2][1]*m[3][2]
        +m[0][0]*m[1][3]*m[2][1]*m[3][2]
        +m[0][1]*m[1][0]*m[2][3]*m[3][2]
        -m[0][0]*m[1][1]*m[2][3]*m[3][2]
        -m[0][2]*m[1][1]*m[2][0]*m[3][3]
        +m[0][1]*m[1][2]*m[2][0]*m[3][3]
        +m[0][2]*m[1][0]*m[2][1]*m[3][3]
        -m[0][0]*m[1][2]*m[2][1]*m[3][3]
        -m[0][1]*m[1][0]*m[2][2]*m[3][3]
        +m[0][0]*m[1][1]*m[2][2]*m[3][3];
	}

	// Recursive method, adapted from http://wiki.answers.com/Q/Determinant_of_matrix_in_java	
	// by Avishek Ghosh. This won't be as fast as it could be, but the non-4x4 case isn't
	// used often anyhow.
	// Given that we're still dealing with integers, we don't have to worry
	// about roundoff, so a naive algorithm is OK.	    
	public long det() {

		int result = 0;

		switch (rows) {
			case 1: return this.det1();
			case 2: return this.det2();
			case 3: return this.det3();
			case 4: return this.det4();
			default:
		}

		for(int i = 0; i < rows; i++) {
			IntMatrix temp = new IntMatrix(rows-1);
			for(int j = 1; j < rows; j++) {
				for(int k = 0; k < rows; k++) {
					if(k < i) {
						temp.m[j-1][k] = m[j][k];
					} else if(k > i) {
						temp.m[j-1][k-1] = m[j][k];
					}
				}
			}
			result += m[0][i] * (i%2==0 ? 1 : -1) * temp.det();
		}
		return result;
	} 

	// The adjugate of a matrix is the inverse multiplied by the determinant (except that it also exists for singular matrices.)
	// Using this rather than inverse allows us to keep everything integer.
	public IntMatrix adjugate1() {
		IntMatrix adj = new IntMatrix(1);
		adj.m[0][0] = m[0][0];
		return adj;
	}
	public IntMatrix adjugate2() {
		IntMatrix adj = new IntMatrix(2);
		adj.m[0][0] =  m[1][1];
		adj.m[0][1] = -m[0][1];
		adj.m[1][0] = -m[1][0];
		adj.m[1][1] =  m[0][0];
		return adj;
	}
	public IntMatrix adjugate3() {
		IntMatrix adj = new IntMatrix(3);
		adj.m[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
		adj.m[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
		adj.m[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
		adj.m[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
		adj.m[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];
		adj.m[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
		adj.m[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
		adj.m[2][1] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
		adj.m[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
		return adj;
	}
	public IntMatrix adjugate4() {
		 IntMatrix adj = new IntMatrix(4);
		 adj.m[0][0]=m[1][2]*m[2][3]*m[3][1]-m[1][3]*m[2][2]*m[3][1]+m[1][3]*m[2][1]*m[3][2]-m[1][1]*m[2][3]*m[3][2]-m[1][2]*m[2][1]*m[3][3]+m[1][1]*m[2][2]*m[3][3];
		 adj.m[0][1]=m[0][3]*m[2][2]*m[3][1]-m[0][2]*m[2][3]*m[3][1]-m[0][3]*m[2][1]*m[3][2]+m[0][1]*m[2][3]*m[3][2]+m[0][2]*m[2][1]*m[3][3]-m[0][1]*m[2][2]*m[3][3];
		 adj.m[0][2]=m[0][2]*m[1][3]*m[3][1]-m[0][3]*m[1][2]*m[3][1]+m[0][3]*m[1][1]*m[3][2]-m[0][1]*m[1][3]*m[3][2]-m[0][2]*m[1][1]*m[3][3]+m[0][1]*m[1][2]*m[3][3];
		 adj.m[0][3]=m[0][3]*m[1][2]*m[2][1]-m[0][2]*m[1][3]*m[2][1]-m[0][3]*m[1][1]*m[2][2]+m[0][1]*m[1][3]*m[2][2]+m[0][2]*m[1][1]*m[2][3]-m[0][1]*m[1][2]*m[2][3];
		 adj.m[1][0]=m[1][3]*m[2][2]*m[3][0]-m[1][2]*m[2][3]*m[3][0]-m[1][3]*m[2][0]*m[3][2]+m[1][0]*m[2][3]*m[3][2]+m[1][2]*m[2][0]*m[3][3]-m[1][0]*m[2][2]*m[3][3];
		 adj.m[1][1]=m[0][2]*m[2][3]*m[3][0]-m[0][3]*m[2][2]*m[3][0]+m[0][3]*m[2][0]*m[3][2]-m[0][0]*m[2][3]*m[3][2]-m[0][2]*m[2][0]*m[3][3]+m[0][0]*m[2][2]*m[3][3];
		 adj.m[1][2]=m[0][3]*m[1][2]*m[3][0]-m[0][2]*m[1][3]*m[3][0]-m[0][3]*m[1][0]*m[3][2]+m[0][0]*m[1][3]*m[3][2]+m[0][2]*m[1][0]*m[3][3]-m[0][0]*m[1][2]*m[3][3];
		 adj.m[1][3]=m[0][2]*m[1][3]*m[2][0]-m[0][3]*m[1][2]*m[2][0]+m[0][3]*m[1][0]*m[2][2]-m[0][0]*m[1][3]*m[2][2]-m[0][2]*m[1][0]*m[2][3]+m[0][0]*m[1][2]*m[2][3];
		 adj.m[2][0]=m[1][1]*m[2][3]*m[3][0]-m[1][3]*m[2][1]*m[3][0]+m[1][3]*m[2][0]*m[3][1]-m[1][0]*m[2][3]*m[3][1]-m[1][1]*m[2][0]*m[3][3]+m[1][0]*m[2][1]*m[3][3];
		 adj.m[2][1]=m[0][3]*m[2][1]*m[3][0]-m[0][1]*m[2][3]*m[3][0]-m[0][3]*m[2][0]*m[3][1]+m[0][0]*m[2][3]*m[3][1]+m[0][1]*m[2][0]*m[3][3]-m[0][0]*m[2][1]*m[3][3];
		 adj.m[2][2]=m[0][1]*m[1][3]*m[3][0]-m[0][3]*m[1][1]*m[3][0]+m[0][3]*m[1][0]*m[3][1]-m[0][0]*m[1][3]*m[3][1]-m[0][1]*m[1][0]*m[3][3]+m[0][0]*m[1][1]*m[3][3];
		 adj.m[2][3]=m[0][3]*m[1][1]*m[2][0]-m[0][1]*m[1][3]*m[2][0]-m[0][3]*m[1][0]*m[2][1]+m[0][0]*m[1][3]*m[2][1]+m[0][1]*m[1][0]*m[2][3]-m[0][0]*m[1][1]*m[2][3];
		 adj.m[3][0]=m[1][2]*m[2][1]*m[3][0]-m[1][1]*m[2][2]*m[3][0]-m[1][2]*m[2][0]*m[3][1]+m[1][0]*m[2][2]*m[3][1]+m[1][1]*m[2][0]*m[3][2]-m[1][0]*m[2][1]*m[3][2];
		 adj.m[3][1]=m[0][1]*m[2][2]*m[3][0]-m[0][2]*m[2][1]*m[3][0]+m[0][2]*m[2][0]*m[3][1]-m[0][0]*m[2][2]*m[3][1]-m[0][1]*m[2][0]*m[3][2]+m[0][0]*m[2][1]*m[3][2];
		 adj.m[3][2]=m[0][2]*m[1][1]*m[3][0]-m[0][1]*m[1][2]*m[3][0]-m[0][2]*m[1][0]*m[3][1]+m[0][0]*m[1][2]*m[3][1]+m[0][1]*m[1][0]*m[3][2]-m[0][0]*m[1][1]*m[3][2];
		 adj.m[3][3]=m[0][1]*m[1][2]*m[2][0]-m[0][2]*m[1][1]*m[2][0]+m[0][2]*m[1][0]*m[2][1]-m[0][0]*m[1][2]*m[2][1]-m[0][1]*m[1][0]*m[2][2]+m[0][0]*m[1][1]*m[2][2];
		 return adj;
	}
	public IntMatrix adjugate() {
		if (rows != columns) {
			throw new IllegalArgumentException("Can only find adjugate of a square matrix");
		}
		switch (rows) {
			case 1: return this.adjugate1();
			case 2: return this.adjugate2();
			case 3: return this.adjugate3();
			case 4: return this.adjugate4();
			default:
		}
		// TODO: implement the general case.
		throw new RuntimeException("Adjugate not yet implemented for matrices larger than size 4");
	} 

	// Matrix multiplication. Special routines for square matrices up to size four.
	public IntMatrix mult1(IntMatrix that) {
		IntMatrix answer = new IntMatrix(rows);
		answer.m[0][0] = this.m[0][0] * that.m[0][0];
		return answer;
	}
	public IntMatrix mult2(IntMatrix that) {
		IntMatrix answer = new IntMatrix(rows);
		long[][] a = this.m;
		long[][] b = that.m;
		long[][] c = answer.m;
		c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0];
		c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1];
		c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0];
		c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1];
		return answer;
	}
	public IntMatrix mult3(IntMatrix that) {
		IntMatrix answer = new IntMatrix(rows);
		long[][] a = this.m;
		long[][] b = that.m;
		long[][] c = answer.m;
		c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
		c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
		c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
		c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
		c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
		c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
		c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
		c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
		c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
		return answer;
	}
	public IntMatrix mult4(IntMatrix that) {
		IntMatrix answer = new IntMatrix(rows);
		long[][] a = this.m;
		long[][] b = that.m;
		long[][] c = answer.m;
		c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0] + a[0][3]*b[3][0];
		c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1] + a[0][3]*b[3][1];
		c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2] + a[0][3]*b[3][2];
		c[0][3] = a[0][0]*b[0][3] + a[0][1]*b[1][3] + a[0][2]*b[2][3] + a[0][3]*b[3][3];
		c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0] + a[1][3]*b[3][0];
		c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1] + a[1][3]*b[3][1];
		c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2] + a[1][3]*b[3][2];
		c[1][3] = a[1][0]*b[0][3] + a[1][1]*b[1][3] + a[1][2]*b[2][3] + a[1][3]*b[3][3];
		c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0] + a[2][3]*b[3][0];
		c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1] + a[2][3]*b[3][1];
		c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2] + a[2][3]*b[3][2];
		c[2][3] = a[2][0]*b[0][3] + a[2][1]*b[1][3] + a[2][2]*b[2][3] + a[2][3]*b[3][3];
		c[3][0] = a[3][0]*b[0][0] + a[3][1]*b[1][0] + a[3][2]*b[2][0] + a[3][3]*b[3][0];
		c[3][1] = a[3][0]*b[0][1] + a[3][1]*b[1][1] + a[3][2]*b[2][1] + a[3][3]*b[3][1];
		c[3][2] = a[3][0]*b[0][2] + a[3][1]*b[1][2] + a[3][2]*b[2][2] + a[3][3]*b[3][2];
		c[3][3] = a[3][0]*b[0][3] + a[3][1]*b[1][3] + a[3][2]*b[2][3] + a[3][3]*b[3][3];		
		return answer;
	}

	// The general case
	public IntMatrix mult(IntMatrix that) {
		if (this.columns != that.rows) {
			throw new IllegalArgumentException("Dimension mismatch");
		}
		// The common square matrix case:
		if (this.rows == this.columns && that.rows == that.columns) {
			switch (this.rows) {
				case 1: return this.mult1(that);
				case 2: return this.mult2(that);
				case 3: return this.mult3(that);
				case 4: return this.mult4(that);
				default:
			}
		}
		IntMatrix answer = new IntMatrix(this.rows, that.columns);
		for (int i = 0; i < this.rows; i++) {
			for (int j = 0; j < that.columns; j++) {
				long sum = 0;
				for (int k = 0; k < this.columns; k++) {
					sum += this.m[i][k] * that.m[k][j];
				}
				answer.m[i][j] = sum;
			}
		}
		return answer;
	}
	
	public int[][] getMatrixAsShortInts() {
		int[][] answer = new int[rows][columns];
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {
				long value = m[i][j];
				if (value < Integer.MIN_VALUE || value > Integer.MAX_VALUE) {
					throw new RuntimeException("Integer overflow");
				} else {
					answer[i][j] = (int)value;
				}
			}
		}
		return answer;
	}
	
	public IntMatrix transpose() {
		IntMatrix t = new IntMatrix(columns,rows);
		for (int i=0; i<rows; i++) {
			for (int j=0; j<columns; j++) {
				t.m[j][i] = this.m[i][j];
			}
		}
		return t;
	}
	
	// Unit tests for this class.
	public static boolean test() {
		return test(0);
	}
	public static boolean test(long seed) {
		boolean success;
		success = test(1,seed);
		success = test(2,seed) && success;
		success = test(3,seed) && success;
		success = test(4,seed) && success;
		// Won't work yet:
		// success = test(5,1) && success;
		System.out.println("Warning: IntMatrix not fully implemented for square matrices of size 5 or larger.");
		return success;
	}
	// Seed = 0 used system-clock based seed. Any other value will give repeatable results.
	public static boolean test(int n, long seed) {
		boolean success = true;
		MersenneTwisterFast rng;
		if (seed == 0) {
			rng = new MersenneTwisterFast();
		} else {
			rng = new MersenneTwisterFast(seed);
		}
		// Make two random matrices filled with values -128-127
		IntMatrix a = new IntMatrix(n);
		IntMatrix b = new IntMatrix(n);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				a.m[i][j] = rng.nextByte();
				b.m[i][j] = rng.nextByte();
			}
		}
		long detA = a.det();
		IntMatrix adjA = a.adjugate();
		// should be identity * determinant:
		IntMatrix id = a.mult(adjA);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if (id.m[i][j] != (i==j ? detA : 0)) {
					System.out.printf("IntMatrix test 1 fails for n=%d, seed =%d, i=%d, j=%d\n",n,seed,i,j);
					success = false;
				}
			}
		}
		// And the other way around:
		id = adjA.mult(a);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if (id.m[i][j] != (i==j ? detA : 0)) {
					System.out.printf("IntMatrix test 2 fails for n=%d, seed =%d, i=%d, j=%d\n",n,seed,i,j);
					success = false;
				}
			}
		}
		IntMatrix ab = a.mult(b);
		id = adjA.mult(ab);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if (id.m[i][j] != detA*b.m[i][j]) {
					System.out.printf("IntMatrix test 3 fails for n=%d, seed =%d, i=%d, j=%d\n",n,seed,i,j);
					success = false;
				}
			}
		}
		
		return success;
	}
}
