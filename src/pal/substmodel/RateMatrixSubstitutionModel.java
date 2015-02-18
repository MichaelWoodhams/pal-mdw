/**
 * 
 */
package pal.substmodel;

/**
 * For substitution models which are capable of returning a RateMatrix
 * 
 * @author woodhams
 *
 */
public interface RateMatrixSubstitutionModel extends SubstitutionModel {
	public void getRateMatrix(double[][] matrixStorage);
}
