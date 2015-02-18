// SingleClassSubstitutionModel.java
//
// (c) 1999-2004 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package pal.substmodel;

/**
 * <p>Title: SingleClassSubstitutionModel </p>
 * <p>Description: A SubstitutionModel class for the new style of rate matrix. Can use getAllParameters() to match up with rate matrix parameters.</p>
 * @author Matthew Goode
 */
import pal.misc.*;
import pal.math.*;
import pal.datatype.*;
import pal.eval.PatternInfo;

import java.io.*;

// MDW: now implements RateMatrixSubstitutionModel instead of SubstitutionModel
public class SingleClassSubstitutionModel extends Parameterized.ParameterizedUser implements RateMatrixSubstitutionModel {
	private RateMatrixHandler handler_;
	private DataType dataType_;

	private static final long serialVersionUID = 4938429359234234234L;

	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(1); //Version number
		out.writeObject(handler_);
		out.writeObject(dataType_);
	}

	private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
		byte version = in.readByte();
		switch(version) {
			default : {
				this.handler_ = (RateMatrixHandler)in.readObject();
				this.dataType_ = (DataType)in.readObject();
				setParameterizedBase(handler_);
				break;
			}
		}
	}

	private SingleClassSubstitutionModel(SingleClassSubstitutionModel toCopy) {
		this.handler_ = toCopy.handler_.getCopy();
		this.dataType_ = toCopy.dataType_;
		setParameterizedBase(handler_);
	}

	// MDW: added 'fixedFrequencies' parameter.
	public SingleClassSubstitutionModel(NeoRateMatrix base, DataType dt, double[] frequencies, boolean fixedFrequencies) {
		this.handler_ = new RateMatrixHandler(base, frequencies, fixedFrequencies);
		this.dataType_ = dt;
		setParameterizedBase(handler_);
	}
	public SingleClassSubstitutionModel(NeoRateMatrix base, DataType dt, double[] frequencies) {
		this(base, dt, frequencies, false);
	}
	// MDW: Added for debugging purposes. Perhaps require this in the interface? Use 'report' instead?
	public void printObject() { this.printObject(""); }
	public void printObject(String ps) {
		System.out.println(ps+"Object: SingleClassSubstititutionModel");
		handler_.printObject("  ");
		System.out.println(ps+"unfinished...");
	}

	// MDW: to implement RateMatrixSubstutionModel
	public void getRateMatrix(double[][] matrixStore) {
		handler_.getRateMatrix(matrixStore);
	}
	
	public DataType getDataType() {		return dataType_;		}
	public int getNumberOfTransitionCategories() {	return 1; 	}
	public double getTransitionCategoryProbability(int category) {	return 1;		}
	public double[] getTransitionCategoryProbabilities() {		return new double[] { 1 };		}

	public double[] getEquilibriumFrequencies() {	return handler_.getEquilibriumFrequencies();			}
	public boolean getFixedFrequencies() { return handler_.getFixedFrequencies(); }


	public void getTransitionProbabilities(double branchLength, ModelEdgeParameters edgeParams, double[][][] store) {
		handler_.getTransitionProbabilities(branchLength, store[0]);
	}
	public void getTransitionProbabilities(double branchLength, ModelEdgeParameters edgeParams, int category, double[][] store) {
		handler_.getTransitionProbabilities(branchLength, store);
	}
	public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParams, double[][][] store) {
		handler_.getTransitionProbabilitiesTranspose(branchLength, store[0]);
	}
	public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParams, int category, double[][] store) {
		handler_.getTransitionProbabilitiesTranspose(branchLength, store);
	}
	
	public void addPalObjectListener(PalObjectListener l) {
		throw new RuntimeException("Sorry, NeoRateMatrix stuff does not work with old likelihood calculators!");
	}
	public void removePalObjectListener(PalObjectListener l) {
		throw new RuntimeException("Sorry, NeoRateMatrix stuff does not work with old likelihood calculators!");
	}
	
	public boolean mustOptimiseMultipleParametersOnEdge() { return false; }
	
	public OrthogonalHints getOrthogonalHints() {		return null; 	 }
	public double likelihoodPenalty(PatternInfo pi, ModelEdgeParameters mep) { return 0; }
	public boolean isTimeReversible() { return true; }

	// interface Report
	public void report(PrintWriter out) {		handler_.report(out);		}
	public String toString() {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw,true);
		report(pw);
		return "Single Class Substitution Model:\n"+sw.toString();
	}
	public Object clone() {
		return new SingleClassSubstitutionModel(this);
	}
	public SubstitutionModel getCopy() {
		return new SingleClassSubstitutionModel(this);
	}
}