// SubstitutionModel.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.substmodel;

import pal.misc.*;
import pal.math.*;
import pal.datatype.*;
import pal.eval.PatternInfo;

import java.io.*;


/**
 * <b>model of sequence substitution (rate matrix + rate variation)</b>.
 * provides a convenient interface for the computation of transition probabilities
 *
 * @version $Id: SubstitutionModel.java,v 1.33 2004/05/19 04:05:21 matt Exp $
 *
 * @author Alexei Drummond
 * @author Matthew Goode
 */
public interface SubstitutionModel extends Parameterized, Report, java.io.Serializable {

	public DataType getDataType();
	public int getNumberOfTransitionCategories();
	public double getTransitionCategoryProbability(int category);
	/**
	 * @return all the category probabilites for each category respectively.
	 * @note Applications should not alter the returned array in
	 * any way!
	 */
	public double[] getTransitionCategoryProbabilities();
	/**
	 * Table is organized as [transition_group][from][to]
	 */
	//public void getTransitionProbabilities(double branchLength, double[][][] tableStore);
	/** 
	 * MDW added version of routine. 
	 * If your model doesn't have edge-specific parameters (beyond edge length), then this should just be a call to
	 * getTransitionProbabilities(double, double[][][])
	 */
	 
	public void getTransitionProbabilities(double branchLength, ModelEdgeParameters edgeParameters, double[][][] tableStore);
	/**
	 * Table is organized as [transition_group][to][from]
	 * i.e. columns sum to one.
	 */
	//public void getTransitionProbabilitiesTranspose(double branchLength, double[][][] tableStore);
	/**
	 * MDW: For most models, this will just call the corresponding routine without ModelEdgeParameters.
	 */
	public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParameters, double[][][] tableStore);

	/**
	 * Table is organized as [transition_group][from][to]
	 * i.e. rows sum to one.
	 */
	//public void getTransitionProbabilities(double branchLength, int category, double[][] tableStore);
	/**
	 * MDW: For most models, this will just call the corresponding routine without ModelEdgeParameters.
	 */
	public void getTransitionProbabilities(double branchLength, ModelEdgeParameters edgeParameters, int category, double[][] tableStore);

	/**
	 * Table is organized as [transition_group][to][from]
	 */
	//public void getTransitionProbabilitiesTranspose(double branchLength, int category, double[][] tableStore);
	/**
	 * MDW: For most models, this will just call the corresponding routine without ModelEdgeParameters.
	 */
	public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParameters, int category, double[][] tableStore);

	/**
	 * Should return a double[] array of the related equilibrium frequencies. As a rule, callers should not alter the returned array (it may be used internally)
	 */
	public double[] getEquilibriumFrequencies();

	public void addPalObjectListener(PalObjectListener l);
	public void removePalObjectListener(PalObjectListener l);

	/*
	 * If your model doesn't have edge-specific parameters (beyond edge length), then this should return false.
	 * If it does, it should return true only if they parameters should get optimized at the same time as edge lengths.
	 */
	public boolean mustOptimiseMultipleParametersOnEdge();
	
	/**
	 * May return null
	 */
	public OrthogonalHints getOrthogonalHints();

	public Object clone();
	
	/**
	 * For almost all models, this should return zero.
	 * The number returned is a penalty applied to likelihood calculations. This is to allow for AIC or BIC
	 * scores from models with variable numbers of parameters.
	 * @param pi (use this to find number of sites if trying to implement BIC)
	 * @param mep (used by the MosaicSubstitutionModel to count how many edge model weight parameters there are.)
	 * @return
	 */
	public double likelihoodPenalty(PatternInfo pi, ModelEdgeParameters mep);
	/**
	 * Is this a time reversible (and stationary) model? 
	 * (stationary implies a priori base probabilities equals equilibrium base frequencies.) 
	 * @return
	 */
	public boolean isTimeReversible();



	//===========================
	//===== Utils
	/**
	 * A small Utility class for things relating to Substitution Models in general
	 */
	public static class Utils {
		public static final double[][][] generateTransitionProbabilityTables(SubstitutionModel model) {
			int numberOfStates = model.getDataType().getNumStates();
			return new double[model.getNumberOfTransitionCategories()][numberOfStates][numberOfStates];
		}

		/**
		 * @return a substitution model base on a rate matrix. There is only one transition category.
		 * There is no independent distribution access (as not distribution across one transition category)
		 */
		public static final SubstitutionModel createSubstitutionModel(RateMatrix rm) {
			return new SimpleSubstitutionModel(rm);
		}
		/**
		 * @return a substitution model base on a rate matrix. There is only one transition category.
		 * There is no independent distribution access (as not distribution across one transition category)
		 */
		public static final SubstitutionModel createSubstitutionModel(NeoRateMatrix rm, DataType dt, double[] equilibriumFrequencies) {
			return new SingleClassSubstitutionModel(rm,dt, equilibriumFrequencies, true);
		}

		/**
		 * @return a substitution model base on a rate matrix, and a rate distribution. There are as many transition categories as there are rate categories in the rate distribution.
		 * There is no independent distribution access (as rate distributions don't normally allow changing of probabilities for a category without changing the rate of a category - at least none in pal do)
		 */
		public static final SubstitutionModel createSubstitutionModel(RateMatrix rm, RateDistribution rd) {
			return new RateDistributionSubstitutionModel(rm,rd);
		}

		/**
		 * @param parameteriseDistribution If true will include the distribution parameters as part of the
		 * substitution modle parameters, if false then the distribution parameters are set from the
		 * point of view of the substitution model.
		 * @return a substitution model base on a rate matrix, and a rate distribution. There are as many transition categories as there are rate categories in the rate distribution.
		 * There is no independent distribution access (as rate distributions don't normally allow changing of probabilities for a category without changing the rate of a category - at least none in pal do)
		 */
		public static final SubstitutionModel createSubstitutionModel(RateMatrix rm, RateDistribution rd, boolean parameteriseDistribution) {
			return new RateDistributionSubstitutionModel(rm,rd,parameteriseDistribution);
		}



//======== Private Inner classes
//==============================
		//========= SimpleSubstitutionModel ===============
		//=================================================
		private static class SimpleSubstitutionModel extends Parameterized.ParameterizedUser implements SubstitutionModel {
			private RateMatrix matrixBase_;

			private static final long serialVersionUID = 3054360219040005677L;

			private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
				out.writeByte(1); //Version number
				out.writeObject(matrixBase_);
			}

			private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
				byte version = in.readByte();
				switch(version) {
					default : {
						matrixBase_ = (RateMatrix)in.readObject();
						setParameterizedBase(matrixBase_);
						break;
					}
				}
			}


			private SimpleSubstitutionModel(SimpleSubstitutionModel toCopy) {
				this.matrixBase_ = (RateMatrix)toCopy.matrixBase_.clone();
				setParameterizedBase(matrixBase_);
			}

			public SimpleSubstitutionModel(RateMatrix base) {
				super(base);
				this.matrixBase_ = base;
			}

			public DataType getDataType() {
				return matrixBase_.getDataType();
			}
			public int getNumberOfTransitionCategories() {
				return 1;
			}
			public double getTransitionCategoryProbability(int category) {
				return 1;
			}
			public double[] getTransitionCategoryProbabilities() {
				return new double[] { 1 };
			}
			public double[] getEquilibriumFrequencies() {	return matrixBase_.getEquilibriumFrequencies();			}


			public void getTransitionProbabilities(double branchLength, ModelEdgeParameters edgeParams, double[][][] store) {
				matrixBase_.setDistance(branchLength);
				matrixBase_.getTransitionProbabilities(store[0]);
			}
			public void getTransitionProbabilities(double branchLength, ModelEdgeParameters edgeParams, int category, double[][] store) {
				matrixBase_.setDistance(branchLength);
				matrixBase_.getTransitionProbabilities(store);
			}
			public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParams, double[][][] store) {
				matrixBase_.setDistanceTranspose(branchLength);
				matrixBase_.getTransitionProbabilities(store[0]);
			}
			public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParams, int category, double[][] store) {
				matrixBase_.setDistanceTranspose(branchLength);
				matrixBase_.getTransitionProbabilities(store);
			}
			

			public void addPalObjectListener(PalObjectListener l) {
				matrixBase_.addPalObjectListener(l);
			}
			public void removePalObjectListener(PalObjectListener l) {
				matrixBase_.removePalObjectListener(l);
			}
			public OrthogonalHints getOrthogonalHints() {		return null; 	 }
			public double likelihoodPenalty(PatternInfo pi, ModelEdgeParameters mep) { return 0; }
			public boolean isTimeReversible() { return true; }
			
			public boolean mustOptimiseMultipleParametersOnEdge() { return false; }

			// interface Report
			public void report(PrintWriter out) {
				matrixBase_.report(out);
			}
			public String toString() {
				StringWriter sw = new StringWriter();
				PrintWriter pw = new PrintWriter(sw,true);
				report(pw);
				return "Simple Substitution Model:\n"+sw.toString();
			}
			public Object clone() {
				return new SimpleSubstitutionModel(this);
			}
			public SubstitutionModel getCopy() {
				return new SimpleSubstitutionModel(this);
			}
		}
		//========= SimpleSubstitutionModel ===============


		//============ RateDistributionSubstitutionModel ===================
		//======================================
		/**
		 * MDW: I think that this class is obsolete and you should use GeneralRateMatrixSubstitutionModel instead.
		 * (But I'm not so sure as to deprecate it.)
		 */
		private static class RateDistributionSubstitutionModel extends Parameterized.ParameterizedUser implements SubstitutionModel {
			private RateMatrix matrixBase_;
			private RateDistribution distribution_;
			private int numberOfDistributionCategories_;
			private boolean parameteriseDistribution_;

			private static final long serialVersionUID = -3530291767049646272L;

			private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
				out.writeByte(2); //Version number
				out.writeObject(matrixBase_);
				out.writeObject(distribution_);
				out.writeBoolean(parameteriseDistribution_);
			}

			private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
				byte version = in.readByte();
				switch(version) {
					default : {
						matrixBase_ = (RateMatrix)in.readObject();
						distribution_ = (RateDistribution)in.readObject();
						numberOfDistributionCategories_ = distribution_.getNumberOfRates();
						parameteriseDistribution_ = in.readBoolean();
						setParameterizedBase(parameteriseDistribution_ ? Parameterized.Utils.combine(new Parameterized[] {matrixBase_,distribution_}) : matrixBase_);
						break;
					}
					case 1 : {
						matrixBase_ = (RateMatrix)in.readObject();
						distribution_ = (RateDistribution)in.readObject();
						numberOfDistributionCategories_ = distribution_.getNumberOfRates();
						parameteriseDistribution_ = true;
						setParameterizedBase(parameteriseDistribution_ ? Parameterized.Utils.combine(new Parameterized[] {matrixBase_,distribution_}) : matrixBase_);
						break;
					}
				}
			}


			private RateDistributionSubstitutionModel(RateDistributionSubstitutionModel toCopy ) {
				this.matrixBase_ = (RateMatrix)toCopy.matrixBase_.clone();
				this.distribution_ = (RateDistribution)toCopy.distribution_.clone();
				this.parameteriseDistribution_ = toCopy.parameteriseDistribution_;
				this.numberOfDistributionCategories_ = distribution_.getNumberOfRates();
				setParameterizedBase(parameteriseDistribution_ ? Parameterized.Utils.combine(new Parameterized[] {matrixBase_,distribution_}) : matrixBase_);
			}
			public RateDistributionSubstitutionModel( RateMatrix base, RateDistribution distribution) {
				this(base,distribution,true);
			}
			public RateDistributionSubstitutionModel( RateMatrix base, RateDistribution distribution, boolean parameteriseDistribution) {
				super(parameteriseDistribution ? Parameterized.Utils.combine(new Parameterized[] {base,distribution}) : base);
				this.matrixBase_ = base;
				this.parameteriseDistribution_ = parameteriseDistribution;
				this.distribution_ = distribution;
				this.numberOfDistributionCategories_ = distribution_.getNumberOfRates();
			}

			public double[] getTransitionCategoryProbabilities() {
				return distribution_.probability;
			}
			public DataType getDataType() {
				return matrixBase_.getDataType();
			}
			public int getNumberOfTransitionCategories() {
				return distribution_.getNumberOfRates();
			}
			public double getTransitionCategoryProbability(int category) {
				return distribution_.probability[category];
			}
			public double[] getEquilibriumFrequencies() {
				return matrixBase_.getEquilibriumFrequencies();
			}

			public void getTransitionProbabilities(double branchLength, ModelEdgeParameters edgeParams, double[][][] store) {
				for(int i = 0 ; i < numberOfDistributionCategories_ ; i++) {
					matrixBase_.setDistance(branchLength*distribution_.rate[i]);
					matrixBase_.getTransitionProbabilities(store[i]);
				}

			}
			public void getTransitionProbabilities(double branchLength, ModelEdgeParameters edgeParams, int category, double[][] store) {
				matrixBase_.setDistance(branchLength*distribution_.rate[category]);
				matrixBase_.getTransitionProbabilities(store);
			}
			public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParams, double[][][] store) {
				for(int i = 0 ; i < numberOfDistributionCategories_ ; i++) {
					matrixBase_.setDistanceTranspose(branchLength*distribution_.rate[i]);
					matrixBase_.getTransitionProbabilities(store[i]);
				}
			}
			public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParams, int category, double[][] store) {
				matrixBase_.setDistanceTranspose(branchLength*distribution_.rate[category]);
				matrixBase_.getTransitionProbabilities(store);
			}
			public void addPalObjectListener(PalObjectListener l) {
				matrixBase_.addPalObjectListener(l);
				distribution_.addPalObjectListener(l);
			}
			public void removePalObjectListener(PalObjectListener l) {
				matrixBase_.removePalObjectListener(l);
				distribution_.removePalObjectListener(l);
			}
			public OrthogonalHints getOrthogonalHints() {		return null; 	 }
			public double likelihoodPenalty(PatternInfo pi, ModelEdgeParameters mep) { return 0; }
			public boolean isTimeReversible() { return true; }
			
			public boolean mustOptimiseMultipleParametersOnEdge() { return false; }

			public boolean isParameterBaseIncludingDistribution() { return parameteriseDistribution_; }

			// interface Report
			public void report(PrintWriter out) {
				matrixBase_.report(out);
				out.println();
				distribution_.report(out);
			}

			public String toString() {
				StringWriter sw = new StringWriter();
				PrintWriter pw = new PrintWriter(sw,true);
				report(pw);
				return "Substitution Model (with Rate Distribution):\n"+sw.toString();
			}
			public Object clone() {
				return new RateDistributionSubstitutionModel(this);
			}
			public SubstitutionModel getCopy() {
				return new RateDistributionSubstitutionModel(this);
			}
		}
	}
}
