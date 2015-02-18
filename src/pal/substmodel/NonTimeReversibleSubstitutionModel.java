package pal.substmodel;

import java.util.Arrays;

import pal.datatype.DataType;
import pal.misc.*;

/**
 * For non-time-reversible models, the appropriate 'equilibrium frequencies' to use will vary
 * according to where in the tree you are evaluating from. This interface adds a method
 * to set those frequencies.
 * 
 * @author woodhams
 *
 * TODO: We might want different root frequencies for different rate classes - this is not
 * currently supported, but wouldn't be difficult. 
 */

@SuppressWarnings("serial")
public abstract class NonTimeReversibleSubstitutionModel extends Parameterized.ParameterizedUser implements SubstitutionModel {
	private final double[] ones_; // vector of ones. Used as 'equilibrium frequency'.
	
	public NonTimeReversibleSubstitutionModel(DataType dt) {
		ones_ = new double[dt.getNumStates()];
		Arrays.fill(ones_, 1);
	}
	public final boolean isTimeReversible() { return false; }
	
	public abstract void setRootFrequencies(double[] freqs);
	public abstract double[] getRootFrequencies();
	public abstract void getRootFrequencies(double[] store);
	/**
	 * From now on, root frequencies will be considered parameters to the model.
	 * WARNING: This increases the number of model parameters. 
	 * Do not call this method if the NTRSModel is referenced by any
	 * sort of optimizer, as it will crash or fail to optimize the new parameters.
	 */
	public abstract void makeRootFrequenciesVariable();
	/**
	 * From now on, root frequencies will no longer be considered parameters to the model.
	 * WARNING: This decreases the number of model parameters. 
	 * Do not call this method if the NTRSModel is referenced by any
	 * sort of optimizer, as it will crash.
	 */
	public abstract void makeRootFrequenciesFixed();
	public abstract Object clone();
	/*
	 * For NTRSubModel, 'getEquilibriumFrequencies' must return a vector of ones.
	 * In the final likelihood calculation for an NTR model, we perform a 3-way dot-product
	 * of two conditional probability vectors with the equilibrium frequencies. For
	 * NTR models, the equivalent of the equilibrium frequencies are the root frequencies,
	 * but they must be applied at a specific place, which is taken care of elsewhere.
	 * By returning a vector of ones, we effectively nullify the 'multiply by equilibrium
	 * frequencies' part of the calculation for NTR models but not for TR models.
	 */ 
	public double[] getEquilibriumFrequencies() { 
		return ones_;
	}
	/**
	 * For models which use ModelEdgeParameters *and* are non-time-reversible:
	 * Allows for the option of the Root 'edge' (in PossiblyRootedMLSearcher) to have
	 * a special ModelEdgeParameter. Used by MosaicSubstitutionModel. For most other
	 * models, should return null. 
	 * @return
	 */
	public abstract ModelEdgeParameters getRootModelEdgeParameters();
}
