package pal.substmodel;

import java.io.PrintWriter;

import pal.datatype.DataType;
import pal.eval.PatternInfo;
import pal.math.OrthogonalHints;
import pal.math.StochasticVector;
import pal.misc.MultiParameterized;
import pal.misc.PalObjectListener;
import pal.misc.Parameterized;

/**
 * A simple example of a non-time-reversible model:
 * This is just a wrapper around a TR model, but with root frequency not equal
 * to equilibrium frequency. 
 * 
 * @author M Woodhams
 *
 */


@SuppressWarnings("serial")
public class GeneralNonEquilibriumModel extends NonTimeReversibleSubstitutionModel {
	private SubstitutionModel base_;
	private StochasticVector rootFreq_;
	private Parameterized parameterization_;
	
	/**
	 * Normal constructor
	 * @param baseModel Any substitution model (usually time reversible.)  
	 * @param rootFreq
	 * @param fixedFrequencies: if false, root frequencies are optimisable parameters of the model.
	 */
	public GeneralNonEquilibriumModel(SubstitutionModel baseModel, double[] rootFreq, boolean fixedFrequencies) {
		super(baseModel.getDataType());
		base_ = baseModel;
		rootFreq_ = new StochasticVector(rootFreq);
		if (fixedFrequencies) {
			makeRootFrequenciesFixed();
		} else {
			makeRootFrequenciesVariable();
		}
	}
	/*
	 * Simplified constructor: assume root frequencies are model parameters, start off with flat frequencies.
	 */
	public GeneralNonEquilibriumModel(SubstitutionModel baseModel) {
		super(baseModel.getDataType());
		base_ = baseModel;
		int n = base_.getDataType().getNumStates();
		double[] flat = new double[n];
		for (int i=0; i<n; i++) {
			flat[i] = 1./n;
		}
		rootFreq_ = new StochasticVector(flat);
		makeRootFrequenciesVariable();
	}
	private GeneralNonEquilibriumModel(GeneralNonEquilibriumModel toCopy) {
		super(toCopy.base_.getDataType());
		this.base_ = (SubstitutionModel)toCopy.base_.clone();
		this.rootFreq_ = new StochasticVector(toCopy.getRootFrequencies());
		if (toCopy.parameterization_ == toCopy.base_) {
			makeRootFrequenciesFixed();
		} else {
			makeRootFrequenciesVariable();
		}
	}
	
	
	public void makeRootFrequenciesVariable() {
		parameterization_ = new MultiParameterized(base_,rootFreq_);
		setParameterizedBase(parameterization_);
	}
	public void makeRootFrequenciesFixed() {
		parameterization_ = base_;
		setParameterizedBase(parameterization_);
	}
	
	public ModelEdgeParameters getRootModelEdgeParameters() { return null; }
	public void getRootFrequencies(double[] store) { rootFreq_.getVector(store); }
	public double[] getRootFrequencies()           { return rootFreq_.getVector(); }
	public void setRootFrequencies(double[] freqs) { rootFreq_.setVector(freqs); }
	public double likelihoodPenalty(PatternInfo pi, ModelEdgeParameters mep) { return 0; }
	public boolean mustOptimiseMultipleParametersOnEdge() { return false; }
	public Object clone() {
		return new GeneralNonEquilibriumModel(this);
	}
	
	/*
	 * Methods which pass through to the base substitution model:
	 */
	public void addPalObjectListener(PalObjectListener l)    { base_.addPalObjectListener(l); }
	public void removePalObjectListener(PalObjectListener l) { base_.removePalObjectListener(l); }
	public DataType getDataType()                { return base_.getDataType(); }
	public int getNumberOfTransitionCategories() { return base_.getNumberOfTransitionCategories(); }
	// Is this valid when parameterization might include root freqs? Orthogonal hints never seem to be used anyhow.
	public OrthogonalHints getOrthogonalHints()  { return base_.getOrthogonalHints(); } 
	public double[] getTransitionCategoryProbabilities()         { return base_.getTransitionCategoryProbabilities(); }
	public double getTransitionCategoryProbability(int category) { return base_.getTransitionCategoryProbability(category); }
	public void getTransitionProbabilities(double branchLength,	ModelEdgeParameters edgeParameters, double[][][] tableStore) 
	{
		base_.getTransitionProbabilities(branchLength, edgeParameters, tableStore);
	}
	public void getTransitionProbabilities(double branchLength,	ModelEdgeParameters edgeParameters, int category, double[][] tableStore) {
		base_.getTransitionProbabilities(branchLength, edgeParameters, category, tableStore);
	}
	public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParameters, double[][][] tableStore) {
		base_.getTransitionProbabilitiesTranspose(branchLength, edgeParameters, tableStore);
	}
	public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters edgeParameters, int category, double[][] tableStore) {
		base_.getTransitionProbabilitiesTranspose(branchLength, edgeParameters, category, tableStore);
	}


	public void report(PrintWriter out) { 
		out.printf("GeneralNonEquilibriumModel, base freqs = %s", pal.misc.Utils.toString(rootFreq_.getVector()));
		out.print(", base model = ");
		base_.report(out);
	}
}
