package pal.substmodel;

import java.io.*;
import java.util.*;

import pal.algorithmics.StoppingCriteria;
import pal.algorithmics.UndoableAction;
import pal.alignment.*;
import pal.datatype.*;
import pal.eval.PatternInfo;
import pal.math.*;
import pal.misc.*;
import pal.tree.*;
import pal.treesearch.*;


/**
 * Unofficial modification to PAL by Michael Woodhams (School of IT, University of Sydney)
 * 
 * This is the scenario the MosaicSubstitutionModel is intended to address:
 * We are attempting to build a deep phylogeny. Not all taxa and not all genes evolve under the
 * same model. We select (usually overlapping) subsets of the taxa which appear to evolve under
 * a single model - giving us a set of (model, taxon subset) pairs. We now want to calculate the
 * maximum likelihood of some tree topology on the full set of taxa. The tree topology induces
 * a subtree topology for each taxon subset. The full tree can then be considered a supertree
 * of these subtrees. However, a given edge on the full tree may correspond to edges on several
 * of the subtrees, each with its own evolution model. To calculate the likelihood for such an
 * edge, we need to weight each of the models (rate matrices) which apply to that edge. The MosaicSubstitutionModel
 * allows such weights to be applied to each edge, and to calculate the resulting likelihood.
 * The weights are held in MosaicModelEdgeParameters object for each edge, and typically
 * are to be optimised. (Model optimisation for the MosaicSubstitutionModel optimises these
 * but leaves the submodels alone.) 
 * 
 * The MosaicSubstitutionModel holds the set of models from the subtrees and the list of the
 * MosaicModelEdgeParameters from all the edges in the full tree. (This is an undesirable
 * violation of encapsulation: the SubstitutionModel object is no longer divorced from the tree.
 * However, avoiding this is difficult.) An upshot of this is that you must create a new
 * MosaicSubstitutionModel for each tree you want to analyse.
 * 
 * Model 0 is special: when searching tree topologies, any edge which otherwise would have no model
 * apply to it gets model 0. 
 * 
 * I'm using GeneralRateDistributionSubstitutionModel as a guide for writing this class.
 * 
 * IMPORTANT ASSUMPTION: 
 * The submodels are not allowed to change once they have been put into a MosaicSubstitutionModel, except
 * via the methods supplied with this class. Otherwise we have no way of knowing whether the rate matrices are up to date.
 * 
 * MosaicSubstitutionModel and the MosaicModelEdgeParameters are interlinked. The current strategy is
 * to first create the MosaicSubstitutionModel, then add MosaicModelEdgeParameters to it one at a time.
 * 
 * TODO: Update this comment with some instructions on typical use cases.
 * 
 * @author Michael Woodhams
 *
 */
public class MosaicSubstitutionModel extends NonTimeReversibleSubstitutionModel {
	private int numberOfSubmodels_;
	private RateMatrixSubstitutionModel[] models_;
	private long[] submodelCurrentVersion_; // increment every time a rate matrix changes.
	private double[][][] submodelRateMatrices_; //MosaicModelEdgeParameters objects will have pointers into this array (pointers to a double[][]).
	private boolean[] submodelChanged_; // do we need to update submodelRateMatrices_?
	private StochasticVector rootFrequencies_;
	private DataType dataType_;
	private MultiParameterized parameterization_; // collects all the model parameters together.
	private SubmodelParameterAccessWatcher watcher_; // to recognize when submodelChanged_ needs to be updated.
	private boolean fixedRootFrequencies_; // if false, root frequencies are parameters of the model, else they only change by setRootFrequencies() 
	private RateDistribution rateDistribution_; // can be null, which indicates uniform distribution.
	
	// if data members change, need to update writeObject, readObject, constructor
	
	private static final long serialVersionUID = 287293574181L;

	/**
	 * Copies a MosaicSubstitutionModel, except for any attached tree. I.e. makes a new treeless MSM 
	 * with the same submodels.
	 * @param toCopy
	 */
	private MosaicSubstitutionModel(MosaicSubstitutionModel toCopy) {
		this((SingleClassSubstitutionModel[])toCopy.models_.clone(),toCopy.getRootFrequencies(),toCopy.fixedRootFrequencies_);
	}
	
	/**
	 * Simplified constructor:
	 * Root frequencies are not model parameters, are equal to model 0's root frequencies
	 * (at the time that this constructor is run: they do not change if model 0 changes.)
	 * @param models The list of submodels.
	 */
	public MosaicSubstitutionModel(SingleClassSubstitutionModel[] models) {
		this(models,models[0].getEquilibriumFrequencies(),true);
	}

	/**
	 * Basic constructor for uniform rates-over-sites Mosaic model
	 * @param models The submodels. (Any rates-over-sites in the submodels are ignored)
	 * @param initialRootFreq Initial root frequencies.
	 * @param fixedRootFrequencies If false, root frequencies are parameters of the model, else they only change by setRootFrequencies()
	 */
	public MosaicSubstitutionModel(RateMatrixSubstitutionModel[] models, double[] initialRootFreq, boolean fixedRootFrequencies) {
		this(models,new UniformRate(),initialRootFreq,fixedRootFrequencies);
	}
	
	/**
	 * Basic constructor for rates over sites Mosaic model. Note that the rate distribution given to this constructor
	 * is the only one that matters - the submodels may or may not have rate distributions, any they do have
	 * are ignored.
	 * @param models The submodels. 
	 * @param rateDist: the RateDistribution object
	 * @param initialRootFreq Initial root frequencies.
	 * @param fixedRootFrequencies If false, root frequencies are parameters of the model, else they only change by setRootFrequencies()
	 */
	// MDW: Change argument to 'RateMatrixSubstitutionClass[]' when I've coded support for this 
	public MosaicSubstitutionModel(RateMatrixSubstitutionModel[] models, RateDistribution rateDist, double[] initialRootFreq, boolean fixedRootFrequencies) {
		super(models[0].getDataType()); // NonTimeReversibleSubstitutionModel
		models_ = models;
		numberOfSubmodels_ = models.length;
		dataType_ = models[0].getDataType();
		submodelCurrentVersion_ = new long[numberOfSubmodels_];
		submodelChanged_ = new boolean[numberOfSubmodels_];
		for (int i=0; i<numberOfSubmodels_; i++) {
			submodelChanged_[i] = true; // all submodelRateMatrices are initially 'out of date'
			submodelCurrentVersion_[i] = 0;
		}
		int numStates = dataType_.getNumStates();
		// to start with a reasonable default for root frequencies, but these should be manually set later, 
		// or fitted as parameters.
		rootFrequencies_ = new StochasticVector(initialRootFreq);
		fixedRootFrequencies_ = fixedRootFrequencies;
		submodelRateMatrices_ = new double[numberOfSubmodels_][numStates][numStates];
		updateSubmodelRateMatrices();
		watcher_ = new SubmodelParameterAccessWatcher();
		rateDistribution_ = rateDist;
		if (fixedRootFrequencies) {
			makeRootFrequenciesFixed();
		} else {
			makeRootFrequenciesVariable();
		}
	}
	
	public void makeRootFrequenciesVariable() {
		Parameterized[] parameterizedBases = new Parameterized[models_.length+2];
		System.arraycopy(models_, 0, parameterizedBases, 0, models_.length);
		parameterizedBases[models_.length] = rateDistribution_;
		parameterizedBases[models_.length+1] = rootFrequencies_;

		parameterization_ = new MultiParameterized(parameterizedBases, watcher_);
		setParameterizedBase(parameterization_);
	}
	public void makeRootFrequenciesFixed() {
		Parameterized[] parameterizedBases = new Parameterized[models_.length+1];
		System.arraycopy(models_, 0, parameterizedBases, 0, models_.length);
		parameterizedBases[models_.length] = rateDistribution_;
		parameterization_ = new MultiParameterized(parameterizedBases, watcher_);
		setParameterizedBase(parameterization_);
	}
	
	/*
	 * Uses AIC correction: penalize log likelihood by 1 for every free edge weight parameter.
	 */
	public double likelihoodPenalty(PatternInfo pi, ModelEdgeParameters edgeParams) {
		if (!(edgeParams instanceof MosaicModelEdgeParameters)) {
			throw new RuntimeException("Bad edgeParams");
		}
		return ((MosaicModelEdgeParameters)edgeParams).getTreeEdgeParametersCount();
	}
	
	/*
	 * The methods and inner class relating to keeping the submodel rate matrices up to date.
	 * Because we create 'watcher_' and pass it to the MultiParameterized constructor, it will call
	 * watcher_'s parameterSet model whenever a parameter gets set, which lets us know which rate
	 * matrices are out of date.
	 */
	
	public void updateSubmodelRateMatrices() {
		for (int model = 0; model < numberOfSubmodels_; model++) {
			if (submodelChanged_[model]) {
				models_[model].getRateMatrix(submodelRateMatrices_[model]);
				submodelChanged_[model] = false;
				submodelCurrentVersion_[model]++;
			}
		}
	}
	
	public long getSubmodelCurrentVersion(int submodel) {
		return submodelCurrentVersion_[submodel];
	}
	
	private class SubmodelParameterAccessWatcher implements MultiParameterized.ParameterAccessWatcher {
		public void parameterSet(Parameterized baseParameterized, double param, int localParameter) {
			boolean found = false;
			for (int model = 0; !found && model < numberOfSubmodels_; model++) {
				if (baseParameterized == models_[model]) {
					found = true;
					submodelChanged_[model] = true;
				}
			}
			if (!found && baseParameterized instanceof RateMatrixSubstitutionModel) {
				// can't happen.
				throw new RuntimeException("Unrecognized submodel passed to SubmodelParameterAccessWatcher.parameterSet");
			}
		}
	}
	
	// used by MosaicModelEdgeParameters
	public double[][] getSubmodelRateMatrix(int index) {
		return submodelRateMatrices_[index];
	}
		
	/**
	 * Apply labels to each taxon saying which evolution models apply to it - in preparation for
	 * a MosaicSubstitutionModel analysis.
	 * 
	 * A front end to labelLeavesWithModels(Tree,SetOfSmallIntegers[]) or to labelLeavesWithWeightedModels(Tree,SetOfSmallIntegers)
	 * depending on useAllModelsWeighted flag.
	 * 
	 * @param tree
	 * @param taxaSubsets taxaSubset[submodel][taxon] is whether 'submodel' applies to 'taxon', where taxon are 
	 * in the order returned by tree.getExternalNode
	 * @param useAllModelsWeighted If true, all leaves get all models, but the initial edge weights will favour the 
	 * models listed in taxaSubsets.
	 */
	public void labelLeavesWithModels(Tree tree, boolean[][] taxaSubsets, boolean useAllModelsWeighted) {
		SetOfSmallIntegers[] newSubsets = new SetOfSmallIntegers[taxaSubsets.length];
		for (int model=0; model<taxaSubsets.length; model++) {
			newSubsets[model] = new SetOfSmallIntegers(taxaSubsets[model]);
		}
		labelLeavesWithModels(tree, newSubsets, useAllModelsWeighted);
	}
	
	/**
	 * Apply labels to each taxon saying which evolution models apply to it - in preparation for
	 * a MosaicSubstitutionModel analysis.
	 * 
	 * @param tree
	 * @param taxaSubsetsIDs taxaSubsetsIDs[submodel] is the list of taxa to which submodel applies.
	 * @param useAllModelsWeighted If true, all leaves get all models, but the initial edge weights will favour the 
	 * models listed in taxaSubsets.
	 */	public void labelLeavesWithModels(Tree tree, IdGroup[] taxaSubsetsIDs, boolean useAllModelsWeighted) {
		int nModels = taxaSubsetsIDs.length;
		SetOfSmallIntegers taxaSubsets[] = new SetOfSmallIntegers[nModels];
		for (int model=0; model < nModels; model++) {
			taxaSubsets[model] = new SetOfSmallIntegers();
			for (int taxon=0; taxon < taxaSubsetsIDs[model].getIdCount(); taxon++) {
				Identifier id = taxaSubsetsIDs[model].getIdentifier(taxon);
				int taxonNumber = tree.whichIdNumber(id.getName());
				taxaSubsets[model].addToSet(taxonNumber);
			}
		}
		labelLeavesWithModels(tree,taxaSubsets, useAllModelsWeighted);
	}
	
	/**
	 * Apply labels to each taxon saying which evolution models apply to it - in preparation for
	 * a MosaicSubstitutionModel analysis.
	 * 
	 * @param tree
	 * @param taxaSubsets taxaSubset[submodel] is the set of taxa which 'submodel' applies to, where 
	 * taxa are in the order returned by tree.getExternalNode
	 * @param useAllModelsWeighted If true, all leaves get all models, but the initial edge weights will favour the 
	 * models listed in taxaSubsets.
	 */
	 public void labelLeavesWithModels(Tree tree, SetOfSmallIntegers[] taxaSubsets, boolean useAllModelsWeighted) {
		 if (useAllModelsWeighted) {
			 labelLeavesWithWeightedModels(tree,taxaSubsets);
		 } else {
			 labelLeavesWithModels(tree,taxaSubsets);
		 }
	 }
		/**
		 * Apply labels to each taxon saying which evolution models apply to it - in preparation for
		 * a MosaicSubstitutionModel analysis.
		 * 
		 * @param tree
		 * @param taxaSubsets taxaSubset[submodel] is the set of taxa which 'submodel' applies to, where 
		 * taxa are in the order returned by tree.getExternalNode
		 */
	public void labelLeavesWithModels(Tree tree, SetOfSmallIntegers[] taxaSubsets) {
		int nLeaves = tree.getExternalNodeCount();
		int nModels = taxaSubsets.length;
		for (int leaf=0; leaf<nLeaves; leaf++) {
			SetOfSmallIntegers models = new SetOfSmallIntegers();
			for (int model=0; model < nModels; model++) {
				if (taxaSubsets[model].isMember(leaf)) {
					models.addToSet(model);
				}
			}
			AttributeNode node = (AttributeNode)tree.getExternalNode(leaf);
			ModelNodeParameters mnp = new MosaicModelNodeParameters(models);
			node.setAttribute(ModelNodeParameters.ATTRIBUTE_LABEL, mnp);
		}
	}
	
	/**
	 * Sets all leaf nodes to have all models apply to it, and uses taxaSubsets to give inital
	 * weights to the edge leading to this leaf node: if taxaSubsets[submodel] includes this node,
	 * give larger weight to submodel in the parent edge. 
	 * @param tree
	 * @param taxaSubsets
	 */
	public void labelLeavesWithWeightedModels(Tree tree, SetOfSmallIntegers[] taxaSubsets) {
		final double WEIGHT_RATIO = 4; // How much more heavily to weight members of taxaSubsets. 
		int nLeaves = tree.getExternalNodeCount();
		int nModels = taxaSubsets.length;
		SetOfSmallIntegers modelsSet = new SetOfSmallIntegers(nModels);
		int[] modelsList = new int[nModels];
		for (int i=0; i<nModels; i++) {modelsList[i]=i;}
		double[] modelWeights = new double[nModels];
		double[] parameters   = new double[nModels-1];
		for (int leaf=0; leaf<nLeaves; leaf++) {
			AttributeNode node = (AttributeNode)tree.getExternalNode(leaf);
			ModelNodeParameters mnp = new MosaicModelNodeParameters(modelsSet);
			node.setAttribute(ModelNodeParameters.ATTRIBUTE_LABEL, mnp);
			
			MosaicModelEdgeParameters mep = new MosaicModelEdgeParameters(this);
			mep.setNonZeroWeights(modelsList); // all models can have non-zero weight.
			double sum = 0;
			for (int model=0; model<nModels; model++) {
				modelWeights[model] = taxaSubsets[model].isMember(leaf) ? WEIGHT_RATIO : 1;
				sum += modelWeights[model];
			}
			// normalize
			for (int model=0; model<nModels; model++) { modelWeights[model] /= sum; }
			// convert to stochastic vector parameterized form:
			StochasticVector.makeParameters(modelWeights, parameters, 0);
			// set the MMEP edge weights via the parameters
			for (int param=0; param<nModels-1; param++) {
				mep.setParameter(parameters[param],param);
			}
			// and finally attach the edge params as an attribute
			node.setAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL, mep);
		}
	}

	/**
	 * Set up a tree to use this MosaicSubstitutionModel. (Alters the edges to allow them to have edge weights for the 
	 * submodels in this Mosaic model.)
	 * @param tree
	 */
	public void attachToTree(Tree tree) {
		// Attach an MMEP for this model to every node
		// TODO: Exception: only one of root's children should have an MMEP.
		int nLeaves = tree.getExternalNodeCount();
		for (int i=0; i<nLeaves; i++) {
			AttributeNode leaf = (AttributeNode)tree.getExternalNode(i);
			if (leaf.getAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL) == null) {
				MosaicModelEdgeParameters mmep = new MosaicModelEdgeParameters(this);
				leaf.setAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL, mmep);
			}
		}
		int nInternal = tree.getInternalNodeCount();
		for (int i=0; i<nInternal; i++) {
			AttributeNode internal = (AttributeNode)tree.getInternalNode(i);
			if (internal.getAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL) == null) {
				MosaicModelEdgeParameters mmep = new MosaicModelEdgeParameters(this);
				internal.setAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL, mmep);
			}
		}
	}
	
	public ModelEdgeParameters getRootModelEdgeParameters() {
		return new MosaicModelEdgeParameters(this);
	}
	
	/**
	 * Takes a tree which has had labels attached (labelLeavesWithModels). Returns the
	 * subtree left by keeping only the leaves to which a given model applies.
	 * @param masterTree
	 * @param submodel
	 * @return a tree.
	 */
	
	public static Tree treeRestrictedToModel(Tree masterTree, int submodel) {
		// Find which leaves which use 'model'
		int nLeaves = masterTree.getExternalNodeCount();
		Vector<Integer> leafList = new Vector<Integer>();
		for (int leaf = 0; leaf < nLeaves; leaf++) {
			AttributeNode leafNode = (AttributeNode)masterTree.getExternalNode(leaf);
			MosaicModelNodeParameters mmnp = (MosaicModelNodeParameters)leafNode.getAttribute(ModelNodeParameters.ATTRIBUTE_LABEL);
			if (mmnp.getSubmodelsInSubtree(0).isMember(submodel)) {
				// This leaf uses 'model'
				leafList.add(leaf);
			}
		}
		// Feed their names to TreeRestricter and get the restricted tree
		String[] leafNames = new String[leafList.size()];
		for (int i = 0; i<leafNames.length; i++) {
			leafNames[i] = masterTree.getExternalNode(leafList.elementAt(i)).getIdentifier().toString();
		}
		TreeRestricter tr = new TreeRestricter(masterTree, leafNames, true);
		Tree restrictedTree = tr.generateTree();
		// TODO: one day, it might be useful to copy the MosaicModelNodeParameters over to the restricted tree.
		return restrictedTree;
	}
	
	/**
	 * Requires masterTree to have had leaves labeled with labelLeavesWithModels
	 * @param alignment
	 * @param masterTree
	 * @param submodel
	 */
	public void optimiseSubmodel(Alignment alignment, Tree masterTree, int submodel) {
		Tree submodelTree = treeRestrictedToModel(masterTree, submodel);
		TreeSearchTool.optimiseUnrootedFixed(submodelTree, alignment, this.models_[submodel], true);
		submodelChanged_[submodel] = true;
	}

	
	public DataType getDataType() {		return dataType_;		}
	public int getNumberOfTransitionCategories() {return rateDistribution_.getNumberOfRates();}
	public double getTransitionCategoryProbability(int category) {	return rateDistribution_.getCategoryProbability(category); }
	public double[] getTransitionCategoryProbabilities() { return rateDistribution_.getCategoryProbabilities(); }
	public double[] getRateDistributionParameters() { return parameterization_.getBaseParameters(rateDistribution_); }
	
	public void setRootFrequencies(double[] freqs) {
		rootFrequencies_.setVector(freqs);
	}
	
	public double[] getRootFrequencies() {
		return rootFrequencies_.getVector();
	}
	public void getRootFrequencies(double[] store) {
		System.arraycopy(rootFrequencies_.getVector(), 0, store, 0, rootFrequencies_.length());
	}
	public int getNumberSubmodels() { return numberOfSubmodels_;}

	public void getTransitionProbabilities(double branchLength, ModelEdgeParameters params, int category, double[][] tableStore) {
		MosaicModelEdgeParameters meParams;

		updateSubmodelRateMatrices();
		if (params instanceof MosaicModelEdgeParameters) {
			meParams = (MosaicModelEdgeParameters)params;
		} else {
			throw new RuntimeException("I need to be passed a MosaicModelEdgeParameters object");
		}
		meParams.getTransitionProbabilities(branchLength*rateDistribution_.getRate(category), tableStore);
	}
	
	public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters params, int category, double[][] tableStore) {
		MosaicModelEdgeParameters meParams;

		updateSubmodelRateMatrices();
		if (params instanceof MosaicModelEdgeParameters) {
			meParams = (MosaicModelEdgeParameters)params;
		} else {
			throw new RuntimeException("I need to be passed a MosaicModelEdgeParameters object");
		}
		meParams.getTransitionProbabilitiesTranspose(branchLength*rateDistribution_.getRate(category), tableStore);
	}
	
	
	public void getTransitionProbabilities(double branchLength, ModelEdgeParameters params, double[][][] tableStore) {
		int nCategories = rateDistribution_.getNumberOfRates();
		MosaicModelEdgeParameters meParams;
		
		if (params instanceof MosaicModelEdgeParameters) {
			meParams = (MosaicModelEdgeParameters)params;
		} else {
			throw new RuntimeException("I need to be passed a MosaicModelEdgeParameters object");
		}
		
		updateSubmodelRateMatrices();
		for (int category=0; category<nCategories; category++) {
			meParams.getTransitionProbabilities(branchLength*rateDistribution_.getRate(category), tableStore[category]);
		}
	}
	
	public void getTransitionProbabilitiesTranspose(double branchLength, ModelEdgeParameters params, double[][][] tableStore) {
		int nCategories = rateDistribution_.getNumberOfRates();
		MosaicModelEdgeParameters meParams;
		
		if (params instanceof MosaicModelEdgeParameters) {
			meParams = (MosaicModelEdgeParameters)params;
		} else {
			throw new RuntimeException("I need to be passed a MosaicModelEdgeParameters object");
		}
		
		updateSubmodelRateMatrices();
		for (int category=0; category<nCategories; category++) {
			meParams.getTransitionProbabilitiesTranspose(branchLength*rateDistribution_.getRate(category), tableStore[category]);
		}
	}

	// (SingleClass|GeneralRateDistribuiton)SubstitutionModel do this, so I'm just blindly copying it
	public void addPalObjectListener(PalObjectListener l) {
		throw new RuntimeException("Sorry, NeoRateMatrix stuff does not work with old likelihood calculators!");
	}
	public void removePalObjectListener(PalObjectListener l) {
		throw new RuntimeException("Sorry, NeoRateMatrix stuff does not work with old likelihood calculators!");
	}

	public boolean mustOptimiseMultipleParametersOnEdge() { return true; }
	
	// Again, blindly following:
	public OrthogonalHints getOrthogonalHints() { return null; }

	/* It is lower maintenance and safer to let the default methods take care of this:
	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(1); //Version number
		out.writeInt(numberOfSubmodels_);
		out.writeObject(models_);
		out.writeObject(submodelRateMatrices_);
		out.writeObject(dataType_);
		out.writeObject(parameterization_);
		throw new RuntimeException("This routine not maintained and out of date");
	}

	private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
		byte version = in.readByte();
		switch(version) {
			default : {
				this.numberOfSubmodels_=in.readInt();
				this.models_ = (SingleClassSubstitutionModel[])in.readObject();
				this.dataType_ = (DataType)in.readObject();
				this.parameterization_ = (MultiParameterized)in.readObject();
				break;
			}
		}
	}
	*/

	public void report(PrintWriter out) {
		out.printf("Root base frequencies:%s", pal.misc.Utils.toString(rootFrequencies_.getVector()));
		if (fixedRootFrequencies_) {
			out.print("(fixed frequencies)\n");
		} else {
			out.print("(optimised with the model)\n");
		}
		out.printf("Rate distribution:%s\n", rateDistribution_.toString());
		for (int i=0; i<numberOfSubmodels_; i++) {
			out.printf("Submodel %d: %s\n", i, models_[i].toString());
		}
	}
	// Won't work because 'report' not implemented yet:
	public String toString() {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw,true);
		report(pw);
		return "Mosaic Substitution Model:\n"+sw.toString();
	}
	/**
	 * WARNING! Seems to be buggy.
	 * I have an undiagnosed problem where toplogy optimisations on cloned MosaicSubstitutionModel
	 * drive all the parameters to zero at which point the Markov matrix becomes non-a-number'ed.
	 */
	public Object clone() {
		return new MosaicSubstitutionModel(this);
	}
	public SubstitutionModel getCopy() {
		return new MosaicSubstitutionModel(this);
	}
	
	/*
	 * Put test cases here.
	 * 
	 * Currently, only contains test cases for the 'orphan edge' problem: when one or more
	 * edges would not get any submodel from the basic 'if the edge lies on the path
	 * between two edges with a given submodel' algorithm.
	 */
	public static void test() {
		// The trees for the test cases
		Tree tree1a = TreeUtils.stringToTree("(A1,(A2,(BC1,BC2)));");
		Tree tree1b = TreeUtils.stringToTree("((A1,A2),(BC1,BC2));"); // same as tree1a except for rooting
		// This tree has two generation 2 orphan edges (i.e. where all neighbouring edges are also orphans.)
		Tree tree2a = TreeUtils.stringToTree("(A1,(A2,(((B1,B2),(C1,C2)),((D1,D2),(E1,E2)))));");
		// test 2b relies on modifying the tree from 2a.
		// Taxon subsets for tests 1a, 1b
		IdGroup[] taxaSubsetsIDs1 = new IdGroup[3];
		taxaSubsetsIDs1[0] = new SimpleIdGroup(new String[] {"A1","A2"});
		taxaSubsetsIDs1[1] = new SimpleIdGroup(new String[] {"BC1","BC2"});
		taxaSubsetsIDs1[2] = new SimpleIdGroup(new String[] {"BC1","BC2"});
		IdGroup allIDs1 = new SimpleIdGroup(new String[] {"A1","A2","BC1","BC2"});
		// Taxon subsets for tests 2a, 2b:
		IdGroup[] taxaSubsetsIDs2 = new IdGroup[5];
		taxaSubsetsIDs2[0] = new SimpleIdGroup(new String[] {"A1","A2"});
		taxaSubsetsIDs2[1] = new SimpleIdGroup(new String[] {"B1","B2"});
		taxaSubsetsIDs2[2] = new SimpleIdGroup(new String[] {"C1","C2"});
		taxaSubsetsIDs2[3] = new SimpleIdGroup(new String[] {"D1","D2"});
		taxaSubsetsIDs2[4] = new SimpleIdGroup(new String[] {"E1","E2"});
		IdGroup allIDs2 = new SimpleIdGroup(new String[] {"A1","A2","B1","B2","C1","C2","D1","D2","E1","E2"});
		
		// Trivial alignments for tests 1 and 2	
		Alignment alignment1 = new SimpleAlignment(allIDs1, new String[] {"T","T","T","T"}, Nucleotides.DEFAULT_INSTANCE);
		Alignment alignment2 = new SimpleAlignment(allIDs2, new String[] {"T","T","T","T","T","T","T","T","T","T"}, Nucleotides.DEFAULT_INSTANCE);
		
		// Make submodels for use with test cases.
		SingleClassSubstitutionModel[] models1 = new SingleClassSubstitutionModel[3];
		SingleClassSubstitutionModel[] models2 = new SingleClassSubstitutionModel[5];
		double[] flatFreqDist = new double[] {0.25, 0.25, 0.25, 0.25}; 
		for (int i=0; i<5; i++) {
			// Make each submodel a Jukes-Cantor - (it doesn't really matter what they are.)
			NeoRateMatrix myRateMatrix = new GeneralREVRateMatrix(4,new int[] {0,0,0,0,0,0}, new double[] {});
			models2[i] = new SingleClassSubstitutionModel(myRateMatrix, Nucleotides.DEFAULT_INSTANCE, flatFreqDist, true);
			if (i<3) {models1[i] = models2[i];}
		}
		
		//The Mosaic models:
		MosaicSubstitutionModel mosaicModel1 = new MosaicSubstitutionModel(models1,flatFreqDist,true);
		MosaicSubstitutionModel mosaicModel2 = new MosaicSubstitutionModel(models2,flatFreqDist,true);
		
		testSingle("1a", mosaicModel1, taxaSubsetsIDs1, tree1a, alignment1);
		testSingle("1b", mosaicModel1, taxaSubsetsIDs1, tree1b, alignment1);
		PossiblyRootedMLSearcher mls = testSingle("2a", mosaicModel2, taxaSubsetsIDs2, tree2a, alignment2);
		
		mosaicModel2.labelLeavesWithModels(tree2a, taxaSubsetsIDs2, false); 
		mosaicModel2.attachToTree(tree2a);
		
		// Test case 2b: do a topology change on an existing tree and confirm that
		// branch weights still agree with the algorithm.
		// We want the post-topology-change tree to be 'interesting', having orphan-level 2 edges
		// like the 2a tree did. This is 'ensured' simply by having set the rng seed in
		// 'testSingle' to a value which happens to give an appropriate tree. 
		// NOTE!!! That means this test is fragile to changes in the RNG or in the SPR action.
		
		UndoableAction spr = mls.getSPRAction(StoppingCriteria.Utils.getIterationCount( 1 ));
		spr.doAction(-9999, 0);
		mls.calculateLogLikelihood(); 
		Tree optTree = mls.buildPALTree();
		Node root = (SimpleNode) optTree.getRoot();
		StringWriter sw = new StringWriter();
		sw.write("Test case 2b: the tree with submodels is:\n");
		NodeUtils.printNH(new PrintWriter(sw), root, false, false, true,  false, false, 0, false); // no edge lengths, internal labels, line breaks. With edge params.
		sw.write('\n');
		NodeUtils.printNH(new PrintWriter(sw), root, false, false, false, false, false, 0, false);
		System.out.println(sw.toString());
	}
	
	// Perform a single test case
	private static PossiblyRootedMLSearcher testSingle(String testCase, MosaicSubstitutionModel model, IdGroup[] taxaSubsetsIDs, Tree tree, Alignment alignment) {
		
		model.labelLeavesWithModels(tree, taxaSubsetsIDs, false); 
		model.attachToTree(tree);
		// seed of '3' is chosen so that test case 2b has an SPR action which does what we want
		// (I.e. shift a cherry to be next to another cherry.)
		PossiblyRootedMLSearcher mls = new PossiblyRootedMLSearcher(tree,alignment,model,3);
		mls.calculateLogLikelihood(); // we want the side-effects of having done the calculation, ignore the actual the likelihood.
		Tree optTree = mls.buildPALTree();
		SimpleNode root = (SimpleNode) optTree.getRoot();
		StringWriter sw = new StringWriter();
		sw.write("Test case "+testCase+": the tree with submodels is:\n");
		NodeUtils.printNH(new PrintWriter(sw), root, false, false, true, false, false, 0, false); // no edge lengths, internal labels, line breaks. With edge params.
		System.out.println(sw.toString());
		return mls;
	}
}
