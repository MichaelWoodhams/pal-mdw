/**
 * Unofficial modification to PAL by Michael Woodhams (School of IT, University of Sydney)
 * 
 * For use by MosaicModel. In MosaicModel, we have a selection of models. Each edge
 * will have a rate matrix which is a weighted average of a subset of the available models.
 * Each edge differs in which models may be included in the weighted average, and the
 * weights of each model. This object stores the required information.
 */

package pal.substmodel;

import java.io.StringWriter;
import java.text.DecimalFormat;

import pal.math.SetOfSmallIntegers;
import pal.math.StochasticVector;

/**
 * Note important assumption of the MosaicSubstitutionModel: once the submodels have been set, they're not allowed to change.
 * Hence I can get the rate matrices once and not worry about them again.
 * 
 * @author woodhams
 *
 */
/*
 * 'parameters_' and 'weights_' are in a (nearly) 1:1 relationship. 
 * 'weights_' has n dimensions (n=numberNonZeroWeights) and a constraint (sums to one) so it describes
 * a n-1 dimensional space. 'parameters_' is n-1 dimensions, each parameter from 0 to 1 (inclusive) 
 * weights_ is derived from parameters_:
 * weights_[i] = parameters_[i]*(1- sum(weights_[0..i-1]))
 * weights_[n] = (1- sum(weights_[0..n-1])
 * weights_ can only be changed indirectly via setParameter
 * 
 * TODO: The above scheme has since got its own class - StochasticVector. Refactor to use this.
 * 
 * Determining which submodels apply on this edge:
 * The primary rule is that a single submodel applies on a connected subgraph: if taxa on either side 
 * of an edge have a submodel, that submodel applies to the edge. This is calculated via 
 * MosaicModelNodeParameters.getSubmodelsInSubtree(). However, this procedure can result in 
 * 'orphan' edges, which have no submodels. So we apply a second rule: if an edge gets no
 * submodels from the first rule, it finds the nearest edges which do have submodels,
 * and uses those. If there are several equidistant edges, it uses the submodels of all
 * of them.
 * 
 * 'orphanGeneration_' is how far we are from the nearest edge which 'naturally' has
 * submodels. orphanGeneration_=0 for a normal edge. This calculation is done through 
 * the fixOrphanGeneration() and getOrphanGeneration() methods and also uses
 * numberNonZeroWeightsUpToDate_.
 */
public class MosaicModelEdgeParameters implements ModelEdgeParameters {
	private MosaicSubstitutionModel model_; 
	private int numberNonZeroWeights_; // how many submodels apply to this edge
	private boolean numberNonZeroWeightsUpToDate_;
	private int[] submodelNumbers_; // which submodels apply to this edge
	private double[] parameters_;   // see above: mapping of 'weights' to (0-1) range independent variables.
	private double[] parametersSE_;
	private double[] weights_; // weight of each model on this edge - to be optimised.
	private boolean edgeWeightsChanged_;
	private long[] submodelRateMatrixLastVersion_; // to detect when submodels change their rate matrices 
	private double[][][] submodelsRateMatrices_;
	private double[][] rateMatrix_; // TODO: Change to a 'RateMatrix'?
	private int numStates_; // cached from data type
	private MatrixExponential matrixExp_; // TODO: use EigenSystem instead
	private int  leftSubtreeEdgeParametersCount_; // how many free parameters are there for edge weights in the subtree on 'left' of this edge?
	private int rightSubtreeEdgeParametersCount_; // set to -1 to indicate number is currently unknown.
	private MosaicModelNodeParameters leftNodeParameters_, rightNodeParameters_; // neighbouring nodes
	private int orphanGeneration_; // distance (in edges) to nearest edge which has submodels due to the standard submodel application rule.
	
	
	// for interface Markable
	private int[]  markSubmodelNumbers_;
	private double[] markParameters_;
	private double[] markParametersSE_;
	private int  markLeftSubtreeEdgeParametersCount_,  markRightSubtreeEdgeParametersCount_;
	private MosaicModelNodeParameters markLeftNodeParameters_, markRightNodeParameters_;
	
	
	// Constructor for when the models applying to this edge are not yet determined.
	public MosaicModelEdgeParameters(MosaicSubstitutionModel model) {
		model_ = model;
		submodelNumbers_ = null;
		numberNonZeroWeights_ = -1;
		numberNonZeroWeightsUpToDate_ = false;
		weights_ = null;
		submodelsRateMatrices_ = null;
		edgeWeightsChanged_ = true;
		numStates_ = model_.getDataType().getNumStates();
		rateMatrix_ = new double[numStates_][numStates_]; 
		matrixExp_ = new MatrixExponential(numStates_); //TODO: use EigenSystem instead of buggy MatrixExponential 
		parameters_ = null;
		parametersSE_ = null;
		leftSubtreeEdgeParametersCount_ = -1;
		rightSubtreeEdgeParametersCount_ = -1;
		leftNodeParameters_ = null;
		rightNodeParameters_ = null;
		
		markSubmodelNumbers_ = null;
		markParameters_ = null;
		markParametersSE_ = null;
		markLeftSubtreeEdgeParametersCount_ = -1;
		markRightSubtreeEdgeParametersCount_ = -1;
		markLeftNodeParameters_ = null;
		markRightNodeParameters_ = null;
	}
	
	public void setNumberNonZeroWeights(int n) {
		// Nothing to do if we already have n non-zero weights
		// otherwise (re)allocate space.
		if (numberNonZeroWeights_ != n) {
			numberNonZeroWeights_ = n;
			weights_ = new double[numberNonZeroWeights_];
			submodelRateMatrixLastVersion_ = new long[numberNonZeroWeights_];
			submodelsRateMatrices_ = new double[numberNonZeroWeights_][][];
			edgeWeightsChanged_ = true;
			if (numberNonZeroWeights_ >1) {
				// parameters initialized to zero is suitable.
				parameters_ = new double[numberNonZeroWeights_-1];
				parametersSE_ = new double[numberNonZeroWeights_-1]; // I'm not sure these ever get used.
			} else {
				parameters_ = null;
				parametersSE_ = null;
			}
		}
		numberNonZeroWeightsUpToDate_ = true;
	}
	
	/**
	 * Reset weights to be equal
	 */
	public void initializeParameters() {
		for (int i=0; i<numberNonZeroWeights_; i++){
			weights_[i]=1./numberNonZeroWeights_;
		}
		StochasticVector.makeParameters(weights_, parameters_, 0);
	}
	
	public void setNonZeroWeights(int[] nonZeroWeights) {
		setNumberNonZeroWeights(nonZeroWeights.length);
		submodelNumbers_ = nonZeroWeights;
		for (int i=0; i<numberNonZeroWeights_; i++){
			submodelsRateMatrices_[i] = model_.getSubmodelRateMatrix(submodelNumbers_[i]);
			// will always be 'out of date' as MosaicSubstitutionModel increments version to 1 on first evaluation of rate matrix.
			submodelRateMatrixLastVersion_[i] = 0; 
		}
		// By default, start with equal weights
		initializeParameters();
	}
	
	private boolean isSameSubmodelSet(SetOfSmallIntegers comparison) {
		return submodelNumbers_ != null && comparison.equals(new SetOfSmallIntegers(submodelNumbers_));
	}
	
	// Return total number of free model weight parameters over the entire tree.
	// This routine is the point of all of the parameter counting stuff.
	public int getTreeEdgeParametersCount() {
		if (!numberNonZeroWeightsUpToDate_) {
			fixOrphanGeneration(); // this is an orphan. We know have enough info to find its submodels.
		}
		return getLeftSubtreeEdgeParametersCount() + getRightSubtreeEdgeParametersCount() + (numberNonZeroWeights_-1);
	}
	private int getLeftSubtreeEdgeParametersCount() {
		if (leftSubtreeEdgeParametersCount_ == -1) {
			// We have not got a cached value, so need to ask neighbouring node to calculate it for us
			leftSubtreeEdgeParametersCount_ = leftNodeParameters_.getSubtreeEdgeParametersCount(this);
		}
		return leftSubtreeEdgeParametersCount_;
	}
	private int getRightSubtreeEdgeParametersCount() {
		if (rightSubtreeEdgeParametersCount_ == -1) {
			// We have not got a cached value, so need to ask neighbouring node to calculate it for us
			if (rightNodeParameters_ != null) {
				rightSubtreeEdgeParametersCount_ = rightNodeParameters_.getSubtreeEdgeParametersCount(this);
			} else {
				// if this MMEP is attached to Root, we'll execute this
				rightSubtreeEdgeParametersCount_ = 0;
			}
		}
		return rightSubtreeEdgeParametersCount_;
	}

	// A neighbouring MMNP is asking us for the count in out subtree not including that node
	public int getSubtreeEdgeParametersCount(MosaicModelNodeParameters nodeParams) {
		if (!numberNonZeroWeightsUpToDate_) {
			throw new RuntimeException("Bug in Mosaic model edge parameter counting - number not up to date");
		}
		if (nodeParams == leftNodeParameters_) {
			return getRightSubtreeEdgeParametersCount() + numberNonZeroWeights_-1;
		} else if (nodeParams == rightNodeParameters_) {
			return getLeftSubtreeEdgeParametersCount() + numberNonZeroWeights_-1;
		} else {
			throw new RuntimeException("Bug in Mosaic model edge parameter counting - unknown nodeParams");
		}
	}
	
	private void fixOrphanGeneration() {
		SetOfSmallIntegers  leftSubmodelSet = new SetOfSmallIntegers();
		SetOfSmallIntegers rightSubmodelSet = new SetOfSmallIntegers();
		int  leftOrphanGeneration =  leftNodeParameters_.getOrphanGeneration(this,  leftSubmodelSet);
		int rightOrphanGeneration = rightNodeParameters_.getOrphanGeneration(this, rightSubmodelSet);
		if (leftOrphanGeneration < rightOrphanGeneration) {
			// there are non-orphans closer on the left
			orphanGeneration_ = leftOrphanGeneration+1;
			setNonZeroWeights(leftSubmodelSet.getList());	
		} else if (rightOrphanGeneration < leftOrphanGeneration) {
			// there are non-orphans closer on the right
			orphanGeneration_ = rightOrphanGeneration+1;
			setNonZeroWeights(rightSubmodelSet.getList());
		} else {
			// They must be equal
			orphanGeneration_ = leftOrphanGeneration+1;
			leftSubmodelSet.unionEquals(rightSubmodelSet);
			setNonZeroWeights(leftSubmodelSet.getList());	
		}
		numberNonZeroWeightsUpToDate_ = true;
	}
	
	/**
	 * Return what the orphan generation of this edge would be if not for the
	 * connection to 'caller' (i.e. looking in the direction away from 'caller') and
	 * the set of submodels which is closest in the direction away from 'caller'.  
	 * @param caller The node we're looking at this edge from
	 * @param subset Used to return the set of submodels.
	 * @return
	 */
	public int getOrphanGeneration(MosaicModelNodeParameters caller, SetOfSmallIntegers subset) {
		if (numberNonZeroWeightsUpToDate_) {
			// orphanGeneration_ is up to date, so no calculation needed
			subset.set(submodelNumbers_);
			return orphanGeneration_;
		} else {
			MosaicModelNodeParameters other = (caller == leftNodeParameters_) ? rightNodeParameters_ : leftNodeParameters_;
			if (other == null) {
				// special case: this edge is the Root
				return Integer.MAX_VALUE;
			}
			return 1+other.getOrphanGeneration(this,subset);
		}
	}
	
	
	/**
	 * Gets called during the two passes over the 'shadow tree' in UnrootedMaximumLikelihoodSearcher.
	 * leftMNP and rightMNP are the neighbouring ModelNodeParameters - order is not significant.
	 * (In particular, left/right labels have no necessary link to the left/right of UnrootedMaximumLikelihoodSearcher.Connection.)
	 * left/rightBackPointer indicate which edge number we are from the point of view of that node.
	 */
	public void update(
			ModelNodeParameters leftMNP, 
			int leftBackPointer,
			ModelNodeParameters rightMNP,
			int rightBackPointer,
			boolean firstPass) 
	{
		if (firstPass) {
			// Nothing useful to do in first pass yet. Can't even memorize leftNodeParameters_ etc,
			// as neighbouring nodes may not yet have their MNPs.
		} else {
			leftSubtreeEdgeParametersCount_ = -1;
			rightSubtreeEdgeParametersCount_ = -1;
			leftNodeParameters_ = (MosaicModelNodeParameters)leftMNP;
			rightNodeParameters_ = (MosaicModelNodeParameters)rightMNP;
			if (rightNodeParameters_ == null) {
				// special case: this MMEP is attached to 'Root'
				setNonZeroWeights(new int[]{});
				orphanGeneration_ = Integer.MAX_VALUE;
				numberNonZeroWeightsUpToDate_ = true;
			} else {
				// normal case
				SetOfSmallIntegers applicableModels = leftNodeParameters_.getSubmodelsInSubtree(leftBackPointer);
				applicableModels.intersectionEquals(rightNodeParameters_.getSubmodelsInSubtree(rightBackPointer));
				if (!isSameSubmodelSet(applicableModels)) {
					// If the set of models is the same, leave the model weights undisturbed.
					if (applicableModels.isEmpty()) {
						// If no other models apply to this edge, we need to fall back on the 'orphan' method 
						// for allocating submodels - but surrounding edges aren't all yet up to date during
						// 'update', so we can't do it yet. 
						numberNonZeroWeightsUpToDate_ = false;
					} else {
						// Normal case
						setNonZeroWeights(applicableModels.getList());
						orphanGeneration_ = 0; // we are not an orphan
						numberNonZeroWeightsUpToDate_ = true;
					}
				}
			}
		}
	}
	
	
	public double[] getWeights() {
		return weights_.clone();
	}
	
	public MosaicSubstitutionModel getModel() {
		return model_;
	}
	public void setModel(MosaicSubstitutionModel model) {
		model_ = model;
	}
	
	/*public double[][] getRateMatrix() {
		calculateRateMatrix();
		return rateMatrix_;
	}
	*/
	
	private void calculateRateMatrix() {
		if (!numberNonZeroWeightsUpToDate_) {
			fixOrphanGeneration(); // this is an orphan. We know have enough info to find its submodels.
		}
		// Have the rate matricies of any of the submodels changed?
		boolean submodelChanged = false;
		for (int i=0; i<numberNonZeroWeights_; i++) {
			long ver = model_.getSubmodelCurrentVersion(submodelNumbers_[i]);
			if (ver != submodelRateMatrixLastVersion_[i]) {
				submodelChanged = true;
				submodelRateMatrixLastVersion_[i] = ver;
				// submodelsRateMatrices_[] is pointers into MosaicSubstitutionModel.submodelRateMatrices_, so will
				// already be changed - we just need to know it has changed to recalculate rateMatrix_.
			}
		}
		if (!submodelChanged && !edgeWeightsChanged_) {
			// Rate matrix is up to date - nothing to do.
			return;
		}
		int i, j, k;
		for (j=0; j<numStates_; j++) {
			for (k=0; k<numStates_; k++) {
				rateMatrix_[j][k] = weights_[0]*submodelsRateMatrices_[0][j][k];
			}
		}
		for (i=1; i<numberNonZeroWeights_; i++) {
			for (j=0; j<numStates_; j++) {
				for (k=0; k<numStates_; k++) {
					rateMatrix_[j][k] += weights_[i]*submodelsRateMatrices_[i][j][k];
				}
			}
		}
		matrixExp_.updateByRelativeRates(rateMatrix_);
		edgeWeightsChanged_ = false;
	}
	
	public void getTransitionProbabilities(double distance, double[][] store ) {
		calculateRateMatrix();
		matrixExp_.setDistance(distance);
		matrixExp_.getTransitionProbabilities(store);
	}
	public void getTransitionProbabilitiesTranspose(double distance, double[][] store ) {
		calculateRateMatrix();
		matrixExp_.setDistanceTranspose(distance);
		matrixExp_.getTransitionProbabilities(store);
	}

	/**
	 * 'Prints' the edge weights into a string.
	 * Items are comma separated, weights which are not applicable to this edge are printed as a blank.
	 * @param precision The number of decimal places to print
	 * @return
	 */
	public String toString() {
		StringWriter sw = new StringWriter();
		DecimalFormat df = new DecimalFormat("0.000");
		int numModels = model_.getNumberSubmodels();
		double[] allWeights = new double[numModels];
		for (int i=0; i<numModels; i++) {allWeights[i] = -1;}
		for (int i=0; i<numberNonZeroWeights_; i++) {
			allWeights[submodelNumbers_[i]] = weights_[i];
		}
		sw.write("[");
		for (int i=0; i<numModels; i++) {
			if (i>0) { sw.write(","); }
			if (allWeights[i]==-1) {
				sw.write(" ");
			} else {
				sw.write(df.format(allWeights[i]));
			}
		}
		sw.write("]");
		return sw.toString();
	}
	
	// For interface Parameterized
	// TODO: I now have math.StochasticVector class to do this calculation. Use it.
	public int getNumParameters() {
		if(!numberNonZeroWeightsUpToDate_) { fixOrphanGeneration(); }
		return numberNonZeroWeights_-1;	
	}
	public void setParameter(double param, int n) { 
		param = Math.max(0,param);
		param = Math.min(param, 1);
		double sum = 0;
		// Assume parameters_, weights_ were correctly synchronized before this call
		parameters_[n] = param;
		for (int i=0; i<n; i++) { sum += weights_[i]; }
		for (int i=n; i<numberNonZeroWeights_-1; i++) {
			weights_[i]=parameters_[i]*(1-sum);
			sum += weights_[i];
		}
		weights_[numberNonZeroWeights_-1] = 1-sum;
		edgeWeightsChanged_ = true; 
	}
	public double getParameter(int n) { return parameters_[n]; }
	public void setParameterSE(double paramSE, int n) { this.parametersSE_[n] = paramSE; }
	public double getLowerLimit(int n) { return 0; }
	public double getUpperLimit(int n) { return 1; }
	public double getDefaultValue(int n) { return 1/(numberNonZeroWeights_-n); } // gives equal weights_[]
	
	public void mark() {
		// TODO: Make more efficient: if 'mark' matrix already exists and is right size, and is not just a pointer to the actual 
		// (non-mark) matrix, then reuse it.
		markSubmodelNumbers_ = submodelNumbers_.clone();
		if (parameters_ != null) {
			markParameters_   = parameters_.clone();
			markParametersSE_ = parametersSE_.clone();
		} else {
			markParameters_   = null;
			markParametersSE_ = null;
		}
		markLeftSubtreeEdgeParametersCount_ = leftSubtreeEdgeParametersCount_;
		markRightSubtreeEdgeParametersCount_ = rightSubtreeEdgeParametersCount_;
		markLeftNodeParameters_ = leftNodeParameters_;
		markRightNodeParameters_ = rightNodeParameters_;
	}
	public void undoToMark() {
		setNonZeroWeights(markSubmodelNumbers_);
		for (int i=0; i<getNumParameters(); i++) {
			setParameter(markParameters_[i],i);
			setParameterSE(markParametersSE_[i],i);
		}
		leftSubtreeEdgeParametersCount_ = markLeftSubtreeEdgeParametersCount_;
		rightSubtreeEdgeParametersCount_ = markRightSubtreeEdgeParametersCount_;
		leftNodeParameters_ = markLeftNodeParameters_;
		rightNodeParameters_ = markRightNodeParameters_;
	}
}
