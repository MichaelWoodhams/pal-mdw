// PossiblyRootedMLSearcher.java
//
// (c) 1999-2003 PAL Development Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package pal.treesearch;

/**
 * <p>Title: PossiblyRootedMLSearcher</p>
 * <p>Description:
 * A tool for searching for a maximum likelihood unrooted tree under a general substitution model (which may be optimised).
 * Much of the code is delegated to other classes (such as LHCalculator classes) which in turn manage datatype related optimisations, high accuracy  computation, low memory computation, and cached calculation.
 * Even given the offsourcing of code this class is rather large!
 * Includes the algorithm of [1]
 * </p>
 * [1]  Guindon, Stephane  Gascuel, Olivier (2003) A Simple, Fast, and Accurate Algorithm to Estimate Large Phylogenies by Maximum Likelihood. Systematic Biology 52:5 pages  696 - 704 / October 2003
 * @author Matthew Goode
 * @version 1.0
 */

/*
 * HOW IT ALL WORKS:
 * 
 * (Extra note added Dec 2009: most of the other classes in the treesearch package look to be a major but incomplete tidy up of this file.
 * Best long term plan would be to make this file redundant, but it still has significant functionality not yet in the new files.)
 * 
 * These notes added by Michael Woodhams, Feb 2009. This is my understanding of how it works, but I didn't write the program
 * so they are not authoritative.
 * I've added comments within the code. All such comments are prefaced "MDW" and are similarly non-authoritative. 
 *
 * (Also: Anything to do with 'ModelEdgeParameters' or 'ModelNodeParameters' objects is stuff I've added.)
 * 
 * The "end user" should generally be using methods from TreeSearchTool, rather than doing anything with this class directly.
 * This is also a good place to look to see how the class is used.
 * 
 * We start with a tree to be analysed, as a structure of Nodes. This gets copied into a 'shadow tree' (my term) which keeps the
 * same structure, but holds much more information. The shadow tree consists of Connections (edges) and UNodes (nodes). A UNode is
 * either an InternalNode or a LeafNode. Only
 * binary trees are supported. (We could fake a multifurcating tree by constraining some edges to be zero length, but I haven't
 * seen any code to support this.) The structure is doubly linked throughout. 
 *   
 * The internal nodes hold pattern information for each (of 3) adjoining Connections - being all the sequence patterns exhibited by 
 * the subtree at the far end of that Connection. The node also can hold partial likelihoods for each adjacent subtree. (These
 * get recalculated only when necessary.) While profligate on storage, it has efficiency advantages: when any changes are made to the
 * tree (branch length changes, prune and regraft...) any partial likelihood results which don't need recalculating are not
 * recalculated.  
 * 
 * Actual calculations are done with the aid of an LHCalculator object. We can provide LHCalculators which optimise for the common
 * cases (esp. four states (DNA/RNA.) 
 * 
 * The important classes:
 * PossiblyRootedMLSearcher: constructs the shadow tree, calculates likelihoods, allows various optimisation searches. Contains important inner classes:
 *     Connection: (described above)
 *     UNode: (described above)
 *     many UndoableAction subclasses, for manipulations to be optimised over. Work with the SearchEngine class.
 *     ConstructionTool: aids in building the shadow tree
 *     OptimisationHandler: exposes a single parameter (e.g. edge length) and hides the rest, so that it can be optimised
 *         by a standard 1D optimisation method.
 * PatternInfo: site patterns within a subtree
 * LHCalculator: (described above). Each UNode has its own LHCalculator which contains a ConditionalProbabilityStore.
 * ConditionalProbabilityStore: conditional probabilities for a subtree for each site pattern 
 * SubstitutionModel: e.g Jukes Cantor, HKY+I+Gamma etc. 
 * AlgorithmCallback: Has something to do with monitoring progress, aborting calculations which are too slow etc.
 * 
 * 'flat' and 'extended' conditional probabilities:
 * 'flat' probabilities have not been evolved over the length of an edge, 'extended' ones have been.
 * E.g. consider a branch with a 'top' end (closest to root) and a 'bottom' end. The node at the
 * bottom has two 'child' edges descending from it. For this node, we get the 'flat' probabilities
 * by taking the (extended) conditional probability vectors from its two child edges and component-wise
 * multiplying the two vectors. To get the 'extended' conditional probability, we multiply the 'flat'
 * vector by an appropriate (derived from the model and edge length) Markov matrix.
 * 
 */
import pal.tree.*;
import pal.datatype.*;
import pal.math.*;
import pal.util.*;

import pal.substmodel.*;
import pal.misc.*;
import pal.alignment.*;

import pal.algorithmics.*;
import pal.eval.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.util.*;

/*
 * MDW:
 * Class PossiblyRootedMLSearcher.
 * Important data members:
 * + A 'shadow tree' (via treeAccess_)
 * + A substitution model
 * Data members primarily used for during shadow tree construction:
 * + ConstructionTool
 * + calcFactory_. This 'knows' whether the model is 4 state or not and provides optimised LHCalculators.
 * Important public methods:
 * + calculateLogLikelihood (for interface Markable)
 * + simpleOptimiseLikelihood (for interface Markable) optimises edge lengths
 * + buildPALTree (converts shadow tree back to standard PAL Node tree.)
 */
public class PossiblyRootedMLSearcher implements Markable, StateProvider, UnrootedTreeInterface.Instructee {
	
	public static boolean DEBUG = false; // public so can turn on debugging from external program. 
	//The initial value of the update id (when not updated)
	private static final int STARTING_UPDATE_ID = 0;
	private static final double CONSTRUCTED_BRANCH_LENGTH = 0.01;

	private static final double MINIMUM_BRANCH_LENGTH = 0;
	private static final double MAXIMUM_BRANCH_LENGTH = 10;

	private ConnectionOrRoot treeAccess_; // MDW: the 'shadow tree'. As this is doubly linked, we just need a pointer to an arbitrary edge or node within it.
	private final SubstitutionModel model_;
	private final ConstructionTool tool_;
	private final Connection[] allConnections_; // MDW: allConnections_[0] = treeAccess_. 
	private final UNode[] orderedNodes_;

	private final MersenneTwisterFast random_; // MDW:  = new MersenneTwisterFast();

	private final OptimisationHandler optimisationHandler_;

	private final LHCalculator.Factory calcFactory_;
	/**
	 * Build an unconstrained optimiser based on a randomly generated tree.
	 * @param alignment the alignment used to represent each OTU
	 * @param model the substitution model that is used for calcuation. If optimisation on the model occurs than
	 * this model will be altered
	 */
	public PossiblyRootedMLSearcher(Alignment alignment, SubstitutionModel model) {
		this(alignment,model,SimpleModelFastFourStateLHCalculator.getFactory(),0);
	}
	
	public PossiblyRootedMLSearcher(Alignment alignment, SubstitutionModel model, LHCalculator.Factory calcFactory) {
		this(alignment, model, calcFactory, 0);
	}

	/*
	 * MDW: Never use a pseudo-random number generator without at least the option to set the seed.
	 * If nothing else, this allows you to reliably replicate bugs.
	 * Use seed = 0 to get a seed based off the system time.
	 */
	public PossiblyRootedMLSearcher(Alignment alignment, SubstitutionModel model, LHCalculator.Factory calcFactory, long seed) {
		if (seed != 0) {
			random_ = new MersenneTwisterFast(seed);
		} else {
			random_ = new MersenneTwisterFast();
		}
		this.model_ = model;
		int numberOfStates = model_.getDataType().getNumStates();
		this.calcFactory_ = calcFactory;
		tool_ = new ConstructionTool(alignment,numberOfStates,model.getNumberOfTransitionCategories(),calcFactory);
		if (model.isTimeReversible()) {
			this.treeAccess_ = new Connection(Identifier.getNames(alignment),tool_,random_);
		} else {
			this.treeAccess_ = new Root(Identifier.getNames(alignment),tool_,(NonTimeReversibleSubstitutionModel)model,random_);
		}
		ArrayList<Connection> v = new ArrayList<Connection>();
		this.treeAccess_.getAllConnections(v);

		this.allConnections_ = new Connection[v.size()];
		v.toArray(allConnections_);
		if (model_.mustOptimiseMultipleParametersOnEdge()) {
			this.optimisationHandler_ = new MultiOptimisationHandler(model_,tool_);
		} else {
			this.optimisationHandler_ = new UniOptimisationHandler(model_,tool_);
		}
		this.treeAccess_.setup(tool_,allConnections_);
		optimisationHandler_.setup(allConnections_[0],true); // allConnections_[0] == treeAccess_ except when treeAccess_ is a Root.

		this.orderedNodes_ = tool_.getOrderedNodes();
	}

	// MDW: The two most useful, 'user friendly' constructors
	public PossiblyRootedMLSearcher(Tree t, Alignment alignment, SubstitutionModel model) {
		this(t.getRoot(),alignment,model,SimpleModelFastFourStateLHCalculator.getFactory(),0);
	}
	public PossiblyRootedMLSearcher(Node root, Alignment alignment, SubstitutionModel model) {
		this(root,alignment,model,SimpleModelFastFourStateLHCalculator.getFactory(),0);
	}
	// Seeded versions of the above.
	public PossiblyRootedMLSearcher(Tree t, Alignment alignment, SubstitutionModel model, long seed) {
		this(t.getRoot(),alignment,model,SimpleModelFastFourStateLHCalculator.getFactory(),seed);
	}
	public PossiblyRootedMLSearcher(Node root, Alignment alignment, SubstitutionModel model, long seed) {
		this(root,alignment,model,SimpleModelFastFourStateLHCalculator.getFactory(),seed);
	}

	/**
	 * Create a searcher based on a given tree, that has no alignment specified (useful as backbone tree for attaching new nodes)
	 * @param root the root of the tree to base things on (doesn't matter if it's rooted)
	 * @param model the substitution model to be used
	 */
	public PossiblyRootedMLSearcher(Node root,  SubstitutionModel model) {
		this(root,null,model,SimpleModelFastFourStateLHCalculator.getFactory(),0);
	}
	/**
	 * Create a searcher based on a given tree, that has no alignment , or model, specified (useful as backbone tree for attaching new nodes)
	 * @param root the root of the tree to base things on (doesn't matter if it's rooted)
	 */
	public PossiblyRootedMLSearcher(Node root) {
		this(root,null,null,SimpleModelFastFourStateLHCalculator.getFactory(),0);
	}
	
	public PossiblyRootedMLSearcher(Node root, Alignment alignment, SubstitutionModel model, LHCalculator.Factory calcFactory) {
		this(root, alignment, model, calcFactory, 0);
	}
	// MDW: Most constructors are based on this one.
	// seed = 0 uses the system time to seed the random number generator.
	// I'm not sure that random_ is actually used in this case anyhow - no random topology changes.
	public PossiblyRootedMLSearcher(Node root, Alignment alignment, SubstitutionModel model, LHCalculator.Factory calcFactory, long seed) {
		if (seed != 0) {
			random_ = new MersenneTwisterFast(seed);
		} else {
			random_ = new MersenneTwisterFast();
		}
		this.calcFactory_ = calcFactory;
		this.model_ = model;
		if(model_==null) {
		  tool_ = new ConstructionTool(alignment,0,1,calcFactory);
		} else {
		  int numberOfStates = model_.getDataType().getNumStates();
		  tool_ = new ConstructionTool(alignment,numberOfStates,model.getNumberOfTransitionCategories(),calcFactory);
		}
		// MDW: Recursively creates the 'shadow tree' of UNodes and Connections to copy the structure of the tree under 'root'.
		// Allocate objects for patterns, partial probabilities. For the leaf nodes, set up partial probabilities 
		// as pointers into the corresponding edge's Markov matrix (but Markov matrices not yet calculated.)
		if (model.isTimeReversible()) {
			this.treeAccess_ = new Connection(root,tool_);
		} else {
			this.treeAccess_ = new Root(root,tool_,(NonTimeReversibleSubstitutionModel)model);
		}
		ArrayList<Connection> v = new ArrayList<Connection>();
		this.treeAccess_.getAllConnections(v);

		this.allConnections_ = new Connection[v.size()];
		v.toArray(allConnections_);
		// Recursively set up all of the patterns within the shadow tree
		this.treeAccess_.setup(tool_,allConnections_);
		if(model_==null) {
		  this.optimisationHandler_ = null;
		} else {
		  if (model_.mustOptimiseMultipleParametersOnEdge()) {
			this.optimisationHandler_ = new MultiOptimisationHandler(model_,tool_);
		  } else {
			this.optimisationHandler_ = new UniOptimisationHandler(model_,tool_);
		  }
		  //  Calculate partial probabilities through the tree
		  optimisationHandler_.setup(allConnections_[0],true); // allConnections_[0] == treeAccess_ except when treeAccess_ is a Root.
		}
		this.orderedNodes_ = tool_.getOrderedNodes();
	}
	private PossiblyRootedMLSearcher(PossiblyRootedMLSearcher base, Connection attachmentPoint, Node newSubtree, Alignment newSequences, SubstitutionModel model) {
		this(base, attachmentPoint, newSubtree, newSequences, model, 0);
	}
	private PossiblyRootedMLSearcher(PossiblyRootedMLSearcher base, Connection attachmentPoint, Node newSubtree, Alignment newSequences, SubstitutionModel model, long seed) {
		if (seed != 0) {
			random_ = new MersenneTwisterFast(seed);
		} else {
			random_ = new MersenneTwisterFast();
		}
		this.model_ = model;
		this.calcFactory_  = base.calcFactory_;
		int numberOfStates = model_.getDataType().getNumStates();
		tool_ = new ConstructionTool(newSequences,numberOfStates,model_.getNumberOfTransitionCategories(),calcFactory_);
		// MDW: being lazy and not implementing code I won't use. I have not made the modifications
		// to bring this method up to date with the fact that treeAccess_ is now a ConnectionOrRoot
		// instead of a Connection.
		if (!(base.treeAccess_ instanceof Connection)) {
			throw new RuntimeException("Not implemented error: This method not updated Code for rooted tree");
		}
		this.treeAccess_ = new Connection((Connection)base.treeAccess_, attachmentPoint, newSubtree,tool_);
		ArrayList<Connection> v = new ArrayList<Connection>();
		this.treeAccess_.getAllConnections(v);

		this.allConnections_ = new Connection[v.size()];
		v.toArray(allConnections_);
		if (model_.mustOptimiseMultipleParametersOnEdge()) {
			this.optimisationHandler_ = new MultiOptimisationHandler(model_,tool_);
		} else {
			this.optimisationHandler_ = new UniOptimisationHandler(model_,tool_);
		}
		this.treeAccess_.setup(tool_,allConnections_);
		optimisationHandler_.setup(allConnections_[0],true); // allConnections_[0] == treeAccess_ except when treeAccess_ is a Root.

		this.orderedNodes_ = tool_.getOrderedNodes();
	}
	public BranchAccess[] getAccessToBranches() {
	  BranchAccess[] bas = new BranchAccess[allConnections_.length];
		for(int i = 0 ; i < bas.length ; i++) {
		  bas[i] = new BranchAccessImpl(allConnections_[i],this);
		}
		return bas;
	}
	public NodeAccess[] getAccessToNodes() {
	  NodeAccess[] bas = new NodeAccess[orderedNodes_.length];
		for(int i = 0 ; i < bas.length ; i++) {
		  bas[i] = new NodeAccessImpl(orderedNodes_[i],this);
		}
		return bas;
	}
	// Returns a string summarizing all of the parameters relevant to likelihood calculation
	// (including tree topology). Intended for use by debugger.
	public String toString() {
		return unrootedMLSearcherToString(model_,allConnections_) + "\n" + buildPALTree().toString();
	}
	// static method for when we have model and connections but no access to an PossiblyRootedMLSearcher object
	// Unlike the non-static method, does not output the tree. (Could be done, but messy.)
	private static String unrootedMLSearcherToString(SubstitutionModel model, Connection[] allConnections) {
		StringBuffer buffer = new StringBuffer("Model params (");
		for (int i=0; i<model.getNumParameters(); i++) {
			if (i!=0) { buffer.append(", "); }
			buffer.append(Double.toString(model.getParameter(i)));
		}
		buffer.append("), branches (");
		for (Connection c: allConnections) {
			buffer.append(c.debugString()).append(", ");
		}
		buffer.append(")");
		return buffer.toString();
	}
// -=-=-=-=-=-=-==-=--=--==--=-=-=-=-==-=--==-=
// -=-=-=-= State Provider stuff -=-==-=-=-=-=-

// Implement correctly! 
// MDW: Thanks for that useful comment. Any hints on how?

	public Object getStateReference() {
		return new StateObject(allConnections_);
	}

	public void restoreState(Object stateReference) {
		StateObject so = (StateObject)stateReference;
		so.rebuildTree(allConnections_,orderedNodes_);
	}

// -=-=-=-=-=-=-==--=-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-==-=-=-=-=-=--==--

	public void instruct(UnrootedTreeInterface treeInterface) {
	  UnrootedTreeInterface.BaseBranch base = treeInterface.createBase();
		treeAccess_.instructBase(base);
	}


	public UndoableAction getNNIAction(StoppingCriteria.Factory stopper) {
		Assessor a = getSimpleAssessor(stopper);
		return new NNIAction(allConnections_,a,random_,tool_);
	}
	public UndoableAction getSweepNNIAction(StoppingCriteria.Factory stopper) {
		Assessor a = getSimpleAssessor(stopper);
		NNIAction base = new NNIAction(allConnections_,a,random_,tool_);
		return new SweepNNIAction(allConnections_,this,base,random_,tool_);
	}
	public UndoableAction getFullSweepNNIAction(StoppingCriteria.Factory stopper) {
		Assessor a = getSimpleAssessor(stopper);
		NNIAction base = new NNIAction(allConnections_,a,random_,tool_);
		return new FullSweepNNIAction(allConnections_,this,base,random_,tool_);
	}
	
	public UndoableAction getBranchLengthOptimiseAction(StoppingCriteria.Factory stopper) {
		return new BranchLengthOptimiseAction(allConnections_,stopper.newInstance(),AlgorithmCallback.Utils.getNullCallback());
	}
		/**
		 *
		 * @param stopper The means for determining when a set of round should be stopped
		 * @return An undoable action that does the Simulataneous NNI/Branch length of Stephan Guindon
		 * @note this action cannot undo (well, it could but it hasn't been implemented). This is okay as it should always find a better or equal valued state
		 */
	public UndoableAction getNNIBranchLengthOptimiseAction(StoppingCriteria.Factory stopper) {
		return new NNIBranchLengthOptimiseAction(allConnections_,model_, stopper.newInstance(),AlgorithmCallback.Utils.getNullCallback(),tool_);
	}
	public UndoableAction getBranchLengthWithModelOptimiseAction(StoppingCriteria.Factory stopper, MultivariateMinimum minimiser, int fxFracDigits, int xFracDigits, MersenneTwisterFast random) {
		return new BranchLengthWithModelOptimiseAction(tool_, allConnections_,stopper.newInstance(),AlgorithmCallback.Utils.getNullCallback(), minimiser, MinimiserMonitor.Utils.createNullMonitor(),model_,fxFracDigits,xFracDigits,random);
	}
	public UndoableAction getModelOptimiseAction(MultivariateMinimum minimiser, int fxFracDigits, int xFracDigits) {
		return new ModelOptimiseAction(treeAccess_,minimiser, MinimiserMonitor.Utils.createNullMonitor(),model_,fxFracDigits,xFracDigits, tool_);
	}
	public UndoableAction getModelOptimiseAction(MultivariateMinimum minimiser, MinimiserMonitor monitor, int fxFracDigits, int xFracDigits) {
		return new ModelOptimiseAction(treeAccess_,minimiser, monitor,model_,fxFracDigits,xFracDigits, tool_);
	}



	public UndoableAction getSPRAction(StoppingCriteria.Factory stopper) {
		Assessor a = getSimpleAssessor(stopper);
		return new SPRAction(allConnections_,this, a,random_,tool_);
	}
	public UndoableAction getSweepSPRAction(StoppingCriteria.Factory stopper) {
		Assessor a = getSimpleAssessor(stopper);
		SPRAction base = new SPRAction(allConnections_,this, a,random_,tool_);
		return new SweepSPRAction(allConnections_,base,random_);
	}
	public UndoableAction getFullSweepSPRAction(StoppingCriteria.Factory stopper) {
		Assessor a = getSimpleAssessor( stopper );
		SPRAction base = new SPRAction( allConnections_, this, a, random_, tool_ );
		return new FullSweepSPRAction( allConnections_, base );
}



//=--=-=-=-=-=-=-==--==--=
// For Markable interface
	public final void mark() {
		for(int i = 0 ; i < allConnections_.length ; i++) {
			allConnections_[i].mark();
		}
		for(int i = 0 ; i < orderedNodes_.length ; i++) {
			orderedNodes_[i].mark();
		}
	}
	public final void undoToMark() {
		for(int i = 0 ; i < allConnections_.length ; i++) {
			allConnections_[i].undoToMark();
		}
		for(int i = 0 ; i < orderedNodes_.length ; i++) {
			orderedNodes_[i].undoToMark();
		}
		treeAccess_.setup(tool_,allConnections_);
	}

	private final Connection getRandomConnection() {
		return allConnections_[random_.nextInt(allConnections_.length)];
	}
	/**
	 * Recurse through the tree and evaluate likelihoods on all edges in all possible ways.
	 * For testing only.
	 */
	public void testLikelihood() {
		Workspace workspace = new Workspace(30,tool_.getNumberOfSites(),tool_);
		OptimisationHandler oh;
		if (model_.mustOptimiseMultipleParametersOnEdge()) {
			oh = new MultiOptimisationHandler(model_,tool_);
		} else {
			oh = new UniOptimisationHandler(model_,tool_);
		}
		oh.setup(allConnections_[0],true); // allConnections_[0] == treeAccess_ except when treeAccess_ is a Root.
		treeAccess_.testLikelihood(model_,tool_);
	}
	 /**
		* Likelihood calculation method (not optimisation)
		* @return the log likelihood, based on current model, branchlengths and topology
		* @note not valid if no alignment/model given in construction
		*/
	public double calculateLogLikelihood() {
		 return treeAccess_.calculateLogLikelihood(model_,true, tool_.allocateNewExternalCalculator(), tool_);
	}
	 /**
		* An alternative likelihood calculation method (should give same results as other method, and in same time)
		* @return the log likelihood, based on current model, branchlengths and topology
		* @note not valid if no alignment/model given in construction
		*/
	public double calculateLogLikelihood2() {
		return treeAccess_.calculateLogLikelihood2(model_,true, tool_.allocateNewExternalCalculator(), tool_);
	}

	public SiteDetails calculateSiteDetails() {
		 return treeAccess_.calculateSiteDetails(model_,true, tool_.allocateNewExternalCalculator(), tool_);
	}

 /**
		* Optimise the branch lengths of the tree to obtain the maximum likelihood. Does not change the model
		* or the topology
		* @param epsilon the tolerance places for convergence (on the likelihood score)
		* @param callback a callback to monitor progress
		* @return the resulting likelihood
		*/
	public double simpleOptimiseLikelihood(double epsilon, AlgorithmCallback callback) {
	  return simpleOptimiseLikelihood(StoppingCriteria.Utils.getNonExactUnchangedScore(2,true,epsilon).newInstance(),callback);
	}

	 /**
		* Optimise the branch lengths of the tree to obtain the maximum likelihood. Does not change the model
		* or the topology
		* @param stopper the stopping criteria (on the likelihood score)
		* @param callback a callback to monitor progress
		* @return the resulting likelihood
		*/
	public double simpleOptimiseLikelihood(StoppingCriteria stopper, AlgorithmCallback callback) {
		stopper.reset();
		double maximum = Double.NEGATIVE_INFINITY;
		boolean firstTime = true;
		while(!stopper.isTimeToStop()) {
			for (Connection c : allConnections_) {
				maximum = optimisationHandler_.optimiseBranchLength(c,firstTime);
				firstTime = false;
			}
			stopper.newIteration(maximum,maximum,true,true,callback);
		}
		return maximum;
	}


	public Tree buildPALTree() {
		return new SimpleTree(buildPALNode());
	}
	public Node buildPALNode() {
		Node n = treeAccess_.buildPALNode();
		NodeUtils.lengths2Heights(n);
		return n;
	}
// =--=-=-=-=-==--=-=-=-=-=-=-=-=-=-=-=-=-=-==--==--=-=-=-=-=-==--=-=-=-=-==--=
// === State stuff
// =--=-=-=-=-==--=-=-=-=-=-=-=-=-=-=-=-=-=-==--==--=-=-=-=-=-==--=-=-=-=-==--=
	public static final class StateObject {
		private final double[] branchLengths_;
		private final int[] connectionInfo_;
		public StateObject(Connection[] allConnections) {
			this.branchLengths_ = new double[allConnections.length];
			this.connectionInfo_ = new int[allConnections.length*2];
			for(int i = 0 ; i < allConnections.length ; i++) {
				this.branchLengths_[i] = allConnections[i].getBranchLength();
				allConnections[i].fillInConnectionState(connectionInfo_,i*2);
			}
		}
		private final int fillIn(Connection[] store, int nodeIndex, Connection[] allConnections) {
			int found = 0;

			for(int i = 0 ; i < allConnections.length ; i++ ) {
				if(connectionInfo_[i*2]==nodeIndex||connectionInfo_[i*2+1]==nodeIndex) {
					store[found++]=allConnections[i];
				}
			}
			return found;
		}
		public final void rebuildTree(Connection[] allConnections, UNode[] orderedNodes) {
			// MDW: Is this ever used? If so, I need to also remember ModelEdgeParameters, if any.
			if (true) throw new RuntimeException("Need to implement remembering ModelEdgeParameters in UMLS.ObjectState"); 
			for(int i = 0 ; i < allConnections.length ; i++) {
				allConnections[i].setNodes(orderedNodes[connectionInfo_[i*2]],orderedNodes[connectionInfo_[i*2+1]]);
			}
			Connection[] store = new Connection[3];
			for(int nodeIndex = 0 ; nodeIndex < orderedNodes.length ; nodeIndex++) {
				int found = fillIn(store,nodeIndex, allConnections);
				orderedNodes[nodeIndex].setConnections(store,found);

			}
		}
		public boolean equals(Object o) {
			if(o instanceof StateObject) {
				StateObject so = (StateObject)o;
				int[] other = so.connectionInfo_;
				for(int i = 0 ; i < connectionInfo_.length ; i++) {
					if(connectionInfo_[i]!=other[i]) {
						return false;
					}
				}
				return true;
			}
			return false;
		}
	}
// =--=-=-=-=-==--=-=-=-=-=-=-=-=-=-=-=-=-=-==--==--=-=-=-=-=-==--=-=-=-=-==--=
// === Algorithmics stuff Stuff
// =--=-=-=-=-==--=-=-=-=-=-=-=-=-=-=-=-=-=-==--==--=-=-=-=-=-==--=-=-=-=-==--=
	public final Assessor getSimpleAssessor(StoppingCriteria.Factory stopper) {
		return new SimpleAssessor(stopper.newInstance(),AlgorithmCallback.Utils.getNullCallback());
	}

	//-=-=-==--==-

	private final class SimpleAssessor implements Assessor {
		private final StoppingCriteria stopper_;
		private final AlgorithmCallback callback_;
		public SimpleAssessor(StoppingCriteria stopper, AlgorithmCallback callback) {
			this.stopper_ = stopper;
			stopper_.reset();
			this.callback_ = callback;
		}
		public double getCurrentValue() {
			return simpleOptimiseLikelihood(stopper_, callback_);
		}
	}

// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
		/**
		 *
		 * <p>Title: NNIBranchLengthOptimiseAction</p>
		 * <p>Description: Implements the simultaneous NNI/branch length optimisation of Stephan Guindon et al</p>
		 */
	/*
	 * MDW: This action seems little used, so I haven't updated it to work with non-time-reversible models (rooted trees) 
	 * or extra edge parameters. 
	 */
	private static final class NNIBranchLengthOptimiseAction implements UndoableAction {
		private static final String DESCRIPTION_ = "NNIBranchLengthOptimise";
		private final StoppingCriteria stopper_;
		private final AlgorithmCallback callback_;
//    private final double[] branchLengths_;
		private final Connection[] allConnections_;
		private NNIOptimisationHandler handler_;
		public NNIBranchLengthOptimiseAction(Connection[] allConnections, SubstitutionModel model,  StoppingCriteria stopper, AlgorithmCallback callback, ConstructionTool tool) {
			if (model.mustOptimiseMultipleParametersOnEdge() || model.isTimeReversible()) {
				throw new RuntimeException("NNIBranchLengthOptimiseAction doesn't work with rooted trees or extra edge parameters");
				// It could be made to do so, but it is tricky. Could subclass NNIOptimisationHandler into 
				// single and multi parameter versions (like OptimisationHandler), but need to be careful in 'multi'
				// version that node and edge parameters get rebuilt correctly. Could do this by using SPR actions
				// choosing the edges to SPR on so that the result is an NNI.
			}
			this.stopper_ = stopper;
			stopper_.reset();
			this.callback_ = callback;
			this.allConnections_ = allConnections;
			this.handler_ = new NNIOptimisationHandler(allConnections, model, tool);
//      this.branchLengths_ = new double[allConnections.length];
		}
		public double doAction(double currentScore, double desparationValue) {
			if (DEBUG) System.out.println("+++NNIBranchLengthOptimise"); Date start = new Date();
//      for(int i = 0 ; i < allConnections_.length ; i++) {
//        branchLengths_[i] = allConnections_[i].getBranchLength();
//      }
			stopper_.reset();
			double maximum = Double.NEGATIVE_INFINITY;
			boolean firstTime = true;
			while(!stopper_.isTimeToStop()) {
				for(int i = 0 ; i < allConnections_.length ; i++) {
					Connection c = allConnections_[i];
					//Do unchanged way
					maximum = handler_.optimiseSimulataneousNNIBranchLength(c,firstTime);
					firstTime = false;
					}
				stopper_.newIteration(maximum,maximum,true,true,callback_);
			}
			if (DEBUG) System.out.printf("---NNIBranchLengthOptimise %d ms; %.4f\n",(new Date()).getTime()-start.getTime(),maximum); 
			return maximum;
		}
		public boolean isActionSuccessful() {
			return true;
		}
		/**
		 *
		 * @return false
		 */
		public boolean undoAction() {
//      for(int i = 0 ; i < allConnections_.length ; i++) {
//        allConnections_[i].setBranchLength(branchLengths_[i]);
//      }
			return false;
		}
		/**
		 * @return false
		 */
		public boolean isActionDeterministic() {
			return false;
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	}

	// Two helper routines for BranchLengthOptimiseAction and BranchLengthWithModelOptimiseAction
	private static void storeBranchData(Connection[] allConnections, double[] branchLengths, double[][] edgeParameters) {
		for(int i = 0 ; i < allConnections.length ; i++) {
			Connection c = allConnections[i];
			branchLengths[i] = c.getBranchLength();
			if (c.edgeParams_ == null || c.edgeParams_.getNumParameters() == 0) {
				edgeParameters[i] = null;
			} else {
				int n = c.edgeParams_.getNumParameters();
				// test is for efficiency to avoid needless object creation.
				if (edgeParameters[i] == null || edgeParameters[i].length != n) {
					edgeParameters[i] = new double[n];
				}
				for (int j=0; j<n; j++) {
					edgeParameters[i][j] = c.edgeParams_.getParameter(j);
				}
			}
		}
	}
	private static void restoreBranchData(Connection[] allConnections, double[] branchLengths, double[][] edgeParameters) {
		for(int i = 0 ; i < allConnections.length ; i++) {
			allConnections[i].setBranchLength(branchLengths[i]);
			if (edgeParameters[i] != null) {
				for (int j=0; j<edgeParameters[i].length; j++) {
					allConnections[i].edgeParams_.setParameter(edgeParameters[i][j], j);
				}
			}
		}
	}
// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
	private final  class BranchLengthOptimiseAction implements UndoableAction {
		private static final String DESCRIPTION_ = "BranchLengthOptimise";
		private final StoppingCriteria stopper_;
		private final AlgorithmCallback callback_;
		private final double[] branchLengths_;
		private final double[][] edgeParameters_;
		private final Connection[] allConnections_;
		public BranchLengthOptimiseAction(Connection[] allConnections, StoppingCriteria stopper, AlgorithmCallback callback) {
			this.stopper_ = stopper;
			stopper_.reset();
			this.callback_ = callback;
			this.allConnections_ = allConnections;
			this.branchLengths_ = new double[allConnections.length];
			this.edgeParameters_ = new double[allConnections.length][];
		}
		public double doAction(double currentScore, double desparationValue) {
			if (DEBUG) System.out.println("+++BranchLengthOptimise"); Date start = new Date();
			storeBranchData(allConnections_,branchLengths_,edgeParameters_);
			double result = simpleOptimiseLikelihood(stopper_, callback_);
			if (DEBUG) System.out.printf("---BranchLengthOptimise %d ms; %.4f\n",(new Date()).getTime()-start.getTime(),result);
			return result;
		}
		public boolean isActionSuccessful() {
			return true;
		}
		/**
		 *
		 * @return true
		 */
		public boolean undoAction() {
			restoreBranchData(allConnections_,branchLengths_,edgeParameters_);
			if (DEBUG) System.out.println("***BranchLengthOptimise undone");
			return true;
		}
		/**
		 * @return false
		 */
		public boolean isActionDeterministic() {
			return false;
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	}

// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
	/*
	 * Simultaneously optimises model and branch lengths (and edge parameters, if any.) 
	 * It alternates between optimising model and optimising branches,
	 * with reduced accuracy the first few iterations, tightening accuracy as it goes.
	 * If given a random number generator, it will randomizer the model before starting.
	 */
	private final class BranchLengthWithModelOptimiseAction implements UndoableAction, MultivariateFunction {
		private static final String DESCRIPTION_ = "BranchLengthWithModelOptimise";
		private final StoppingCriteria stopper_;
		private final AlgorithmCallback callback_;
		private final double[] branchLengths_;
		private final double[][] edgeParameters_;
		private final Connection[] allConnections_;
		private final MultivariateMinimum minimiser_;
		private final MinimiserMonitor monitor_;
		private final SubstitutionModel model_;
		private final double[] modelParameterStore_;
		private final double[] xvec_;
		private final int fxDigits_;
		private final int xDigits_;
		private final LHCalculator.External calculator_;
		private final MersenneTwisterFast random_; // MAY BE NULL! In which case we don't randomize at start.

		private final UnivariateMinimum um_;
		private final UniOptimisationHandler optimisationHandler_;
		private final ConstructionTool tool_;

		// Might make these into constructor parameters
		// We slowly ramp up optimizations to full accuracy. How long til this happens?
		// NOTE: I haven't experimented to find optimial values for these two parameters.
		private final int iterationsAtReducedAccuracy_ = 3;
		// How many of these iterations at each number of digits accuracy?
		private final int iterationsPerDigit_ = 1;
		// The first iteration will have iterationsAtReducedAccuracy/iterationsPerDigit digits less
		// accuracy than the final (xDigits, fxDigits) requested accuracy.
		
		public BranchLengthWithModelOptimiseAction(
				ConstructionTool tool, 
				Connection[] allConnections, 
				StoppingCriteria stopper, 
				AlgorithmCallback callback, 
				MultivariateMinimum minimiser, 
				MinimiserMonitor monitor, 
				SubstitutionModel model, 
				int fxDigits, 
				int xDigits,
				MersenneTwisterFast random
		) {
			this.stopper_ = stopper;
			stopper_.reset();
			this.tool_ = tool;
			this.callback_ = callback;
			this.allConnections_ = allConnections;
			this.branchLengths_ = new double[allConnections.length];
			this.edgeParameters_ = new double[allConnections.length][];
			this.minimiser_ = minimiser;
			this.monitor_ = monitor;
			this.model_ = model;
			this.calculator_ = tool.allocateNewExternalCalculator();
			this.modelParameterStore_ = new double[model.getNumParameters()];
			this.xvec_ = new double[model.getNumParameters()];
			this.xDigits_ = xDigits;
			this.fxDigits_ = fxDigits;
			this.um_ = new UnivariateMinimum();
			this.optimisationHandler_ = new UniOptimisationHandler(model_,tool_);
			this.random_ = random;
		}
		// 'evaluate' is only used in optimizing model parameters, not edge lengths (& parameters)
		public double evaluate(double[] xvec) {
			for(int i = 0 ; i < xvec.length ; i++) {
				model_.setParameter(xvec[i],i);
			}
			return -allConnections_[0].calculateLogLikelihood(model_,true,calculator_,tool_);
		}
		public int getNumArguments() { return xvec_.length; }
		public double getLowerBound(int n) { return model_.getLowerLimit(n); }
		public double getUpperBound(int n) { return model_.getUpperLimit(n); }
		public OrthogonalHints getOrthogonalHints() { return model_.getOrthogonalHints(); }

		public double doAction(double currentScore, double desparationValue) {
			if (DEBUG) System.out.println("+++BranchLengthWithModelOptimise"); Date start = new Date();
			storeBranchData(allConnections_,branchLengths_,edgeParameters_);
			// TODO: if Parameterized had the getAllParameters() method of ParameterizedUser, this would be one liner.
			for(int i = 0 ; i < modelParameterStore_.length ; i++) {
				modelParameterStore_[i] = model_.getParameter(i);
			}
			
			
			if (random_ == null) {
				// start search from current model
				System.arraycopy(modelParameterStore_,0,xvec_,0,xvec_.length);
			} else {
				// randomize model 
				int n = model_.getNumParameters();
				for (int i=0; i<n; i++) {
					// This is good for zero-to-one parameters (e.g.) but for zero-to-1000
					// I might prefer something which gives values near 1.
					xvec_[i] = random_.nextDouble(model_.getLowerLimit(i),model_.getUpperLimit(i));
				}
				// Arguably randomizing model and reinitializing edge parameters should be independently selectable.
				for (Connection c : allConnections_) {
					if (c.edgeParams_ != null) { c.edgeParams_.initializeParameters(); }
				}
			}
			
			stopper_.reset();
			double minimum = Double.POSITIVE_INFINITY;
					
			/*
			boolean firstTime = true;
			while(!stopper_.isTimeToStop()) {
				for(int i = 0 ; i < allConnections_.length ; i++) {
					Connection c = allConnections_[i];
					optimisationHandler_.setup(c,firstTime);
					firstTime = false;
					um_.findMinimum(c.getBranchLength(),optimisationHandler_);
					c.setBranchLength(um_.minx);
					//minimum = um_.fminx;
					minimum = minimiser_.findMinimum(this,xvec_,fxDigits_,xDigits_,monitor_);
					for(int x = 0 ; x < xvec_.length ; x++) {
						model_.setParameter(xvec_[x],x);
					}
					c.setBranchLength(um_.minx);
			    }
				stopper_.newIteration(minimum,minimum,true,true,callback_);		
			}
			*/
			
			int remainingIterationsAtReducedAccuracy = iterationsAtReducedAccuracy_;
			while (remainingIterationsAtReducedAccuracy>0 && !stopper_.isTimeToStop()) {
				int fewerDigits = (iterationsAtReducedAccuracy_ + iterationsPerDigit_ - 1)/ iterationsPerDigit_;
				// Step 1: optimise model
				// On rare occasions, ConjugateDirectionSearch fails with 'QR failed' or generating NaNs which kill EigenSystem.
				// As it is a pseudorandom process, if this occurs, try several restarts before giving up.
				boolean success = false;
				RuntimeException error = null;
				for (int i=0; !success && i<5; i++) {
					try {
						minimum = -minimiser_.findMinimum(this,xvec_,fxDigits_-fewerDigits,xDigits_-fewerDigits,monitor_);
					} catch (RuntimeException e) {
						if (minimiser_ instanceof ConjugateDirectionSearch) {
							System.out.printf("Caught ConjugateDirectionSearch error '%s' %d times, trying again\n",e.getMessage(),i+1);
							error = e;
						} else {
							throw e;
						}
						continue;
					}
					success = true;
				}
				if (!success) {
					System.out.println("ConjugateDirectionSearch failed too many times.");
					throw error;
				}
				// Step 2: optimise branches
				final int MAX_IT_AT_SAME_SCORE = 20;
				final double tolerance = Math.pow(10, -xDigits_+fewerDigits);
				StoppingCriteria branchStopper = StoppingCriteria.Utils.getNonExactUnchangedScore(MAX_IT_AT_SAME_SCORE, true, tolerance).newInstance();
				minimum = simpleOptimiseLikelihood(branchStopper, callback_);
				remainingIterationsAtReducedAccuracy--;
			}
			if (DEBUG) System.out.printf("---BranchLengthWithModelOptimise %d; %.4f ms\n",(new Date()).getTime()-start.getTime(),minimum); 
			return minimum;
		}
		public boolean isActionSuccessful() {
			return true;
		}
		/**
		 *
		 * @return true
		 */
		public boolean undoAction() {
			if (DEBUG) System.out.println("***BranchLengthWithModelOptimise undone"); 
			restoreBranchData(allConnections_,branchLengths_,edgeParameters_);
			for(int i = 0 ; i < modelParameterStore_.length ; i++) {
				model_.setParameter(modelParameterStore_[i],i);
			}
			return true;
		}
		/**
		 * @return false
		 */
		public boolean isActionDeterministic() {
			return false;
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	}

	// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
	private final static class ModelOptimiseAction implements UndoableAction, MultivariateFunction {
		private static final String DESCRIPTION_ = "ModelOptimise";
		private final ConnectionOrRoot treeAccess_;
		private final MultivariateMinimum minimiser_;
		private final MinimiserMonitor monitor_;
		private final SubstitutionModel model_;
		private final double[] modelParameterStore_;
		private final double[] xvec_;
		private final int fxDigits_;
		private final int xDigits_;
		private final LHCalculator.External calculator_;
		private final ConstructionTool tool_;
		public ModelOptimiseAction(ConnectionOrRoot treeAccess, MultivariateMinimum minimiser, MinimiserMonitor monitor, SubstitutionModel model, int fxDigits, int xDigits, ConstructionTool tool) {
			this.treeAccess_ = treeAccess;
			this.minimiser_ = minimiser;
			this.monitor_ = monitor;
			this.model_ = model;
			this.tool_ = tool;
			this.calculator_ = tool_.allocateNewExternalCalculator();

			this.modelParameterStore_ = new double[model.getNumParameters()];
			this.xvec_ = new double[model.getNumParameters()];
			this.xDigits_ = xDigits;
			this.fxDigits_ = fxDigits;
		}
		public double evaluate(double[] xvec) {
			for(int i = 0 ; i < xvec.length ; i++) {
				model_.setParameter(xvec[i],i);
			}
			// Line split to make it easier to follow progress with debugger
			double fx = -treeAccess_.calculateLogLikelihood(model_,true,calculator_, tool_);
			return fx;
		}
		public int getNumArguments() { return xvec_.length; }
		public double getLowerBound(int n) { return model_.getLowerLimit(n); }
		public double getUpperBound(int n) { return model_.getUpperLimit(n); }
		public OrthogonalHints getOrthogonalHints() { return model_.getOrthogonalHints(); }
		/**
		 * @return true
		 */
		public boolean isActionDeterministic() {
			return true;
		}
		public double doAction(double currentScore, double desparationValue) {
			if (DEBUG) System.out.println("+++ModelOptimise"); Date start = new Date();
			// Following section is for debugging only, so I can see the current state of the model and tree.
			// Entire section can be removed without damaging calculation at all.
			if (false) {
				double initialLogLikelihood = treeAccess_.calculateLogLikelihood(model_, true, calculator_, tool_);
				// need to construct an allConnections array from treeAccess_
				ArrayList<Connection> store = new ArrayList<Connection>();
				treeAccess_.getAllConnections(store);
				Connection[] allConnections = new Connection[store.size()];
				store.toArray(allConnections);
				String dummyString = unrootedMLSearcherToString(model_,allConnections);
				store = null;
			}
			
			for(int i = 0 ; i < modelParameterStore_.length ; i++) {
				modelParameterStore_[i] = model_.getParameter(i);
			}
			boolean success = false;
			double minimum = Double.MAX_VALUE;
			RuntimeException error = null;
			// On rare occasions, ConjugateDirectionSearch fails with 'QR failed' or generating NaNs which kill EigenSystem.
			// As it is a pseudorandom process, if this occurs, try several restarts before giving up.
			for (int i=0; !success && i<5; i++) {
				System.arraycopy(modelParameterStore_,0,xvec_,0,xvec_.length);
				try {
					minimum = -minimiser_.findMinimum(this,xvec_,fxDigits_,xDigits_,monitor_);
				} catch (RuntimeException e) {
					if (minimiser_ instanceof ConjugateDirectionSearch) {
						System.out.printf("Caught ConjugateDirectionSearch error '%s' %d times, trying again\n",e.getMessage(),i+1);
						error = e;
					} else {
						throw e;
					}
					continue;
				}
				success = true;
			}
			if (!success) {
				System.out.println("ConjugateDirectionSearch failed too many times.");
				throw error;
			}
			for(int i = 0 ; i < xvec_.length ; i++) {
				model_.setParameter(xvec_[i],i);
			}
			if (DEBUG) System.out.printf("---ModelOptimise %d ms; %.4f\n",(new Date()).getTime()-start.getTime(),minimum); 
			return minimum;
		}
		public boolean isActionSuccessful() {
			return true;
		}
		/**
		 *
		 * @return true
		 */
		public boolean undoAction() {
			if (DEBUG) System.out.println("***ModelOptimise undone");
			for(int i = 0 ; i < modelParameterStore_.length ; i++) {
				model_.setParameter(modelParameterStore_[i],i);
			}
			return true;
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	}

// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
	/**
	 *  Runs through all internal edges (in random order) and tries both NNIs on each one as it goes.
	 *  For each edge, keeps the best scoring result (unchanged or one of the two NNIs).
	 *  (Not intended for use with annealing or other methods which can accept down-hill steps. Use NNIAction with these.) 
	 */
	private static final class FullSweepNNIAction implements UndoableAction {
		private static final String DESCRIPTION_ = "FullSweepNNI";
		private final Connection[] allConnections_;
		private final Connection[] internalConnections_;
		private final NNIAction baseAction_;
		private final Markable subject_;
		private final MersenneTwisterFast random_;
		private final ConstructionTool tool_;
		private boolean lastActionSuccessful_ = false;
		private final double[] falseBranchLengths_; // non-persistant preallocated memory used in doSetupAction.
	
		public FullSweepNNIAction(Connection[] allConnections, Markable subject, NNIAction baseAction, MersenneTwisterFast random, ConstructionTool tool) {
			allConnections_ = allConnections;
			subject_ = subject;
			baseAction_ = baseAction;
			random_ = random;
			tool_ = tool;
			falseBranchLengths_ = new double[allConnections.length];
			// unrooted tree with L leafs has E=2L-3 edges, of which I=E-L=(E-3)/2 are internal edges.
			// Rooted   tree with L leafs has E=2L-2 edges, of which I=E-L=(E-2)/2 are internal edges.
			// for unrooted, we get the right answer due to rounding.
			this.internalConnections_ = new Connection[(allConnections.length-2)/2];
		}
		
		private void findInternalEdges() {
			int i=0;
			for (Connection c: allConnections_) {
				if (!c.isPendantEdge()) {
					internalConnections_[i++] = c;
				}
			}
			if (i!=internalConnections_.length) {
				throw new RuntimeException("Can't happen");
			}
			random_.shuffle(internalConnections_);
		}
		
		// TODO: much cut-and-paste from SweepNNIAction here - look to abstract this stuff out.
		// I can't easily just call SweepNNIAction multiple times because it also uses subject_.mark().
		public double doAction(double originalScore, double desparationValue) {
			if (DEBUG) System.out.println("+++FullSweepNNI"); Date start = new Date();
			subject_.mark();
			findInternalEdges();
			lastActionSuccessful_ = false; // haven't found an improvement yet
			double bestScore = originalScore;
			boolean improvementFromFalse = false;
			// Loop through all the internal edges
			for (Connection c : internalConnections_) {
				baseAction_.setTarget(c, false);
				double falseScore = baseAction_.doSetupAction(bestScore);
				if (falseScore > bestScore) {
					bestScore = falseScore;
					improvementFromFalse = true;
					for (int i=0; i<allConnections_.length; i++) {
						falseBranchLengths_[i] = allConnections_[i].getBranchLength();
					}
				}
				baseAction_.undoAction();
				baseAction_.setTarget(c, true);
				double trueScore = baseAction_.doSetupAction(bestScore);
				if (trueScore > bestScore) {
					// second attempted NNI is the best
					bestScore = trueScore;
					lastActionSuccessful_ = true;
				} else {
					// second NNI was not the best
					baseAction_.undoAction();
					if (improvementFromFalse) {
						// but first NNI was still an improvement, so revert to that.
						// baseAction_.doSetupAction would needlessly redo branch optimisation
						lastActionSuccessful_ = true;
						c.doNNI(false);
						c.setup(tool_,allConnections_);
						for (int i=0; i<allConnections_.length; i++) {
							allConnections_[i].setBranchLength(falseBranchLengths_[i]);
						}
					} // if improvemementFromFalse
				} // if trueScore>bestScore else
			} // loop through internalConnections
			if (DEBUG) System.out.printf("---FullSweepNNI %d ms; %.4f\n",(new Date()).getTime()-start.getTime(), bestScore);
			return bestScore;
		} // doAction
		
		public boolean isActionDeterministic() {
			return false;
		}
		public boolean isActionSuccessful() { return lastActionSuccessful_; }
		public boolean undoAction() {
			if (DEBUG) System.out.println("***FullSweepNNI undone"); 
			if(lastActionSuccessful_) {
				subject_.undoToMark();
				return true;
			} else {
				throw new RuntimeException("Illegal operation : undoLast() called when last operation invalid (may already have been undone)");
			}
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	
	} // FullSweepNNIAction
	
// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
	// like NNIAction, except it tries both possible NNIs on an edge instead of a single random one, and 
	// accept best result, if it is better than current score. (Not intended for use with annealing or other
	// methods which can accept down-hill steps. Use NNIAction with these.)
	private static final class SweepNNIAction implements UndoableAction {
		private static final String DESCRIPTION_ = "SweepNNI";
		private final Connection[] allConnections_;
		private final NNIAction baseAction_;
		private final Markable subject_;
		private final MersenneTwisterFast random_;
		private final ConstructionTool tool_;
		private boolean lastActionSuccessful_ = false;
		private Connection connection_; // the target edge to try both NNIs on
		private final double[] falseBranchLengths_; // non-persistant preallocated memory used in doSetupAction. 
		
		public SweepNNIAction(Connection[] allConnections, Markable subject, NNIAction baseAction, MersenneTwisterFast random, ConstructionTool tool) {
			allConnections_ = allConnections;
			subject_ = subject;
			baseAction_ = baseAction;
			random_ = random;
			tool_ = tool;
			falseBranchLengths_ = new double[allConnections.length];
		}
		
		public void setRandomTarget() {
			do {
				connection_ = allConnections_[random_.nextInt(allConnections_.length)];
			} while (connection_.isPendantEdge()); // select only internal edges
		}
		// Following routine currently not ever used.
		public void setTarget(Connection connection) {
			// With sensible planning, should never call on a pendant edge. This error encourages
			// sensible planning. However, removing the error and letting the action fail 
			// would only harm efficiency.
			if (connection_.isPendantEdge()) {
				throw new IllegalArgumentException("Can only NNI on an internal edge");
			}
			connection_ = connection;
		}

		public double doAction(double currentScore, double desparationValue) {
			if (DEBUG) System.out.println("+++SweepNNI"); Date start = new Date();
			setRandomTarget();
			double result =  doSetupAction(currentScore);
			if (DEBUG) System.out.printf("---SweepNNI %d ms; %.4f\n",(new Date()).getTime()-start.getTime(), result); 
			return result;
		}
		
		public double doSetupAction(double currentScore) {
			
			subject_.mark();
			lastActionSuccessful_ = true; // Failure is not an option - we'll take one NNI or the other.
			// First attempt: NNI 'false'
			baseAction_.setTarget(connection_, false);
			double falseScore = baseAction_.doSetupAction(currentScore);
			// And remember the results
			for (int i=0; i<allConnections_.length; i++) {
				falseBranchLengths_[i] = allConnections_[i].getBranchLength();
			}
			// undo for the next attempt
			baseAction_.undoAction();
			
			// Second attempt: NNI 'true'
			baseAction_.setTarget(connection_, true);
			double trueScore = baseAction_.doSetupAction(currentScore);
			if (trueScore > falseScore) {
				// second attempted NNI is the best, 
				if (DEBUG) System.out.println("...SweepNNI takes second option");
				return trueScore;
			} else {
				// second NNI was not the best. Undo the second NNI:
				baseAction_.undoAction();
				// and redo the first one (skipping repeating tedious branch length optimisation.)
				connection_.doNNI(false);
				connection_.setup(tool_,allConnections_);
				for (int i=0; i<allConnections_.length; i++) {
					allConnections_[i].setBranchLength(falseBranchLengths_[i]);
				}
				if (DEBUG) System.out.println("...SweepNNI takes first option");
				return falseScore;
			}
		}
		
		public boolean isActionDeterministic() {
			return false;
		}
		public boolean isActionSuccessful() { return lastActionSuccessful_; }
		public boolean undoAction() {
			if (DEBUG) System.out.println("***SweepNNI undone"); 
			if(lastActionSuccessful_) {
				subject_.undoToMark();
				return true;
			} else {
				throw new RuntimeException("Illegal operation : undoLast() called when last operation invalid (may already have been undone)");
			}
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	} // class SweepNNIAction
	
// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
	// Can be called directly, or also used indirectly by SweepNNIAction and FullSweepNNIAction.
	// Note: may not use 'mark' on Connections etc. as SweepNNIAction and FullSweepNNIAction use marking.
	private static final class NNIAction implements UndoableAction {
		private static final String DESCRIPTION_ = "NNI";
		private boolean bottomLeftSwapsWithTopLeft_; // which of the two possible NNIs to attempt
		private Connection connection_;
		private final Connection[] allConnections_;
		private final double[] branchLengths_;
		private final Assessor assessor_;
		private boolean lastActionSuccessful_ = false;
		private final MersenneTwisterFast random_;
		private final ConstructionTool tool_;
		public NNIAction(Connection[] allConnections, Assessor assessor, MersenneTwisterFast r, ConstructionTool tool) {
			this.allConnections_ = allConnections;
			this.assessor_ = assessor;
			this.random_ = r;
			this.tool_ = tool;
			this.branchLengths_ = new double[allConnections.length];
		}
		
		public void setRandomTarget() {
			do {
				connection_ = allConnections_[random_.nextInt(allConnections_.length)];
			} while (connection_.isPendantEdge()); // select only internal edges
			bottomLeftSwapsWithTopLeft_ = random_.nextBoolean();
		}
		public void setTarget(Connection connection, boolean direction) {
			// With sensible planning, should never call on a pendant edge. This error encourages
			// sensible planning. However, removing the error and letting the action fail 
			// would only harm efficiency.
			if (connection.isPendantEdge()) {
				throw new IllegalArgumentException("Can only NNI on an internal edge");
			}
			connection_ = connection;
			bottomLeftSwapsWithTopLeft_ = direction;
		}
		public double doAction(double currentScore, double desparationValue) {
			if (DEBUG) System.out.println("+++NNI"); Date start = new Date();
			setRandomTarget();
			double result =  doSetupAction(currentScore);
			if (DEBUG) System.out.printf("---NNI %d ms; %.4f\n",(new Date()).getTime()-start.getTime(), result); 
			return result;
		}
		
		public double doSetupAction(double currentScore) {
			lastActionSuccessful_ = connection_.doNNI(bottomLeftSwapsWithTopLeft_);
			if(lastActionSuccessful_) {
				for(int i = 0 ; i < allConnections_.length ; i++) {
					branchLengths_[i] = allConnections_[i].getBranchLength();
				}
				connection_.setup(tool_,allConnections_);
				//return assessor_.getCurrentValue();
				double result = assessor_.getCurrentValue();
				return result;
			} else {
				throw new RuntimeException("NNI action called on illegal edge");
				//return currentScore;
			}
		}
		/**
		 * @return false
		 */
		public boolean isActionDeterministic() {
			return false;
		}
		public boolean isActionSuccessful() { return lastActionSuccessful_; }
		/**
		 * @return true
		 */
		public boolean undoAction() {
			if (DEBUG) System.out.println("***NNI undone"); 
			if(lastActionSuccessful_) {
				connection_.doNNI( bottomLeftSwapsWithTopLeft_);
				for( int i = 0; i<allConnections_.length; i++ ) {
					allConnections_[i].setBranchLength( branchLengths_[i] );
				}
				connection_.setup(tool_,allConnections_);
				lastActionSuccessful_=false;
				return true;
			} else {
				throw new RuntimeException("Illegal operation : undoLast() called when last operation invalid (may already have been undone)");
			}
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	}
// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
	private static final class SweepSPRAction implements UndoableAction {
		private static final String DESCRIPTION_ = "SweepSPRA";
		private Connection toMove_ = null ;
		private SPRAction baseAction_;

		private boolean lastActionSuccessful_ = false;

		private final Connection[] shuffledConnections_;
		private final MersenneTwisterFast random_;

		/**
		 * @param solution a reference to the UnconstrainedOptiser
		 * @param assessor a means of assessing the solution (assumes gives true likelihood)
		 * @note I choose to use a static inner class and have funny references because of personal style (I like writing inner classes I can ship out to separate files if I need to)
		 */
		public SweepSPRAction(Connection[] allConnections,  SPRAction baseAction, MersenneTwisterFast random) {
			this.baseAction_ = baseAction;
			this.random_ = random;
			this.shuffledConnections_ = new Connection[allConnections.length];
			System.arraycopy(allConnections,0,shuffledConnections_,0,allConnections.length);
		}
		public void setRandomTarget() {
			Connection target;
			do {
				target = shuffledConnections_[random_.nextInt(shuffledConnections_.length)];
			} while (target.hasDirectConnectionToRoot());
			setTarget(target);
		}
		public void setTarget(Connection toMove) {
			this.toMove_ = toMove;
			shuffle(shuffledConnections_);
		}
		/**
		 * @return false
		 */
		public boolean isActionDeterministic() {
			return false;
		}
		private final void shuffle(Connection[] cs) {
			for(int i = 0 ; i < cs.length ; i++) {
				int j = random_.nextInt(cs.length-i)+i;
				Connection t = cs[i];
				cs[i] = cs[j];
				cs[j] = t;
			}
		}
		public double doAction(double originalScore, double desparationValue) {
			if (DEBUG) System.out.println("+++SweepSPRA"); Date start = new Date();
			setRandomTarget();
			//return doSetupAction(originalScore);
			double result = doSetupAction(originalScore);
			if (DEBUG) System.out.printf("---SweepSPRA %d ms; %.4f\n",(new Date()).getTime()-start.getTime(),result);
			return result;
		}
		public double doSetupAction(double originalScore) {
			Connection best = null;
			double bestScore = originalScore;
			for(int i = 0 ; i < shuffledConnections_.length; i++) {
				Connection c = shuffledConnections_[i];
				if(c!=toMove_) {
					baseAction_.setTarget(toMove_,c);
					double score = baseAction_.doSetupAction(originalScore);
					if(baseAction_.isActionSuccessful()) {
						if(score>bestScore) {
							best = c;
							bestScore = score;
						}
						baseAction_.undoAction();
					}
				}
			}
			if(best!=null) {
				baseAction_.setTarget(toMove_,best);
				lastActionSuccessful_ = true;
				//We assume that if the action worked before it will work now
				return baseAction_.doSetupAction(originalScore);
			}
			lastActionSuccessful_ = false;
			//When we fail the score does not matter
			return originalScore;
		}
		public boolean isActionSuccessful() { return lastActionSuccessful_; }
		public boolean undoAction() {
			if(lastActionSuccessful_) {
				return baseAction_.undoAction();
			} else {
				throw new RuntimeException("Illegal operation : undoLast() called when last operation invalid (may already have been undone)");
			}
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	}

// - - - - -- - - -  - - -- - - -  - -- - - - - - - - - - - -  -
	private static final class FullSweepSPRAction implements UndoableAction {
		private static final String DESCRIPTION_ = "FullSweepSPRA";
		private SPRAction baseAction_;

		private final Connection[] allConnections_;
		private boolean lastActionSuccessful_ = false;

		/**
		 * @param solution a reference to the UnconstrainedOptiser
		 * @param assessor a means of assessing the solution (assumes gives true likelihood)
		 * @note I choose to use a static inner class and have funny references because of personal style (I like writing inner classes I can ship out to separate files if I need to)
		 */
		public FullSweepSPRAction(Connection[] allConnections,  SPRAction baseAction) {
			this.baseAction_ = baseAction;
			this.allConnections_ = allConnections;
		}

		public double doAction(double originalScore, double desparationValue) {
			if (DEBUG) System.out.println("+++FullSweepSPRA"); Date starttime = new Date();
			Connection bestStart = null;
			Connection bestEnd = null;
			double bestScore = originalScore;
			for(int i = 0 ; i < allConnections_.length; i++) {
				Connection start = allConnections_[i];
				if (start.hasDirectConnectionToRoot()) { continue; } // Connections directly below the root are not valid pruning points.
				for(int j = i+1 ; j < allConnections_.length; j++) {
					Connection end = allConnections_[i];
					baseAction_.setTarget(start,end);
					double score = baseAction_.doSetupAction(originalScore);
					if(baseAction_.isActionSuccessful()) {
						if(score>bestScore) {
							bestStart = start;
							bestEnd = end;
							bestScore = score;
						}
						baseAction_.undoAction();
					}
				}
			}
			if(bestStart!=null) {
				baseAction_.setTarget(bestStart,bestEnd);
				lastActionSuccessful_ = true;
				//We assume that if the action worked before it will work now
				
				double result = baseAction_.doSetupAction(originalScore);
				if (DEBUG) System.out.printf("---FullSweepSPRA %d ms; %.4f\n",(new Date()).getTime()-starttime.getTime(),result); 
				return result;
			}
			lastActionSuccessful_ = false;
			//When we fail the score does not matter
			return originalScore;
		}
		/**
		 * @return true
		 */
		public boolean isActionDeterministic() {
			return true;
		}
		public boolean isActionSuccessful() { return lastActionSuccessful_; }
		public boolean undoAction() {
			if (DEBUG) System.out.println("***FullSweepSPRA undone");
			if(lastActionSuccessful_) {
				return baseAction_.undoAction();
			} else {
				throw new RuntimeException("Illegal operation : undoLast() called when last operation invalid (may already have been undone)");
			}
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	}
// - - -- - - - - - -- - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
	private static final class SPRAction implements UndoableAction {
		private static final String DESCRIPTION_ = "SPRA";
		private Connection toRemove_ = null ;
		private Connection attachmentPoint_ = null;
		private ConnectionOrRoot reattachmentPoint_ = null;
		private final Connection[] allConnections_;
		private final Connection[] store_ = new Connection[3];

		private boolean lastActionSuccessful_ = false;
		private final Assessor assessor_;
		private final Markable subject_;
		private final MersenneTwisterFast random_;
		private final ConstructionTool tool_;
		public SPRAction(Connection[] allConnections, Markable subject, Assessor assessor, MersenneTwisterFast random, ConstructionTool tool) {
			this.allConnections_ = allConnections;
			this.subject_ = subject;
			this.assessor_ = assessor;
			this.tool_ = tool;
			this.random_ = random;
		}
		public void setRandomTargets() {
			do {
				int i = random_.nextInt(allConnections_.length);
				int j = random_.nextInt(allConnections_.length-1);
				if(j>=i) { j++; }
				this.toRemove_ = allConnections_[i];
				this.attachmentPoint_ = allConnections_[j];
			// bad choice of Connections if they share an InternalNode, or if toRemove is adjacent to Root. 
			} while (toRemove_.hasDirectConnection(attachmentPoint_) || toRemove_.hasDirectConnectionToRoot());
		}
		public void setTarget(Connection toRemove, Connection attachmentPoint) {
			if (toRemove.hasDirectConnectionToRoot()) { throw new RuntimeException("In SPR: Invalid prune edge: direct child of root"); }
			this.toRemove_ = toRemove;
			this.attachmentPoint_ = attachmentPoint;
		}
		public boolean isActionSuccessful() { return lastActionSuccessful_; }
		/**
		 * Peform an action based on setup connections (ie must have called setTarget() already)
		 * @param currentScore
		 * @return the new score if successful
		 */
		public double doSetupAction(double currentScore) {
			subject_.mark();
			reattachmentPoint_ = toRemove_.attachTo(attachmentPoint_,store_);
			lastActionSuccessful_ = (reattachmentPoint_!=null);
			if(lastActionSuccessful_) {
				toRemove_.setup(tool_,allConnections_);
			}

			return(lastActionSuccessful_? assessor_.getCurrentValue() : currentScore);
		}
		public double doAction(double currentScore, double desparationValue) {
			if (DEBUG) System.out.println("+++SPR"); Date start = new Date();
			setRandomTargets();
			//return doSetupAction(currentScore);
			double result =  doSetupAction(currentScore);
			if (DEBUG) System.out.printf("---SPR %d ms; %.4f\n",(new Date()).getTime()-start.getTime(), result);
			return result;
		}
		/**
		 * @return false
		 */
		public boolean isActionDeterministic() {
			return false;
		}
		public boolean undoAction() {
			// PossiblyRootedMLSearcher.undoToMark() restores the edge lengths to old values.
			if (DEBUG) System.out.println("***SPR undone"); 
			if(lastActionSuccessful_) {
				subject_.undoToMark();
				return true;
			} else {
				throw new RuntimeException("Undo last called when last operation not successful");
			}
		}
		public void printStats(PrintWriter out) {};
		public String getDescription() { return DESCRIPTION_; }
	}
// -=-=-=-=-=-==--==--==--=-=-==-=--=-=-=-==-=-=--==-=-=--==--==-=-=--==--==--=
// == Static Utility Methods ===
// =--=-==-=--=-==--=-==--=-=-=-=-==-=-=--=-==-=-=-=--==-=-=-=-=--==-=-=-=-=-=-
	
	 // MDW: Recursive routine used in creating the shadow tree.
	private static final UNode createUNode(Node n, Connection parentConnection, ConstructionTool tool) {
		if(n.isLeaf()) {
			return new LeafNode(n,parentConnection, tool);
		}
		return new InternalNode(n,parentConnection,tool);
	}

	private static final UNode createUNode(String[] leafNames, Connection parentConnection, ConstructionTool tool, MersenneTwisterFast r) {
			if(leafNames.length==1) {
				return new LeafNode(leafNames[0],parentConnection, tool);
			}
			return new InternalNode(leafNames,parentConnection,tool,r);
		}



// -=-=-=-=-=-==--==--==--=-=-==-=--=-=-=-==-=-=--==-=-=--==--==-=-=--==--==--=
// == Workspace ===
// =--=-==-=---==--=-==--=-=-=-=-==-=-=--=-==-=-=-=--==-=-=-=-=--==-=-=-=-=-=-
	private static final class Workspace {
		private final double[][][][] conditionalLikelihoodTables_;
		private final double[] endStateProbabilityStore_;
		private final double[][][] transitionProbabilityStore_;
		private final int numberOfCategories_;
		private final int numberOfStates_;
		private final boolean[] locks_;
		public Workspace(int capacity, int maxNumberOfPatterns, ConstructionTool tool) {
			this.numberOfCategories_ = tool.getNumberOfTransitionCategories();
			numberOfStates_ = tool.getNumberOfStates();

			this.conditionalLikelihoodTables_ =
					new double
							[capacity]
							[numberOfCategories_]
							[maxNumberOfPatterns]
							[tool.getNumberOfStates()];
		 this.endStateProbabilityStore_ = new double[numberOfStates_];
		 this.transitionProbabilityStore_ = new double[numberOfCategories_][numberOfStates_][numberOfStates_];
		 this.locks_ = new boolean[capacity];
		}
		public final int getNumberOfCategories() { return numberOfCategories_; }
		public final int getNumberOfStates() { return numberOfStates_; }
		public final double[] getEndStateProbabilityStore() { return endStateProbabilityStore_; }
		public final int obtainLock() {
			 for(int i = 0 ; i < locks_.length ; i++) {
				 if(!locks_[i]) {
					 locks_[i]=true;
					 return i;
				 }
			 }
			 throw new RuntimeException("Assertion error : no locks available");
		 }
		 public final double[][][] getTransitionProbabilityStore() {
			 return transitionProbabilityStore_;
		 }
		 public final void returnLock(int lockNumber, double[][][] probabilities) {
			 if(locks_[lockNumber]) {
				 if(conditionalLikelihoodTables_[lockNumber]!=probabilities) {
					 throw new IllegalArgumentException("Table mismatch");
				 }
				 locks_[lockNumber]=false;
			 } else {
				 throw new IllegalArgumentException("Lock already unlocked:"+lockNumber);
			 }
		 }
		 public final void freeLock(int lockNumber) {
			 if(locks_[lockNumber]) {
				 locks_[lockNumber]=false;
			 } else {
				 throw new IllegalArgumentException("Lock already unlocked:"+lockNumber);
			 }
		 }
		 public final double[][][] getConditionalProbabilityTable(int lockNumber) {
			 if(locks_[lockNumber]) {
				 return conditionalLikelihoodTables_[lockNumber];
			 }else {
				 throw new IllegalArgumentException("Accessing unlocked table:"+lockNumber);
			 }
		 }
	}
// -=-=-=-=-=-==--==--==--=-=-==-=--=-=-=-==-=-=--==-=-=--==--==-=-=--==--==--=
// == UNode ===
// =--=-==-=--=-==--=-==--=-=-=-=-==-=-=--=-==-=-=-=--==-=-=-=-=--==-=-=-=-=-=-
	/*
	 * MDW: A node in the shadow tree. Implementations are InternalNode and LeafNode.
	 */
	private static interface UNode extends Markable {

		public String toString(Connection caller);
		public PatternInfo getPatternInfo(ConnectionOrRoot caller);
		public Node buildPALNode(double branchLength, ConnectionOrRoot caller);
		public void testLikelihood(ConnectionOrRoot caller, SubstitutionModel model, ConstructionTool tool);
		public void getAllConnections(ArrayList<Connection> store, ConnectionOrRoot caller);
		public int getIndex();
		public UNode createAlteredCopy(Connection attachmentPoint, Node newSubtree, Connection originalParentConnection, Connection parentConnection, ConstructionTool tool);
		public UNode createAlteredCopy(Connection originalParentConnection, Connection parentConnection, ConstructionTool tool);
		public boolean hasDirectConnection(Connection c);
		public boolean hasConnection(Connection c, ConnectionOrRoot caller);
		public ModelNodeParameters getNodeParameters();
		public void setNodeParameters(ModelNodeParameters mnp);
		public ModelNodeParameters[] getAdjacentNodeParameters();
		public ModelEdgeParameters[] getAdjacentEdgeParameters();
		public int[] getAdjacentNodeBackpointers();
		public int[] getAdjacentEdgeBackpointers();
		public int getBackpointer(ConnectionOrRoot caller);
		
		public void mark();
		public void undoToMark();

		public void instruct(UnrootedTreeInterface.UNode node, Connection callingConnection);

		public boolean isLeaf();
		public boolean hasLabel(String label);
		public void setAnnotation(Object annotation);
		public Object getAnnotation();
		public String getLabel();

		/**
		 * Instruct the node to extract itself from the two connections that aren't the caller
		 * One of the other two connections will become redunant.
		 * @return the redundant connection, or null of this node can't extract
		 */
		public Connection extract(Connection caller);

		/**
		 * @return the left connection with reference to the caller
		 * @note can return null if not possible (if leaf)
		 */
		public ConnectionOrRoot getLeft(Connection caller);
		/**
		 * @return the right connection with reference to the caller
		 * @note can return null if not possible (if leaf)
		 */
		public ConnectionOrRoot getRight(Connection caller);

		/**
		 * Set the connections to this node
		 * @param store a temporary store of the connections - node must copy references, do not use store
		 * @param number the number of connections to look at (ignore the length of store)
		 */
		public void setConnections(Connection[] store, int number);
		/**
		 * Should preserver tree integrity
		 */
		public void swapConnection(Connection original, UNode nodeToReplace, Connection newConnection);
		/**
		 * Should not do anything but swap connections around
		 */
		public void swapConnection(Connection original,Connection newConnection);

		/**
		 * Recurse to all neighbours but caller
		 * @return the maximum number of patterns from any neighbour
		 */
		public int rebuildPattern( ConstructionTool tool, Connection caller, boolean firstPass);
		/**
		 * Recurse to all neighbours
		 * @return the maximum number of patterns from any neighbour
		 */
		public int rebuildPattern( ConstructionTool tool);
		/**
		 * To be used by nodes that cannot properly do a rebuildPattern(tool) call, so they redirect to the other end of a connection
		 */
		public int redirectRebuildPattern( ConstructionTool tool);

		public ConditionalProbabilityStore getFlatConditionalProbabilities( SubstitutionModel model, boolean modelChanged,  ConnectionOrRoot callingConnection, int depth);
		public ConditionalProbabilityStore getFlatConditionalProbabilities( SubstitutionModel model, boolean modelChanged,  Connection callingConnection, LHCalculator.External external, ConditionalProbabilityStore resultStore);

		public ConditionalProbabilityStore getLeftExtendedConditionalProbabilities( SubstitutionModel model, boolean modelChanged,  Connection callingConnection, LHCalculator.External external, ConditionalProbabilityStore resultStore);
		public ConditionalProbabilityStore getRightExtendedConditionalProbabilities( SubstitutionModel model, boolean modelChanged,  Connection callingConnection, LHCalculator.External external, ConditionalProbabilityStore resultStore);
		/**
		 *
		 * @param caller
		 * @return Get the pattern info for the relative left (from the caller's perspective), or null if not left pattern info
		 */
		public PatternInfo getLeftPatternInfo(Connection caller);
		/**
		 *
		 * @param caller
		 * @return Get the pattern info for the relative right (from the caller's perspective), or null if not right pattern info
		 */
		public PatternInfo getRightPatternInfo(Connection caller);


		public ConditionalProbabilityStore getExtendedConditionalProbabilities( double distance, SubstitutionModel model, boolean modelChanged,  Connection callingConnection, int depth);
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( double distance, SubstitutionModel model, boolean modelChanged,  Connection callingConnection, LHCalculator.External external, ConditionalProbabilityStore resultStore);

		public void getLeafNames(ArrayList store, Connection caller);

		public void getSplitInformation(int[] splitStore, String[] leafNames, int splitIndex, Connection caller);
		/**
		 * Guarantees that subtree as seen from caller has all Connections correctly pointing towards the root. 
		 * @param caller
		 * @return true if subtree contains the root, false otherwise.
		 */
		public boolean fixConnectionDirections(ConnectionOrRoot caller); 
	}

// -=-=-=-=-=-==--==--==--=-=-==-=--=-=-=-==-=-=--==-=-=--==--==-=-=--==--==--=
// == InternalNode ===
// =--=-==-=--=-==--=-==--=-=-=-=-==-=-=--=-==-=-=-=--==-=-=-=-=--==-=-=-=-=-=-
	/*
	 * MDW:
	 * Internal node in an unrooted binary tree, so it has three adjacent Connections (edges)
	 * 
	 */
	private static final class InternalNode implements UNode {
		private final ConnectionOrRoot[] connections_ = new ConnectionOrRoot[3];
		private final ConnectionOrRoot[] markConnections_ = new ConnectionOrRoot[3]; // MDW: we can remember one past state and revert to it

		// MDW: One PatternInfo for each adjacent Connection. Each Connection sees the node as the
		// root of a subtree. The PatternInfo holds the patterns on this subtree.
		private final PatternInfo[] patterns_; 

		// MDW if we came from Connection 0, then left Connection is #1, right Connection is #2, etc.
		private static final int[] LEFT_LOOKUP = {	1 , 0, 0	}; 
		private static final int[] RIGHT_LOOKUP = {	2 , 2, 1	};

		private final int index_;  // MDW Each node within the tree gets a unique index number
		private final LHCalculator.Internal calculator_;
		private boolean topologyChangedSinceLastFlat_ = true;
		private boolean topologyChangedSincleLastExtended_ = true;
		private ModelNodeParameters nodeParams_; // MDW used by MosaicModel, likely unneeded by anything else.

		private Object annotation_ = null;

		/**
		 * The random tree constructor
		 * @param leafNames The names of the leafs remaining to be created
		 * @param parentConnection The connection that the recursion is "coming from"
		 * @param tool A tool to aid in construction
		 * @param r A method for getting random numbers used in determining branching
		 */
		private InternalNode(String[] leafNames , ConnectionOrRoot parentConnection, ConstructionTool tool, MersenneTwisterFast r) {
			this.connections_[0] = parentConnection;
			this.index_ = tool.allocateNextUNodeIndex(this);
			this.calculator_ = tool.allocateNewInternalCalculator();
			String[][] split = SearcherUtils.split(leafNames,r);
			this.connections_[1] = new Connection(split[0],this, tool,r);
			this.connections_[2] = new Connection(split[1],this, tool,r);

			final int numberOfSites = tool.getNumberOfSites();
			this.patterns_ = new PatternInfo[] {
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true)
				};
			checkConnections();
		}
		/**
		 * The altered tree constructor
		 * @param parentConnection The connection that the recursion is "coming from"
		 * @param attachmentPoint The connection that the new sub tree will be attached to
		 * @param originalSubTree
		 * @param originalSubTreeBranchLength
		 * @param originalSubTreeParentConnection
		 * @param appendedSubTree
		 * @param appendedSubTreeBranchLength
		 * @param tool to aid in construction
		 */
		public InternalNode(Connection parentConnection ,Connection attachmentPoint, UNode originalSubTree, double originalSubTreeBranchLength,  Connection originalSubTreeParentConnection, UNode appendedSubTree, double appendedSubTreeBranchLength, ConstructionTool tool) {
			//The first connection is the parent connection
			this.connections_[0] = parentConnection;

		  //The second connection
			this.connections_[1] = new Connection(appendedSubTree, this, appendedSubTreeBranchLength, tool);
			this.connections_[2] = new Connection(originalSubTree,originalSubTreeParentConnection, this, originalSubTreeBranchLength, tool);

			this.index_ = tool.allocateNextUNodeIndex(this);
			this.calculator_ = tool.allocateNewInternalCalculator();

			final int numberOfSites = tool.getNumberOfSites();
			this.patterns_ = new PatternInfo[] {
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true)
				};
			checkConnections();
		}
		
		// MDW: Recursively create shadow subtree from an original Node 
		public InternalNode(Node i, ConnectionOrRoot parentConnectionOrRoot, ConstructionTool tool) {
			// MDW code change: Attach "ModelEdgeParameters" attribute object to parent Connection
			// (check for parentConnection.edgeParams_ == null is only relevant for the
			// root connection - no other connection gets more than one chance to get an edgeParams.)
			if (parentConnectionOrRoot instanceof Connection) { // could be Root instead
				Connection parentConnection = (Connection)parentConnectionOrRoot;
				if (i instanceof AttributeNode && parentConnection.edgeParams_ == null) {
					parentConnection.edgeParams_ = (ModelEdgeParameters) ((AttributeNode) i).getAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL);
					nodeParams_                  = (ModelNodeParameters) ((AttributeNode) i).getAttribute(ModelNodeParameters.ATTRIBUTE_LABEL);
				}
			}
			
			this.connections_[0] = parentConnectionOrRoot;
			// MDW these next two are recursive calls
			this.connections_[1] = new Connection(i.getChild(0),this, tool);
			this.connections_[2] = new Connection(i.getChild(1),this, tool);
			this.index_ = tool.allocateNextUNodeIndex(this);
			this.calculator_ = tool.allocateNewInternalCalculator();
		  if(calculator_!=null) {
				final int numberOfSites = tool.getNumberOfSites();
				// MDW: Allocate PatternInfo storage (not used yet)
				this.patterns_ = new PatternInfo[] {
												 new PatternInfo( numberOfSites, true ),
												 new PatternInfo( numberOfSites, true ),
												 new PatternInfo( numberOfSites, true )
				};
			} else {
				this.patterns_ = null;
			}
		}
		/**
		 * The cloning with attachment constructor
		 * @param original The node we are replacing
		 * @param attachmentPoint The original connection that will be the attachment point for the new sub tree
		 * @param newSubtree A PAL node representing new sub tree
		 * @param originalParentConnection The orginal parent connection for the original node
		 * @param parentConnection The new parent connection for us
		 * @param tool to aid in construction
		 */
		private InternalNode(InternalNode original, Connection attachmentPoint, Node newSubtree,  Connection originalParentConnection, Connection parentConnection, ConstructionTool tool) {
		  this.connections_[0] = parentConnection;
		    // MDW: I'm not 100% sure these casts are safe: may need to alter the 'Connection' constructor to take a ConnectionOrRoot instead. 
			final Connection originalLeft = (Connection)original.getLeft(originalParentConnection);
			final Connection originalRight = (Connection)original.getRight(originalParentConnection);

			this.connections_[1] = new Connection(originalLeft,attachmentPoint,newSubtree,original, this, tool);
			this.connections_[2] = new Connection(originalRight,attachmentPoint,newSubtree,original, this, tool);
			this.index_ = tool.allocateNextUNodeIndex(this);
			this.calculator_ = tool.allocateNewInternalCalculator();

			final int numberOfSites = tool.getNumberOfSites();
			this.patterns_ = new PatternInfo[] {
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true)
				};
			checkConnections();
		}
		/**
		 * The constructor for the internal node that is added to attach a new sub tree with the altered tree stuff
		 * @param parentConnection The newly create parent connection (the root of the recursion, the caller), forming one of the children from this node
		 * @param newSubtree The PAL node of the sub tree which forms one of the children from this node
		 * @param originalParentConnection The original parent connection
		 * @param originalOtherChild The remaining child, to be cloned from an original
		 * @param otherChildLength the length of the connection to the "other child"
		 * @param tool to aid in construction
		 */
		private InternalNode(Connection parentConnection, Node newSubTree, UNode originalOtherChild,   Connection originalOtherChildParentConnection, double otherChildLength, ConstructionTool tool) {
		  this.index_ = tool.allocateNextUNodeIndex(this);
			this.calculator_ = tool.allocateNewInternalCalculator();
			final int numberOfSites = tool.getNumberOfSites();
			this.patterns_ = new PatternInfo[] {
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true)
				};
			//From parent
			this.connections_[0] = parentConnection;

			//From new sub tree
			this.connections_[1] = new Connection(newSubTree,this, tool);

			//From other child
			this.connections_[2] = new Connection(originalOtherChild,originalOtherChildParentConnection,this, otherChildLength, tool);
			checkConnections();
		}
		/**
		 * Cloning constructor
		 * @param originalNode The original internal node that we are replacing
		 * @param originalParentConnection The original parent connection
		 * @param parentConnection The replacement for the parent connect
		 * @param tool to aid in construction
		 */
		private InternalNode(InternalNode originalNode, Connection originalParentConnection, Connection parentConnection, ConstructionTool tool) {
		  this.connections_[0] = parentConnection;
			final Connection originalLeft = (Connection)originalNode.getLeft(originalParentConnection);
			final Connection originalRight = (Connection)originalNode.getRight(originalParentConnection);

			this.connections_[1] = new Connection(originalLeft,originalNode, this, tool);
			this.connections_[2] = new Connection(originalRight,originalNode, this, tool);
			this.index_ = tool.allocateNextUNodeIndex(this);
			this.calculator_ = tool.allocateNewInternalCalculator();

			final int numberOfSites = tool.getNumberOfSites();
			this.patterns_ = new PatternInfo[] {
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true),
					new PatternInfo(numberOfSites,true)
				};
			checkConnections();
		}
		/**
		 * Sanity check: throws an exception if node has more than one parent.
		 * (for time reversible models, there will be a node with zero parents.)
		 * Has no effect beyond catching bugs.
		 */
		public void checkConnections() {
			int count = 0;
			for (ConnectionOrRoot c : connections_) {
				if (!c.isTopNode(this)) count++;
			}
			if (count > 1) {
				throw new RuntimeException("Internal node failed Connection check: number of parents > 1");
			}
		}
		
		public UNode createAlteredCopy(Connection attachmentPoint, Node newSubtree, Connection originalParentConnection, Connection parentConnection, ConstructionTool tool) {
			return new InternalNode(this, attachmentPoint,newSubtree,originalParentConnection,parentConnection,tool);
		}
		public UNode createAlteredCopy( Connection originalParentConnection, Connection parentConnection, ConstructionTool tool) {
			return new InternalNode(this, originalParentConnection,parentConnection,tool);
		}

		public ModelNodeParameters getNodeParameters() {
			return nodeParams_;
		}
		public void setNodeParameters(ModelNodeParameters params) {
			nodeParams_ = params;
		}
		// Which of my connections is 'caller'?
		public int getBackpointer(ConnectionOrRoot caller) {
			int pointer = -1;
			for (int i=0; i<connections_.length; i++) {
				if (caller == connections_.clone()[i]) {
					pointer = i;
				}
			}
			if (pointer == -1 && caller != null) {
				throw new RuntimeException("Unknown caller!");
			}
			return pointer;
		}
		public ModelNodeParameters[] getAdjacentNodeParameters() {
			ModelNodeParameters[] mnps = new ModelNodeParameters[connections_.length];
			for (int i=0; i<connections_.length; i++) {
				UNode neighbour = connections_[i].getOther(this); // will be null if connections_[i] is Root
				if (neighbour != null) {
					mnps[i] = neighbour.getNodeParameters();
				}
			}
			return mnps;
		}
		public int[] getAdjacentNodeBackpointers() {
			int[] pointers = new int[connections_.length];
			for (int i=0; i<connections_.length; i++) {
				UNode neighbour = connections_[i].getOther(this); // will be null if connections_[i] is Root
				if (neighbour != null) {
					pointers[i] = neighbour.getBackpointer(connections_[i]);
				} else {
					pointers[i] = -1;
				}
			}
			return pointers;
		}
		public ModelEdgeParameters[] getAdjacentEdgeParameters() {
			ModelEdgeParameters[] meps = new ModelEdgeParameters[connections_.length];
			for (int i=0; i<connections_.length; i++) {
				meps[i] = connections_[i].getEdgeParameters();
			}
			return meps;
		}
		public int[] getAdjacentEdgeBackpointers() {
			int[] pointers = new int[connections_.length];
			for (int i=0; i<connections_.length; i++) {
				if (connections_[i].getBottom() == this) {
					pointers[i] = 0;
				} else {
					if (connections_[i].getTop() == this) {
						pointers[i] = 1;
					} else {
						throw new RuntimeException("Inconsistent topology!");
					}
				}
			}
			return pointers;
		}

		public boolean fixConnectionDirections(ConnectionOrRoot caller) {
			boolean subtreeContainsRoot = false;
			for (ConnectionOrRoot c : connections_) {
				// fix directions in each subtree, count how many of the subtrees contain the root.
				if (c != caller) {
					if (c.fixConnectionDirections(this)) {
						// sanity check:
						if (subtreeContainsRoot) { throw new RuntimeException("Tree contains two roots"); }
						subtreeContainsRoot = true;
					}
				}
			}
			return subtreeContainsRoot;
		}
		public boolean isLeaf() { return false; }
		public boolean hasLabel(String label) { return false; }

		public void mark() {
			markConnections_[0] = connections_[0];
			markConnections_[1] = connections_[1];
			markConnections_[2] = connections_[2];
			if (nodeParams_ != null) {
				nodeParams_.mark();
			}
		}
		public void setAnnotation(Object annotation) {		  this.annotation_ = annotation;		}
		public Object getAnnotation() { return annotation_; }
		public String getLabel() { return null; }

		public void undoToMark() {
			connections_[0] = markConnections_[0];
			connections_[1] = markConnections_[1];
			connections_[2] = markConnections_[2];
			if (nodeParams_ != null) {
				nodeParams_.undoToMark();
			}
			topologyChanged();
		}
		private final void topologyChanged() {
		  this.topologyChangedSinceLastFlat_ = true;
		  this.topologyChangedSincleLastExtended_ = true;
		}
		public boolean hasDirectConnection(Connection c) {
			for(int i = 0 ; i < connections_.length ; i++) {
					if(connections_[i]==c) {	return true;	}
				}
				return false;

		}
		/**
		 * Can we get from this node to 'c' without going through 'caller'?
		 */
		public boolean hasConnection(Connection c, ConnectionOrRoot caller) {
			for(int i = 0 ; i < connections_.length ; i++) {
				if((connections_[i]==c)||(connections_[i]!=caller&&connections_[i].hasConnection(c,this))) {
					return true;
				}
			}
			return false;
		}

		public final int getIndex() { return index_; }
		public void testLikelihood(ConnectionOrRoot caller, SubstitutionModel model, ConstructionTool tool) {
			for(int i = 0 ; i < connections_.length ; i++) {
				if(connections_[i]!=caller) {
					connections_[i].testLikelihood(this,model,tool);
				}
			}
		}
		public void setConnections(Connection[] store, int number){
			if(number!=3) {
				throw new IllegalArgumentException("Must be three connections not:"+number);
			}
			System.arraycopy(store,0,connections_,0,3);
			topologyChanged();
		}

		public void getLeafNames(ArrayList store, Connection caller) {
			for(int i = 0 ; i < connections_.length ; i++) {
				if(connections_[i]!=caller) {
					connections_[i].getLeafNames(store,this);
				}
			}
		}

		public void getSplitInformation(int[] splitStore, String[] leafNames, int splitIndex, Connection caller) {
			for(int i = 0 ; i < connections_.length ; i++) {
				if(connections_[i]!=caller) {
					connections_[i].getSplitInformation(splitStore,leafNames,splitIndex,this);
				}
			}
		}

		//Interchange related
		public ConnectionOrRoot getLeft(Connection caller) {
			return connections_[LEFT_LOOKUP[getCallerIndex(caller)]];
		}

		public ConnectionOrRoot getRight(Connection caller) {
			return connections_[RIGHT_LOOKUP[getCallerIndex(caller)]];
		}

		/**
		 * Extract the given node from its current position, leaving it still attached to
		 * 'caller'. The node will come away with one 'dangling' Connection still attached.
		 * The return value is this dangling Connection. The node will have one bad pointer
		 * (the Connection not 'caller' and not the returned one) and the dangling Connection
		 * will have one bad pointer (the end not attached to 'this'.)
		 */
		/*
		 * Further notes
		 * Notation: [...] = an edge, (...) = a node.
		 * The node has three connections, which I will call [caller], [toStay], [toMove]. 
		 * Once we extract the node, [toStay] will remain attaching the nodes which used
		 * to be at the far ends of [toStay] and [toMove], i.e. we transform
		 * (nodeA) <-> [toStay] <-> (this) <-> [toMove] <-> (nodeB)
		 * into
		 * (nodeA) <-> [toStay] <-> (nodeB)
		 * and also end up with [caller] <-> (this) <-> [toMove]
		 * but the other links on (this) and [toMove] are left in an inconsistent state.
		 */
		public Connection extract(Connection caller) {
			int callerIndex = getCallerIndex(caller);
			Connection toStay = (Connection)connections_[LEFT_LOOKUP[callerIndex]];
			Connection toMove = (Connection)connections_[RIGHT_LOOKUP[callerIndex]];
			UNode nodeB = toMove.getOther(this);
			toStay.swapNode(this,nodeB);
			nodeB.swapConnection(toMove,toStay);
			topologyChanged();
			return toMove;
		}

		
		public void swapConnection(Connection original, UNode nodeToReplace, Connection newConnection) {
			int index = getCallerIndex(original);
			connections_[index] = newConnection;
			newConnection.swapNode(nodeToReplace,this);
			original.swapNode(this,nodeToReplace);

			nodeToReplace.swapConnection(newConnection,original);
			topologyChanged();
		}
		public void swapConnection(Connection original,Connection newConnection) {
			int index = getCallerIndex(original);
			connections_[index] = newConnection;
			topologyChanged();
		}
		public PatternInfo getPatternInfo(ConnectionOrRoot caller) {
			return patterns_[getCallerIndex(caller)];
		}

		// MDW: Recursive routine for changing a 'shadow tree' back into a standard PAL Node tree.
		public Node buildPALNode(double branchLength, ConnectionOrRoot caller) {
			final int callerIndex = getCallerIndex(caller);
			final Connection leftConnection = (Connection)connections_[LEFT_LOOKUP[callerIndex]];
			final Connection rightConnection = (Connection)connections_[RIGHT_LOOKUP[callerIndex]];
			Node[] children = new Node[] {
				leftConnection.buildPALNode(this), rightConnection.buildPALNode(this)
			};
			SimpleNode newNode = NodeFactory.createNodeBranchLength(branchLength,children);
			if (caller instanceof Connection) {
				ModelEdgeParameters edgeParams = ((Connection)caller).edgeParams_;
				if (edgeParams != null) {
					newNode.setAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL, edgeParams);
				}
			}
			if (nodeParams_ != null) {
				newNode.setAttribute(ModelNodeParameters.ATTRIBUTE_LABEL, nodeParams_);
			}
			return newNode;
			
		}

		public void instruct(UnrootedTreeInterface.UNode node, Connection callingConnection){
			// MDW: I have no idea what 'instruct' does, there is no documentation. If you need it,
			// uncomment the code and fix it yourself.
			throw new RuntimeException("Disabled because I didn't know how to update it sensibly");
			/*
			if(annotation_!=null) {		  node.setAnnotation(annotation_);			}

			final int callerIndex = getCallerIndex(callingConnection);
			final Connection leftConnection = connections_[LEFT_LOOKUP[callerIndex]];
			final Connection rightConnection = connections_[RIGHT_LOOKUP[callerIndex]];
			final UnrootedTreeInterface.UNode leftNode = node.createUChild();
			final UnrootedTreeInterface.UNode rightNode = node.createUChild();
		  leftConnection.instruct(leftNode.getParentUBranch(),this);
			rightConnection.instruct(rightNode.getParentUBranch(),this);
			*/
		}

		public String toString(Connection caller) {
			StringBuffer sb = new StringBuffer();
			boolean printed = false;
			for(int i = 0 ; i < connections_.length ; i++) {
				if(connections_[i]!=caller) {
					if(printed) {
						sb.append(", ");
					}
					printed = true;
					sb.append(connections_[i].toString(this));
				}
			}
			return sb.toString();
		}
		// --=-=-==--=
		public int rebuildPattern(ConstructionTool tool) {
			return rebuildPattern(tool, null,true);
		}
		/**
		 * We can handle being redirected to use, so we just call rebuild pattern
		 */
		public int redirectRebuildPattern(ConstructionTool tool) {
			return rebuildPattern(tool);
		}
		public void getAllConnections(ArrayList<Connection> store, ConnectionOrRoot caller) {
			for(int i = 0 ; i < connections_.length ; i++) {
				if(connections_[i]!=caller) {
					connections_[i].getAllConnections(store,this);
				}
			}
		}
		private final int getCallerIndex(ConnectionOrRoot caller) {
			if(caller==null) {
				throw new IllegalArgumentException("getCallerIndex() called on null object");
			}
			if(caller==connections_[0]) { return 0; }
			if(caller==connections_[1]) { return 1; }
			if(caller==connections_[2]) { return 2; }
			throw new IllegalArgumentException("Unknown caller");
		}


		public int rebuildPattern(ConstructionTool tool, Connection caller,boolean firstPass) {
			if(firstPass) {
				//First pass
				ModelNodeParameters exampleNodeParams = nodeParams_;
				int callerIndex = -1;
				for(int i = 0 ; i < connections_.length ; i++) {
					if(connections_[i]==caller) {
						callerIndex = i;
					} else if (connections_[i] instanceof Connection) {
						// do nothing if connections_[i] is Root instead. 
						UNode other = connections_[i].getOther(this);
						other.rebuildPattern(tool,(Connection)connections_[i],firstPass);
						if (exampleNodeParams == null) { exampleNodeParams = other.getNodeParameters(); }
						if (connections_[i].getEdgeParameters() != null) {
							// As of time of writing, this does nothing on first pass - this is future-proofing.
							connections_[i].getEdgeParameters().update(this.nodeParams_, i, other.getNodeParameters(), other.getBackpointer(connections_[i]), firstPass);
						}
					}
				}
				if (exampleNodeParams != null) {
					// This tree uses ModelNodeParameters so do updates. This process is a close parallel to rebuilding patterns.
					// As we've 
					// Give this node a MNP if it hasn't got one already, copying from the one found. (3 = number of Connections)
					if (nodeParams_ == null) { nodeParams_ = exampleNodeParams.getInstance(3); }
					nodeParams_.resetForTopologyChanged();
					nodeParams_.update(
								this.getAdjacentNodeParameters(),
								this.getAdjacentNodeBackpointers(),
								this.getAdjacentEdgeParameters(),
								this.getAdjacentEdgeBackpointers(),
								callerIndex,
								firstPass);
				}
				if(caller==null) {
					//Second pass
					return rebuildPattern(tool,null,false);
				} else {
					return rebuildMyPattern(tool, getCallerIndex(caller));
				}
			} else {
				//Second pass
				if (nodeParams_ != null) {
					nodeParams_.update(
							this.getAdjacentNodeParameters(),
							this.getAdjacentNodeBackpointers(),
							this.getAdjacentEdgeParameters(),
							this.getAdjacentEdgeBackpointers(),
							this.getBackpointer(caller),
							firstPass);
				}
				for(int i = 0 ; i < connections_.length ; i++) {
					if (connections_[i].getEdgeParameters() != null) {
						// If we have edge parameters, try to update them for the new tree topology
						UNode other = connections_[i].getOther(this);
						if (other != null) {
							connections_[i].getEdgeParameters().update(this.nodeParams_, i, other.getNodeParameters(), other.getBackpointer(connections_[i]), firstPass);
						} else {
							// if connections_[i] is Root, we'll do this
							connections_[i].getEdgeParameters().update(this.nodeParams_, i, null, 0, firstPass);
						}
					}
					if(connections_[i]!=caller) {
						rebuildMyPattern(tool,i);
					}
				}
				for(int i = 0 ; i < connections_.length ; i++) {
					if(connections_[i]!=caller && connections_[i] instanceof Connection) {
						connections_[i].getOther(this).rebuildPattern(tool,(Connection)connections_[i],firstPass);
					}
				}
				int max = 0;
				for(int i = 0; i < patterns_.length ; i++) {
					int count = patterns_[i].getNumberOfPatterns();
					if(count>max) { max = count; }
				}
				return max;
			}
		}

		public ConditionalProbabilityStore getLeftExtendedConditionalProbabilities( SubstitutionModel model, boolean modelChanged,  Connection callingConnection, LHCalculator.External external, ConditionalProbabilityStore resultStore) {
			final int callerIndex = getCallerIndex(callingConnection);
			final ConnectionOrRoot leftConnection = connections_[LEFT_LOOKUP[callerIndex]];
			return leftConnection.getExtendedConditionalProbabilities(model, modelChanged, this,external, resultStore);
		}
		public ConditionalProbabilityStore getRightExtendedConditionalProbabilities( SubstitutionModel model, boolean modelChanged,  Connection callingConnection, LHCalculator.External external, ConditionalProbabilityStore resultStore) {
			final int callerIndex = getCallerIndex(callingConnection);
			final ConnectionOrRoot rightConnection = connections_[RIGHT_LOOKUP[callerIndex]];
			return rightConnection.getExtendedConditionalProbabilities(model, modelChanged, this, external, resultStore);
		}
		public PatternInfo getLeftPatternInfo(Connection caller) {
			final int callerIndex = getCallerIndex(caller);
			final ConnectionOrRoot leftConnection = connections_[LEFT_LOOKUP[callerIndex]];
			return leftConnection.getFarPatternInfo(this);
		}

		public PatternInfo getRightPatternInfo(Connection caller) {
			final int callerIndex = getCallerIndex(caller);
			final ConnectionOrRoot rightConnection = connections_[RIGHT_LOOKUP[callerIndex]];
			return rightConnection.getFarPatternInfo(this);
		}


		public ConditionalProbabilityStore getFlatConditionalProbabilities( final SubstitutionModel model, final boolean modelChanged, final Connection callerConnection, LHCalculator.External externalCalculator, ConditionalProbabilityStore resultStore) {
			final int callerIndex = getCallerIndex(callerConnection);
			final PatternInfo pi = patterns_[callerIndex];
			final ConnectionOrRoot leftConnection = connections_[LEFT_LOOKUP[callerIndex]];
			final ConnectionOrRoot rightConnection = connections_[RIGHT_LOOKUP[callerIndex]];

			externalCalculator.calculateFlat(
					 pi,
					 leftConnection.getExtendedConditionalProbabilities(model, modelChanged, this,0),
					 rightConnection.getExtendedConditionalProbabilities(model, modelChanged, this,0),
					 resultStore
					 );
			return resultStore;

		}
		public ConditionalProbabilityStore getFlatConditionalProbabilities( final SubstitutionModel model, final boolean modelChanged, final ConnectionOrRoot callerConnection, int depth) {
			final int callerIndex = getCallerIndex(callerConnection);
			final PatternInfo pi = patterns_[callerIndex];
			final ConnectionOrRoot leftConnection = connections_[LEFT_LOOKUP[callerIndex]];
			final ConnectionOrRoot rightConnection = connections_[RIGHT_LOOKUP[callerIndex]];
			topologyChangedSinceLastFlat_ = false;
			return calculator_.calculateFlat( pi,
					 leftConnection.getExtendedConditionalProbabilities(model, modelChanged, this,depth+1),
					 rightConnection.getExtendedConditionalProbabilities(model, modelChanged, this,depth+1)
					 );
		}
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( 
				final double distance,
				final SubstitutionModel model, 
				final boolean modelChanged,  
				final Connection callerConnection, 
				LHCalculator.External externalCalculator, 
				ConditionalProbabilityStore resultStore
		) {
			checkConnections(); // sanity check, has no effect unless there is a bug.
			
			final int callerIndex = getCallerIndex(callerConnection);
			final PatternInfo pi = patterns_[callerIndex];
			final ConnectionOrRoot leftConnection = connections_[LEFT_LOOKUP[callerIndex]];
			final ConnectionOrRoot rightConnection = connections_[RIGHT_LOOKUP[callerIndex]];
			boolean forwardsInTime = callerConnection.isTopNode(this);
			externalCalculator.calculateExtended(
				distance, forwardsInTime, callerConnection.edgeParams_, model, pi,
				leftConnection.getExtendedConditionalProbabilities(model, modelChanged, this,0),
				rightConnection.getExtendedConditionalProbabilities(model, modelChanged, this,0), 
				resultStore
				);
			return resultStore;
		}
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( 
				final double distance,
				final SubstitutionModel model, 
				final boolean modelChanged,  
				final Connection callerConnection, 
				int depth
		) {
			checkConnections(); // sanity check, has no effect unless there is a bug.
			
			final int callerIndex = getCallerIndex(callerConnection);
			final PatternInfo pi = patterns_[callerIndex];
			final ConnectionOrRoot leftConnection = connections_[LEFT_LOOKUP[callerIndex]];
			final ConnectionOrRoot rightConnection = connections_[RIGHT_LOOKUP[callerIndex]];
			topologyChangedSincleLastExtended_ = false;
			boolean forwardsInTime = callerConnection.isTopNode(this);
			return calculator_.calculateExtended(
				distance,
				forwardsInTime,
				callerConnection.edgeParams_, 
				model,
				pi,
				leftConnection.getExtendedConditionalProbabilities(model, modelChanged, this,depth+1),
				rightConnection.getExtendedConditionalProbabilities(model, modelChanged, this,depth+1),   
				modelChanged
				);
		}

		// -=-==--=
		private int rebuildMyPattern(ConstructionTool tool, int index) {
			ConnectionOrRoot  leftConnection = connections_[ LEFT_LOOKUP[index]];
			ConnectionOrRoot rightConnection = connections_[RIGHT_LOOKUP[index]];
			final PatternInfo  leftPattern =  leftConnection.getFarPatternInfo(this);
			final PatternInfo rightPattern = rightConnection.getFarPatternInfo(this);
			return tool.build(patterns_[index],leftPattern,rightPattern);
		} //End of rebuildPatternWeights()
	}

// -=-=-=-=-=-=-==--=-=-=-=-=-=-==--=-=-==-=--=-=-=-==-=--=-=-=-=-=-=-==-=-=-=-
// === LeafNode
// -=-=-=-=-=-=-==--=-=-=-=-=-=-==--=-=-==-=--=-=-=-==-=--=-=-=-=-=-=-==-=-=-=-

	private static final class LeafNode implements UNode {
		private final String id_;
		private Connection parentConnection_;
		private Connection markParentConnection_;
		private final int[] sequence_;
		private final PatternInfo pattern_;
		private final int index_;
		private ModelNodeParameters nodeParams_; 

		private final LHCalculator.Leaf leafCalculator_;

		private Object annotation_;

		public LeafNode(String id, Connection parentConnection, ConstructionTool tool ) {
			this.id_ = id;
			this.parentConnection_ = parentConnection;
			this.index_ = tool.allocateNextUNodeIndex(this);
			this.sequence_ = tool.getSequence(id_);
			if(sequence_==null) {
				this.leafCalculator_ = null; this.pattern_ = null;
			} else {
				final int numberOfSites = sequence_.length;
				final int numberOfStates = tool.getNumberOfStates();

				final int[] patternStateMatchup = new int[numberOfStates+1];
				final int[] sitePatternMatchup = new int[numberOfSites];

				final int uniqueCount = createMatchups( numberOfSites, numberOfStates, sitePatternMatchup, patternStateMatchup );

				this.leafCalculator_ = tool.createNewLeafCalculator( patternStateMatchup, uniqueCount );
				this.pattern_ = new PatternInfo( sitePatternMatchup, uniqueCount );
			}
		}
		public void getLeafNames(ArrayList store, Connection caller) {
		  if(caller!=parentConnection_) {
				throw new RuntimeException("Unknown caller!");
			}
			store.add(id_);
		}
		public boolean isLeaf() { return true; }
		public boolean hasLabel(String label) { return id_.equals(label); }

		public void getSplitInformation(int[] splitStore, String[] leafNames, int splitIndex, Connection caller) {
			if(caller!=parentConnection_) {
				throw new RuntimeException("Unknown caller!");
			}
			for(int i = 0 ; i < leafNames.length ;i++) {
				if(id_.equals(leafNames[i])) {
					splitStore[i] = splitIndex;
				}
			}
		}
		public void setNodeParameters(ModelNodeParameters mnp) {
			nodeParams_ = mnp;
		}
		public ModelNodeParameters getNodeParameters() {
			return nodeParams_;
		}
		public int getBackpointer(ConnectionOrRoot caller) {
			if (caller != parentConnection_) {
				throw new RuntimeException("Unknown caller!");
			}
			return 0;
		}
		public ModelNodeParameters[] getAdjacentNodeParameters() {
			return new ModelNodeParameters[] {parentConnection_.getOther(this).getNodeParameters()};
		}
		public int[] getAdjacentNodeBackpointers() {
			return new int[] {parentConnection_.getOther(this).getBackpointer(parentConnection_)};
		}
		public ModelEdgeParameters[] getAdjacentEdgeParameters() {
			return new ModelEdgeParameters[] {parentConnection_.getEdgeParameters()};
		}
		public int[] getAdjacentEdgeBackpointers() {
			if (parentConnection_.getBottom() == this) {
				return new int[] {0};
			} 
			if (parentConnection_.getTop() == this) {
				return new int[] {1};
			}
			throw new RuntimeException("Inconsistent topology!");
		}

		public void setAnnotation(Object annotation) {		  this.annotation_ = annotation;		}
		public Object getAnnotation() { return annotation_; }
		public String getLabel() { return id_; }
		/**
		 * Fill in matchup arrahs
		 * @param numberOfSites The number of sites
		 * @param sitePatternMatchup Should be of length numberOfSites
		 * @param patternStateMatchup Should be of length numberOfStates+1
		 * @return
		 */
		private final int createMatchups(final int numberOfSites, final int numberOfStates, final int[] sitePatternMatchup, final int[] patternStateMatchup) {
		  final int[] stateCount  = new int[numberOfStates+1];
			// StatePatternMatchup matches a state to it's new pattern (is undefined if state does not occur)
			final int[] statePatternMatchup = new int[numberOfStates+1];
			int uniqueCount = 0;
			for(int site = 0 ; site < numberOfSites ; site++) {
				final int state = sequence_[site];
				if(stateCount[state]==0) {
					stateCount[state] = 1;
					int pattern = uniqueCount++;
					patternStateMatchup[pattern] = state;
					statePatternMatchup[state] = pattern;
				} else {
					stateCount[state]++;
				}
				sitePatternMatchup[site] = statePatternMatchup[state];
			}
			return uniqueCount;
		}

		// MDW: Allow ModelEdgeParameters and ModelNodeParameters to be transfered from 
		// a tree of Nodes: if the Node has an MEP, attach it to the parent connection.
		// If it has an MNP, attach it to this node.
		public LeafNode(Node c, Connection parentConnection, ConstructionTool tool ) {
			this(c.getIdentifier().getName(),parentConnection,tool);
			if (c instanceof AttributeNode) {
				parentConnection.edgeParams_ = (ModelEdgeParameters) ((AttributeNode) c).getAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL);
				nodeParams_                  = (ModelNodeParameters) ((AttributeNode) c).getAttribute(ModelNodeParameters.ATTRIBUTE_LABEL);
			}
		}
		private LeafNode(LeafNode base, Connection parentConnection, ConstructionTool tool ) {
			this.id_ = base.id_;
			this.parentConnection_ = parentConnection;
			this.index_ = tool.allocateNextUNodeIndex(this);
			this.sequence_ = tool.getSequence(id_);

			final int numberOfSites = sequence_.length;
			final int numberOfStates = tool.getNumberOfStates();

			final int[] patternStateMatchup = new int[numberOfStates+1];
			final int[] sitePatternMatchup = new int[numberOfSites];

			final int uniqueCount = createMatchups(numberOfSites,numberOfStates,sitePatternMatchup,patternStateMatchup);

			this.leafCalculator_ = tool.createNewLeafCalculator(patternStateMatchup,uniqueCount);
			this.pattern_ = new PatternInfo(sitePatternMatchup,uniqueCount);
		}


		public UNode createAlteredCopy(Connection attachmentPoint, Node newSubtree, Connection originalParentConnection, Connection parentConnection, ConstructionTool tool) {
			return new LeafNode(this, parentConnection,	tool);
		}
		public UNode createAlteredCopy(Connection originalParentConnection, Connection parentConnection, ConstructionTool tool) {
			return new LeafNode(this, parentConnection,	tool);
		}
		public void instruct(UnrootedTreeInterface.UNode node, Connection callingConnection){
		  node.setLabel(id_);
			if(annotation_!=null) {		  node.setAnnotation(annotation_);			}
			if(callingConnection!=parentConnection_) {
			  throw new IllegalArgumentException("Unknown calling connection!");
			}
		}

		public boolean hasDirectConnection(Connection c) {	return parentConnection_==c;	}
		public void mark() {
			this.markParentConnection_ = parentConnection_;
			if (nodeParams_ != null) {
				nodeParams_.mark();
			}
		}
		public void undoToMark() {	
			this.parentConnection_ = markParentConnection_;
			if (nodeParams_ != null) {
				nodeParams_.undoToMark();
			}
		}
		public boolean hasConnection(Connection c, ConnectionOrRoot caller) {
			if(caller!=parentConnection_) {		throw new IllegalArgumentException("Unknown caller!");	}
			return parentConnection_==c;
		}
		public boolean fixConnectionDirections(ConnectionOrRoot caller) { return false; }
		public Connection extract(Connection caller) {
			if(caller!=parentConnection_) {		throw new IllegalArgumentException("Unknown caller!");	}
			return null;
		}

		public ConditionalProbabilityStore getLeftExtendedConditionalProbabilities( SubstitutionModel model, boolean modelChanged,  Connection callingConnection, LHCalculator.External externalCalculator, ConditionalProbabilityStore resultStore){
			throw new RuntimeException("Assertion error : Not applicable for leaf nodes!");
		}
		public ConditionalProbabilityStore getRightExtendedConditionalProbabilities( SubstitutionModel model, boolean modelChanged,  Connection callingConnection, LHCalculator.External externalCalculator, ConditionalProbabilityStore resultStore){
			throw new RuntimeException("Assertion error : Not applicable for leaf nodes!");
		}
		public PatternInfo getLeftPatternInfo(Connection caller){		return null;	}
		public PatternInfo getRightPatternInfo(Connection caller) {	return null;	}

		public void setConnections(Connection[] store, int number){
			if(number!=1) {		throw new IllegalArgumentException("Must be one connection not:"+number);		}
			this.parentConnection_ = store[0];
		}

		public void testLikelihood(ConnectionOrRoot caller, SubstitutionModel model, ConstructionTool tool) {
			if(caller!=parentConnection_) {	throw new IllegalArgumentException("Unknown caller!");		}
		}
		public final int getIndex() { return index_; }
		public void swapConnection(Connection original,  Connection newConnection) {
			if(original!=parentConnection_) {		throw new IllegalArgumentException("Unknown original");		}
			this.parentConnection_ = newConnection;
		}
		public void swapConnection(Connection original, UNode nodeToReplace, Connection newConnection) {
			swapConnection(original,newConnection);
			newConnection.swapNode(nodeToReplace,this);
			original.swapNode(this,nodeToReplace);
		}

		/**
		 * @return null (as not possible)
		 */
		public Connection getLeft(Connection caller) {	return null;	}
		/**
		 * @return null (as not possible)
		 */
		public Connection getRight(Connection caller) {	return null;	}

		public void getAllConnections(ArrayList store, ConnectionOrRoot caller) {
			if(caller!=parentConnection_) {		throw new IllegalArgumentException("Unknown caller!");		}
		}
		public PatternInfo getPatternInfo(ConnectionOrRoot caller){
			if(caller!=parentConnection_) {		throw new IllegalArgumentException("Unknown caller!");		}
			return pattern_;
		}
		public void rebuildConnectionPatterns(ConstructionTool tool, Connection caller) {
			if(caller!=parentConnection_){	throw new IllegalArgumentException("Unknown caller!");			}
		}

		public int rebuildPattern( ConstructionTool tool, Connection caller, boolean firstPass) {
			if(caller!=parentConnection_) {		throw new IllegalArgumentException("Uknown caller!");		}
			if (nodeParams_ != null) {
				nodeParams_.update(
							this.getAdjacentNodeParameters(),
							this.getAdjacentNodeBackpointers(),
							this.getAdjacentEdgeParameters(),
							this.getAdjacentEdgeBackpointers(),
							0,
							firstPass);
			}
			return pattern_.getNumberOfPatterns();
		}
		public int rebuildPattern(ConstructionTool tool) {	return parentConnection_.getOther(this).redirectRebuildPattern(tool);		}

		/**
		 * This should only be called by another leaf node on the other end of the connection.
		 * In this case we don't have to do much (tree is two node tree)
		 */
		public int redirectRebuildPattern(ConstructionTool tool) {		return pattern_.getNumberOfPatterns();		}

		public final ConditionalProbabilityStore getFlatConditionalProbabilities(SubstitutionModel model, boolean modelChanged, Connection callingConnection, LHCalculator.External external, ConditionalProbabilityStore resultStore) {
			if(callingConnection!=parentConnection_) {		throw new IllegalArgumentException("Unknown calling connection");			}
			return leafCalculator_.getFlatConditionalProbabilities();
		}
		public final ConditionalProbabilityStore getFlatConditionalProbabilities(final SubstitutionModel model, final boolean modelChanged, final ConnectionOrRoot callingConnection, int depth) {
			if(callingConnection!=parentConnection_) {		throw new IllegalArgumentException("Unknown calling connection");			}
			return leafCalculator_.getFlatConditionalProbabilities();
		}

		public ConditionalProbabilityStore getExtendedConditionalProbabilities( double distance, SubstitutionModel model, boolean modelChanged, Connection callingConnection, LHCalculator.External external, ConditionalProbabilityStore resultStore) {
		  if(callingConnection!=parentConnection_) {	throw new IllegalArgumentException("Unknown calling connection");		}
		  	// MDW edgeParams_ added to support MosaicSubstitutionModel
			return leafCalculator_.getExtendedConditionalProbabilities(distance,callingConnection.edgeParams_,model,modelChanged); 			   
		}
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( double distance, SubstitutionModel model, boolean modelChanged, Connection callingConnection, int depth) {
			if(callingConnection!=parentConnection_) {	throw new IllegalArgumentException("Unknown calling connection");			}
		  	// MDW edgeParams_ added to support MosaicSubstitutionModel
			return leafCalculator_.getExtendedConditionalProbabilities(distance,callingConnection.edgeParams_,model,modelChanged); 			   
		}
		public Node buildPALNode(double branchLength, ConnectionOrRoot caller) {
			if(caller!=parentConnection_) {		throw new IllegalArgumentException("Unknown calling connection"); 	}
			SimpleNode newNode = NodeFactory.createNodeBranchLength(branchLength, new Identifier(id_));
			if (caller instanceof Connection) {
				ModelEdgeParameters edgeParams = ((Connection)caller).edgeParams_;
				if (edgeParams != null) {
					newNode.setAttribute(ModelEdgeParameters.ATTRIBUTE_LABEL, edgeParams);
				}
			}
			if (nodeParams_ != null) {
				newNode.setAttribute(ModelNodeParameters.ATTRIBUTE_LABEL, nodeParams_);
			}

			return newNode;
		}
		public String toString(Connection caller) {		return id_; 	}
	} //End of class Leaf

//

	private static final class BranchAccessImpl implements BranchAccess {
	  private final Connection peer_;
		private final PossiblyRootedMLSearcher base_;
		public BranchAccessImpl(Connection peer, PossiblyRootedMLSearcher base) {
		  this.peer_ = peer;
			this.base_ = base;
		}
		public boolean isLeafBranch(String leafLabel) {
		  return peer_.isLeafBranch(leafLabel);
		}
		public PossiblyRootedMLSearcher attach(Node subTree, Alignment alignment) {
	  	return new PossiblyRootedMLSearcher(base_, peer_, subTree, alignment,base_.model_ );
	  }
		public PossiblyRootedMLSearcher attach(String sequence, Alignment alignment) {
	  	return new PossiblyRootedMLSearcher(base_, peer_, NodeFactory.createNode(new Identifier(sequence)), alignment, base_.model_ );
	  }
		public PossiblyRootedMLSearcher attach(Node subTree, Alignment alignment, SubstitutionModel model) {
	  	return new PossiblyRootedMLSearcher(base_, peer_, subTree, alignment, model );
	  }
		public PossiblyRootedMLSearcher attach(String sequence, Alignment alignment, SubstitutionModel model) {
	  	return new PossiblyRootedMLSearcher(base_, peer_, NodeFactory.createNode(new Identifier(sequence)), alignment, model );
	  }
		public void setAnnotation(Object annotation) {
		  peer_.setAnnotation(annotation);
		}
		public Object getAnnotation() { return peer_.getAnnotation(); }
		public String[] getBottomLeafNames() {
			return peer_.getBottomLeafNames();
		}
		public String[] getTopLeafNames() {
		  return peer_.getTopLeafNames();
		}
		public int[] getSplitInformation(String[] leafNames) {
		  return peer_.getSplitInformation(leafNames);
		}
	}

	private static final class NodeAccessImpl implements NodeAccess {
	  private final UNode peer_;
		private final PossiblyRootedMLSearcher base_;
		public NodeAccessImpl(UNode peer, PossiblyRootedMLSearcher base) {
		  this.peer_ = peer;
			this.base_ = base;
		}
		public boolean isLeaf() {	  return peer_.isLeaf();		}
		public void setAnnotation(Object annotation) {		  peer_.setAnnotation(annotation);		}
		public Object getAnnotation() { return peer_.getAnnotation(); }
		public String getLabel() { return peer_.getLabel(); }
	}

	
// -=-=-=-=-=-==--==--==--=-=-==-=--=-=-=-==-=-=--==-=-=--==--==-=-=--==--==--=
// == ConnectionOrRoot ===
// =--=-==-=--=-==--=-==--=-=-=-=-==-=-=--=-==-=-=-=--==-=-=-=-=--==-=-=-=-=-=-
	private interface ConnectionOrRoot {
		public UNode getOther(UNode caller);
		public ModelEdgeParameters getEdgeParameters();
		public UNode getBottom();
		public UNode getTop();
		public boolean isTopNode(UNode caller);
		public boolean hasConnection(Connection c, UNode caller);
		/**
		 * Recurse through the tree evaluating likelihood at every Connection/Root.
		 */
		public void testLikelihood(SubstitutionModel model, ConstructionTool tool);
		public void testLikelihood(UNode caller, SubstitutionModel model, ConstructionTool tool);
		public void getLeafNames(ArrayList<String> store, UNode caller);
		public void getSplitInformation(int[] splitStore, String[] leafNames, int splitIndex, UNode caller);
		public String toString();
		public String toString(UNode caller);
		public void getAllConnections(ArrayList<Connection> store);
		public void getAllConnections(ArrayList<Connection> store, UNode caller);
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( 
				SubstitutionModel model, 
				boolean modelChanged, 
				UNode caller, 
				LHCalculator.External externalCalculator, 
				ConditionalProbabilityStore extendedStore);
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( 
				SubstitutionModel model, 
				boolean modelChanged, 
				UNode caller, 
				int depth);			
		public PatternInfo getFarPatternInfo(UNode caller);
		public void setup( ConstructionTool tool, Connection[] allConnections);
		public void instructBase(UnrootedTreeInterface.BaseBranch base);
		/*
		 * calculateLogLikelihood: pass flat conditionals plus branch length to the calculator.
		 * calculateLogLikelihood2: pass flat bottom conditional and extended top conditional to calculator
		 * calculateLogLikelihood3: pass flat top conditional and extended bottom conditional to calculator
		 * (for Root, they are identical.) Up to rounding, these should all produce the same result.
		 */
		public double calculateLogLikelihood(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool);
		public double calculateLogLikelihood2(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool);
		public double calculateLogLikelihood3(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool);
		public SiteDetails calculateSiteDetails(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool);
		public Node buildPALNode();
		/**
		 * Guarantees that subtree as seen from caller has all Connections correctly pointing towards the root. 
		 * @param caller
		 * @return true if subtree contains the root, false otherwise.
		 */
		public boolean fixConnectionDirections(UNode caller); 
		
	}
	
// -=-=-=-=-=-==--==--==--=-=-==-=--=-=-=-==-=-=--==-=-=--==--==-=-=--==--==--=
// == Root ===
// =--=-==-=--=-==--=-==--=-=-=-=-==-=-=--=-==-=-=-=--==-=-=-=-=--==-=-=-=-=-=-
	/*
	 * The 'Root' class is a simplifed version of 'Connection'. It attaches to an InternalNode
	 * like a Connection, and presents the root base frequencies as the likelihoods on this 'branch'.
	 */
	private static final class Root implements ConnectionOrRoot {
		private InternalNode bottomNode_;
		private PatternInfo topPattern_; // invariant: a single pattern.
		private PatternInfo centerPattern_; // combination of topPattern_ and whatever comes from adjacent node
		private boolean centerPatternValid_;
		private ConditionalProbabilityStore topConditionalProbabilities_; // just the root base frequencies
		private double[] rootFreq_;
		private ModelEdgeParameters edgeParams_; // not normally needed, used for MosaicSubstitutionModel.
		
		/*
		public Root(NonTimeReversibleSubstitutionModel model, InternalNode node, int numberOfSites) {
			this(model.getNumberOfTransitionCategories(), numberOfSites, model.getDataType().getNumStates());
			topNode_ = node;
		}
		*/
		/**
		 * Start point for creating the 'shadow tree' for rooted trees.
		 */
		public Root(Node root, ConstructionTool tool, NonTimeReversibleSubstitutionModel model) {
			this(tool, model);
			bottomNode_ = new InternalNode(root,this,tool);
		}
		public Root(String[] leafNames, ConstructionTool tool, NonTimeReversibleSubstitutionModel model, MersenneTwisterFast r) {
			this(tool, model);
			bottomNode_ = new InternalNode(leafNames,this,tool,r);
		}
		// for use by public constructors. 
		private Root(ConstructionTool tool, NonTimeReversibleSubstitutionModel model) {
			int numberOfCategories = tool.getNumberOfTransitionCategories();
			int numberOfSites      = tool.getNumberOfSites();
			int numberOfStates     = tool.getNumberOfStates();
			rootFreq_ = new double[numberOfStates];
			// All sites have the same pattern (pattern 0)
			int[] sitePatternMatchup = new int[numberOfSites]; // initialized to all zero, as required.
			// The pattern has weight of the number of sites
			int[] patternWeights = new int[] {numberOfSites};
			topPattern_ = new PatternInfo(sitePatternMatchup, patternWeights, 1);
			topConditionalProbabilities_ = new ConditionalProbabilityStore(numberOfCategories, numberOfStates);
			// array[category][pattern][state]
			double[][][] array = topConditionalProbabilities_.getConditionalProbabilityAccess(1, false);
			// make innermost level of array all point to rootFreq_, so changes to rootFreq_ will
			// automagically apply to all categories.
			for (int cat=0; cat<numberOfCategories; cat++) {
				array[cat][0]=rootFreq_;
			}
			this.centerPattern_ = new PatternInfo(numberOfSites,true);
			centerPatternValid_ = false;
			edgeParams_ = model.getRootModelEdgeParameters();
		}
		
		public UNode getOther(UNode caller)            { return null; }
		public ModelEdgeParameters getEdgeParameters() { return edgeParams_; }
		public UNode getBottom()                       { return bottomNode_;}
		public UNode getTop()                          { return null; }
		public boolean isTopNode(UNode caller)         { return false; }
		public boolean hasConnection(Connection c, UNode caller) { return false; } // there is nothing on the far side of Root.
		public void testLikelihood(SubstitutionModel model, ConstructionTool tool)               {
			System.out.printf("Likelihood: %f, %f, %f\n",
					calculateLogLikelihood(model, true,tool.allocateNewExternalCalculator(), tool),
					calculateLogLikelihood2(model, true,tool.allocateNewExternalCalculator(), tool),
					calculateLogLikelihood3(model, true,tool.allocateNewExternalCalculator(), tool)
					);
			bottomNode_.testLikelihood(this,model,tool);
		}
		public void testLikelihood(UNode caller, SubstitutionModel model, ConstructionTool tool) {
			System.out.println("Likleihood:"+calculateLogLikelihood(model, true,tool.allocateNewExternalCalculator(), tool));
		}
		public void getLeafNames(ArrayList<String> store, UNode caller)                          {}
		public void getSplitInformation(int[] splitStore, String[] leafNames, int splitIndex, UNode caller) {}
		public String toString()                       { return "root"; }
		public String toString(UNode caller)           { return "root"; }
		// getAllConnections called without a UNode is trying to get Connections for the entire tree.
		public void getAllConnections(ArrayList<Connection> store) {
			bottomNode_.getAllConnections(store,this);
		}
		// getAllConnections called from a UNode wants all the connections beyond this one, which is nothing.
		public void getAllConnections(ArrayList<Connection> store, UNode caller) {}
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( 
				SubstitutionModel model, 
				boolean modelChanged, 
				UNode caller, 
				LHCalculator.External externalCalculator, 
				ConditionalProbabilityStore extendedStore) 
		{
			if (modelChanged) {  updateRootFrequencies(model); }
			return topConditionalProbabilities_;
		}
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( 
				SubstitutionModel model, 
				boolean modelChanged, 
				UNode caller, 
				int depth)
		{
			if (modelChanged) {  updateRootFrequencies(model); }
			return topConditionalProbabilities_;
		}
		public PatternInfo getFarPatternInfo(UNode caller) { return topPattern_; }
		public void setup( ConstructionTool tool, Connection[] allConnections) {
			if(tool.hasSequences()) {
				bottomNode_.rebuildPattern( tool );
			}
			for(int i = 0 ; i < allConnections.length ; i++) {
				allConnections[i].centerPatternValid_ = false;
			}
			this.centerPatternValid_ = false;
		}
		public void instructBase(UnrootedTreeInterface.BaseBranch base) { throw new RuntimeException("Not implemented"); }
		public double calculateLogLikelihood(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool) {
			if (modelChanged) {  updateRootFrequencies(model); }
			PatternInfo pi = getCenterPatternInfo(tool);
			final ConditionalProbabilityStore bottomConditionalProbabilityProbabilties =
				bottomNode_.getFlatConditionalProbabilities(model, modelChanged, this,0);
			return calculator.calculateLogLikelihood(model, pi, bottomConditionalProbabilityProbabilties, topConditionalProbabilities_, edgeParams_);
		}
		// In this case, there really aren't two ways to do it, so just fall through to the main method
		public double calculateLogLikelihood2(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool) {
			return calculateLogLikelihood(model, modelChanged, calculator, tool);
		}
		public double calculateLogLikelihood3(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool) {
			return calculateLogLikelihood(model, modelChanged, calculator, tool);
		}
		public SiteDetails calculateSiteDetails(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool) {
			PatternInfo pi = getCenterPatternInfo(tool);
			final ConditionalProbabilityStore bottom = bottomNode_.getFlatConditionalProbabilities(model, modelChanged,  this,0);
		  return calculator.calculateSiteDetailsUnrooted(0, null, model, pi, bottom, topConditionalProbabilities_, tool.newConditionalProbabilityStore(false));
		}
		public final PatternInfo getCenterPatternInfo(ConstructionTool tool) {
			if(!centerPatternValid_) {
				tool.build(centerPattern_, bottomNode_.getPatternInfo(this), topPattern_);
				centerPatternValid_ = true;
			}
			return centerPattern_;
		}
		public Node buildPALNode() {
			return bottomNode_.buildPALNode(BranchLimits.MINARC,this);
		}
		private void updateRootFrequencies(SubstitutionModel model) {
			// topConditionalProbabilities_ automagically gets changed by updating rootFreq_, so 
			// now it will be up to date.
			((NonTimeReversibleSubstitutionModel)model).getRootFrequencies(rootFreq_);
		}
		public boolean fixConnectionDirections(UNode caller) { return true; }
	} // class Root
	
// -=-=-=-=-=-==--==--==--=-=-==-=--=-=-=-==-=-=--==-=-=--==--==-=-=--==--==--=
// == Connection ===
// =--=-==-=--=-==--=-==--=-=-=-=-==-=-=--=-==-=-=-=--==-=-=-=-=--==-=-=-=-=-=-	
	private static final class Connection implements ConnectionOrRoot {
		private UNode bottomNode_;
		private UNode topNode_;
		private double branchLength_;
		private final PatternInfo centerPattern_;
		private boolean centerPatternValid_;

		private final int index_;

		private UNode markBottomNode_ = null;
		private UNode marktopNode_ = null;
		private double markBranchLength_;

		private Object annotation_ = null; 
		private ModelEdgeParameters edgeParams_ = null; // Holds any edge-specific parameters except the edge length. Usually null.

		/**
		 * The random tree constructor, for the root node (recursion just started)
		 * @param leafNames the names of all the leafs to be placed in this tree
		 * @param tool an aid for construction
		 * @param r a random number generator to determine branching patterns
		 */
		public Connection(String[] leafNames, ConstructionTool tool, MersenneTwisterFast r) {
			this.index_ = tool.allocateNextConnectionIndex();
			String[][] split = SearcherUtils.split(leafNames,r);
			this.bottomNode_ = createUNode(split[0],this,tool,r);
			this.topNode_ = createUNode(split[1],this,tool,r);
			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
			this.branchLength_ = CONSTRUCTED_BRANCH_LENGTH;
		}

		/**
		 * A random tree constructor, into the recursion
		 * @param leafNames The names of leaves remaining to be created
		 * @param parent the parent UNode (from previous recursion)
		 * @param tool to aid in construction
		 * @param r for determining branching
		 */
		public Connection(String[] leafNames , UNode parent, ConstructionTool tool, MersenneTwisterFast r) {
			this.index_ = tool.allocateNextConnectionIndex();
			this.branchLength_ = CONSTRUCTED_BRANCH_LENGTH;
			this.topNode_ = parent;
			this.bottomNode_ = createUNode(leafNames,this, tool,r);
			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
		}
		/**
		 * The starting constructor for building from a given tree
		 * @param n The normal PAL node structure to base this tree on
		 * @param tool to aid in construction
		 */
		// MDW: When creating the 'shadow tree', this is where we start. A Connection with no 'parent' UNode.
		public Connection(Node n,  ConstructionTool tool) {
			if(n.getChildCount()!=2) {
				throw new IllegalArgumentException("Base tree must be bificating");
			}
			this.index_ = tool.allocateNextConnectionIndex();
			Node l = n.getChild(0);
			Node r = n.getChild(1);
			this.branchLength_ = l.getBranchLength()+r.getBranchLength();

			// MDW: The recursive calls:
			bottomNode_ = createUNode(l, this, tool);
			topNode_ = createUNode(r, this, tool);

			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
		}
		/**
		 * Continuing recurison constructor for a given tree
		 * @param n The PAL node structure to base sub tree on
		 * @param parent The parent node (sub tree in other direction)
		 * @param tool to aid in construction
		 */
		public Connection(Node n, UNode parent, ConstructionTool tool) {
			this.index_ = tool.allocateNextConnectionIndex();
			this.branchLength_ = n.getBranchLength();
			this.topNode_ = parent;
			this.bottomNode_ = createUNode(n,this, tool); // MDW: Recursive
			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
		}

		/**
		 * A generic constructor given two already defined top and bottom children
		 * @param bottom The bottom node
		 * @param top The top node
		 * @param branchLength The length of connection
		 * @param tool to aid in construction
		 */
		public Connection(UNode bottom, UNode top, double branchLength, ConstructionTool tool) {
			this.index_ = tool.allocateNextConnectionIndex();
			this.branchLength_ = branchLength;
			this.topNode_ = top;
			this.bottomNode_ = bottom;
			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
		}
		/**
		 *
		 * @param originalBottom
		 * @param originalBottomParentConnection
		 * @param top
		 * @param branchLength
		 * @param tool
		 */
		public Connection(UNode originalBottom, Connection originalBottomParentConnection, UNode top, double branchLength, ConstructionTool tool) {
			this.index_ = tool.allocateNextConnectionIndex();
			this.branchLength_ = branchLength;
			this.topNode_ = top;
			this.bottomNode_ = originalBottom.createAlteredCopy(originalBottomParentConnection,this,tool);
			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
		}
		/**
		 * An altered tree constructor, initial starting point for recursion
		 * @param original The connection to base tree on
		 * @param attachmentPoint The original connection that is the attachment point for the new sub tree
		 * @param newSubtree The new sub tree as a normal PAL node
		 * @param tool to aid in construction
		 */
		public Connection(Connection original, Connection attachmentPoint, Node newSubtree, ConstructionTool tool) {
			//Allocate index like normal (do not keep indexes of original)
			this.index_ = tool.allocateNextConnectionIndex();

			final UNode baseTopNode = original.topNode_.createAlteredCopy(attachmentPoint,newSubtree,original,this,tool);
			if(original==attachmentPoint) {
			  //If this connection is where the sub tree is attached we will need to create a new Internal node "in the middle"
				// to attach the sub tree
				//This requires the addition of two extra connections (we go to node from right, another connection from left,
				//  and a third from subtree). All three connections meet at internal node (which is now the left node of this
				//  connection).

				//First step is to create subtree node, with this connection as it's parent
				final double splitLength = original.getBranchLength()/2;
			  this.branchLength_ = splitLength;

				this.bottomNode_ = new InternalNode(this,newSubtree, original.bottomNode_, original, splitLength,tool);
				this.topNode_ = baseTopNode;
			} else {
			  this.topNode_ = baseTopNode;
			  this.bottomNode_ = original.bottomNode_.createAlteredCopy(attachmentPoint,newSubtree,original,this,tool);
			  this.branchLength_ = original.getBranchLength();
			}
			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
		}
		/**
		 * Altered copy constructor
		 *
		 * @param original The original connection this one is based on
		 * @param attachmentPoint The connection that is the attachment point for the new sub tree
		 * @param newSubtree The new subtree in normal PAL form (as not part of the original)
		 * @param originalParent The original parent of the original connection (directing recursion)
		 * @param newParent The newly created parent (following same recursion as for original stuff)
		 * @param tool The normal construction tool to help us out
		 */
		public Connection(Connection original, Connection attachmentPoint, Node newSubtree, UNode originalParent, UNode newParent, ConstructionTool tool) {
			this.index_ = tool.allocateNextConnectionIndex();
			this.branchLength_ = original.getBranchLength();
			final boolean parentIsBottom;
			if(originalParent==original.bottomNode_) {
				parentIsBottom = true;
			} else if(originalParent==original.topNode_) {
				parentIsBottom = false;
			} else {
			  throw new RuntimeException("Assertion error : orginalParent does not belong to original connection!");
			}
			if(original==attachmentPoint) {
			  //If this connection is where the sub tree is attached we will need to create a new Internal node "in the middle"
				// to attach the sub tree
				//This requires the addition of two extra connections (we go to node from right, another connection from left,
				//  and a third from subtree). All three connections meet at internal node (which is now the left node of this
				//  connection).

				//First step is to create subtree node, with this connection as it's parent
				final double splitLength = original.getBranchLength()/2;
			  this.branchLength_ = splitLength;
				if(parentIsBottom) {
		  		this.bottomNode_ = newParent;
					this.topNode_ = new InternalNode(this,newSubtree, original.topNode_, original, splitLength,tool);
				} else {
		  		this.topNode_ = newParent;
					this.bottomNode_ = new InternalNode(this,newSubtree, original.bottomNode_, original, splitLength,tool);
				}
			} else {
			  if(parentIsBottom) {
				  this.bottomNode_ = newParent;
			    this.topNode_ = original.topNode_.createAlteredCopy(attachmentPoint,newSubtree,original,this,tool);
				} else {
					//Parent is right
				  this.topNode_ = newParent;
			    this.bottomNode_ = original.bottomNode_.createAlteredCopy(attachmentPoint,newSubtree,original,this,tool);
				}
			}
			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
		}
		/**
		 * The cloning constructor
		 * @param original The original connection that is being cloned
		 * @param originalParent The parent (dictating recursion) of the original connection
		 * @param newParent The new parent to substitute the original parent in the copy
		 * @param tool A tool to help us
		 */
		public Connection(Connection original,  UNode originalParent, UNode newParent, ConstructionTool tool) {
			this.index_ = tool.allocateNextConnectionIndex();
			this.branchLength_ = original.getBranchLength();
			if(originalParent==original.bottomNode_) {
			  this.bottomNode_ = newParent;
			  this.topNode_ = original.topNode_.createAlteredCopy(original,this,tool);
			} else if(originalParent==original.topNode_) {
			  this.topNode_ = newParent;
			  this.bottomNode_ = original.bottomNode_.createAlteredCopy(original,this,tool);
			} else {
			  throw new RuntimeException("Assertion error : orginalParent does not belong to original connection!");
			}
			this.centerPattern_ = new PatternInfo(tool.getNumberOfSites(),true);
			this.centerPatternValid_ = false;
		}
		public void setAnnotation(Object annotation) {	  this.annotation_ = annotation;		}
		public Object getAnnotation() { return this.annotation_; }
		/**
		 *
		 * @return The "bottom" node of this connection.
		 */
		public final UNode getBottom() { return bottomNode_; }
		/**
		 *
		 * @return The "Top" node of this connection.
		 */
		public final UNode getTop() { return topNode_; }
		
		public final boolean isTopNode(UNode caller) {
			if (caller == topNode_)    {return true;}
			if (caller == bottomNode_) {return false;}
			throw new IllegalArgumentException("Calling node not found");
		}
		
		public boolean hasDirectConnection(Connection c) {
			return topNode_.hasDirectConnection(c) || bottomNode_.hasDirectConnection(c);
		}
		public boolean hasDirectConnectionToRoot() {
			return (topNode_.getLeft(this) instanceof Root || topNode_.getRight(this) instanceof Root);
		}
		
		private final String[] getLeafNames(UNode base) {
			ArrayList<String> al = new ArrayList<String>();
			base.getLeafNames(al, this);
			String[] result = new String[al.size()];
			al.toArray(result);
			return result;
		}
		public ModelEdgeParameters getEdgeParameters() {
			return edgeParams_;
		}
		public String[] getBottomLeafNames() {		return getLeafNames(bottomNode_);	}
		public String[] getTopLeafNames()    {		return getLeafNames(topNode_);		}
		public int[] getSplitInformation(String[] leafNames) {
		  int[] result = new int[leafNames.length];
			for(int i = 0 ; i < result.length ; i++) { result[i] = 0;		}
			bottomNode_.getSplitInformation(result,leafNames,-1,this);
			topNode_.getSplitInformation(result,leafNames,1,this);
			return result;
		}
		/**
		 * @return true if this Connection is adjacent to a LeafNode
		 */
		public boolean isPendantEdge() {
			return bottomNode_.isLeaf() || topNode_.isLeaf();
		}
		public boolean isLeafBranch(String leafLabel) {
			if(bottomNode_.isLeaf() ){
				if(bottomNode_.hasLabel(leafLabel) ) { return true; }
			}
			if(topNode_.isLeaf() ){
				if(topNode_.hasLabel(leafLabel) ) { return true; }
			}
			return false;
		}
		public void getLeafNames(ArrayList<String> store, UNode caller) {
			if(caller==bottomNode_) {
				topNode_.getLeafNames(store,this);
			} else if (caller==topNode_) {
				bottomNode_.getLeafNames(store,this);
			} else {
				throw new RuntimeException("Unknown caller!");
			}
		}

		public void getSplitInformation(int[] splitStore, String[] leafNames, int splitIndex, UNode caller) {
			if(caller!=bottomNode_) {
				topNode_.getSplitInformation(splitStore,leafNames,splitIndex,this);
			} else if (caller!=topNode_) {
				bottomNode_.getSplitInformation(splitStore,leafNames,splitIndex,this);
			} else {
				throw new RuntimeException("Unknown caller!");
			}
		}

		/**
		 * Mark this node, or in other words store information on top and bottom nodes and branch length for later retreival (via undoToMark())
		 */
		public final void mark() {
			this.markBranchLength_ = branchLength_;
			this.markBottomNode_ = bottomNode_;		this.marktopNode_ = topNode_;
			if (edgeParams_ != null) {
				edgeParams_.mark();
			}
		}
		
		/**
		 * 
		 * @return the pattern info object for the other node leading to this connection  
		 */
		public PatternInfo getFarPatternInfo(UNode caller) {
			if (caller == bottomNode_) {
				return topNode_.getPatternInfo(this);
			} else if (caller == topNode_) {
				return bottomNode_.getPatternInfo(this);
			} else {
				throw new IllegalArgumentException("Unknown caller");
			}
		}
		
		/**
		 * @return The pattern info object for the bottom node leading to this connection
		 */
		public final PatternInfo getBottomPatternInfo() {		return bottomNode_.getPatternInfo(this);	}
		/**
		 * @return The pattern info object for the right node leading to this connection
		 */
		public final PatternInfo getTopPatternInfo() {
			return topNode_.getPatternInfo(this);
		}
		/**
		 *
		 * @return The pattern info across this connection (for use if this connection is the "root" of the likelihood calculation)
		 */
		public final PatternInfo getCenterPatternInfo(ConstructionTool tool) {
			if(!centerPatternValid_) {
				tool.build(centerPattern_, getBottomPatternInfo(),getTopPatternInfo());
				centerPatternValid_ = true;
			}
			return centerPattern_;
		}
		public final void instructBase(UnrootedTreeInterface.BaseBranch base) {
		  base.setLength(this.branchLength_);
			if(annotation_!=null) { base.setAnnotation(annotation_); }
			bottomNode_.instruct(base.getLeftNode(),this);
			topNode_.instruct(base.getRightNode(),this);
		}
		public final void instruct(UnrootedTreeInterface.UBranch branch, UNode callingNode) {
		  branch.setLength(this.branchLength_);
			if(annotation_!=null) { branch.setAnnotation(annotation_); }
			if(callingNode==bottomNode_) {
			  topNode_.instruct(branch.getFartherNode(),this);
			} else if(callingNode==topNode_) {
			  bottomNode_.instruct(branch.getFartherNode(),this);
			} else {
			  throw new IllegalArgumentException("Unknown calling node!");
			}
		}
		public final void undoToMark() {
			if(markBottomNode_==null) {
				throw new RuntimeException("Assertion error : undo to mark when no mark made");
			}
			this.branchLength_ = markBranchLength_;
			this.bottomNode_ = markBottomNode_;
			this.topNode_ = marktopNode_;
			if (edgeParams_ != null) {
				edgeParams_.undoToMark();
			}
		}

		public String toString() {
			return "("+bottomNode_+", "+topNode_+")";
		}
		public boolean hasConnection(Connection c, UNode caller) {
			if(c==this) { return true; }
			if(caller==bottomNode_) {
				return topNode_.hasConnection(c,this);
			}
			if(caller==topNode_) {
				return bottomNode_.hasConnection(c,this);
			}
			throw new IllegalArgumentException("Unknown caller");
		}

		/**
		 * @return the "left" connection of the bottom node
		 */
		public ConnectionOrRoot getBottomLeftConnection() {
			return bottomNode_.getLeft(this);
		}
		/**
		 * @return the "right" connection of the bottom node
		 */
		public ConnectionOrRoot getBottomRightConnection() {
			return bottomNode_.getRight(this);
		}
		/**
		 * @return the "left" connection of the top node
		 */
		public ConnectionOrRoot getTopLeftConnection() {
			return topNode_.getLeft(this);
		}
		/**
		 * @return the "right" connection of the top node
		 */
		public ConnectionOrRoot getTopRightConnection() {
			return topNode_.getRight(this);
		}

		/**
		 * @return connection that by attaching to we would undo this operation, null if operation not successful
		 * 'this' = connection to be removed.
		 *
		 */
		public ConnectionOrRoot attachTo(Connection attachmentPoint, Connection[] store) {

			// 'used' is the node attached to 'this' which is closer to 'attachment Point'.
			final UNode used = (bottomNode_.hasConnection(attachmentPoint, this) ? bottomNode_ : topNode_ );
			if(used.hasDirectConnection(attachmentPoint)) {
				return null;
			}
			final Connection redundant = used.extract(this);
			final ConnectionOrRoot reattachment;
			final ConnectionOrRoot leftUsed = used.getLeft(this);
			final ConnectionOrRoot rightUsed = used.getRight(this);

			if(leftUsed==redundant) {
				reattachment = rightUsed;
			} else if(rightUsed == redundant) {
				reattachment = leftUsed;
			} else {
				throw new IllegalArgumentException("Assertion error");
			}
			if(redundant==null) {
				throw new RuntimeException("Assertion error : I should be able to extract from one of my nodes!");
			}
			
			// Now we shift one of Connection attachmentPoint's nodes onto the end of Connection redundant.
			// It doesn't matter topologically which one, but we need to maintain the invariant that
			// a leaf node is never detached from its Connection.
			UNode attachmentNodeToMove = (attachmentPoint.topNode_.isLeaf()) ?  attachmentPoint.bottomNode_ : attachmentPoint.topNode_;

			// Make end of Connection attachmentPoint now point to node 'used' (detaching attachmentNodeToMove)
			attachmentPoint.swapNode(attachmentNodeToMove,used);
			// and shift attachmentNodeToMove to the dangling edge 'redundant's spare end 
			redundant.swapNode(redundant.getOther(used),attachmentNodeToMove);
			// and finally tell attachmentNodeToMove where it is now connected.
			attachmentNodeToMove.swapConnection(attachmentPoint,redundant);

			//Fix up the used connections
			store[0] = this;
			store[1] = redundant;
			store[2] = attachmentPoint;
			used.setConnections(store,3);
			
			/*
			 * Fix direction of Connections throughout the tree.
			 * This is a brute force approach - we could figure out exactly which Connections need
			 * flipping. However, it can't be a purely local solution: it is possible that all the
			 * Connections between 'this' and 'attachmentPoint' will need to be flipped.
			 */
			this.fixConnectionDirections();

			return reattachment;
		}
		// MDW: Entry point for convert-shadow-tree-to-PAL-Node-Tree. Returns root node. 
		public Node buildPALNode() {
			Node[] children = new Node[] {
				bottomNode_.buildPALNode(branchLength_/2,this),
				topNode_.buildPALNode(branchLength_/2,this)
			};
			return NodeFactory.createNode(children);
		}
		// MDW: Normal recursive conversion from shadow tree to PAL Node tree routine
		public Node buildPALNode(UNode caller) {
			if(bottomNode_==caller) {
				return topNode_.buildPALNode(branchLength_,this);
			}
			if(topNode_==caller) {
				return bottomNode_.buildPALNode(branchLength_,this);
			}
			throw new IllegalArgumentException("Unknown caller!");
		}

		public boolean fixConnectionDirections(UNode caller) {
			UNode other;
			boolean callerIsBottom;
			if (caller == topNode_) {
				callerIsBottom = false;
				other = bottomNode_;
			} else {
				callerIsBottom = true;
				other = topNode_;
			}
			boolean rootIsBeyondOther = other.fixConnectionDirections(this);
			if (rootIsBeyondOther ^ callerIsBottom) {
				// this edge is the wrong way around, so swap it.
				UNode temp = topNode_;
				topNode_ = bottomNode_;
				bottomNode_ = temp;
				centerPatternValid_ = false;
			}
			return rootIsBeyondOther;
		}
		/**
		 * Fix connection directions (so 'top' always points towards root) throughout the tree
		 * by recursion starting from this Connection. (For unrooted trees, Connections will be 
		 * directed as if this Connection is the root.)
		 */
		public void fixConnectionDirections() {
			topNode_.fixConnectionDirections(this);
			if (bottomNode_.fixConnectionDirections(this)) {
				UNode temp = topNode_;
				topNode_ = bottomNode_;
				bottomNode_ = temp;
				centerPatternValid_ = false;
			}
		}
		/**
		 * @return -1 if null
		 */
		private final static int getIndex(Connection c) {
			if(c==null) { return -1;}
			return c.index_;
		}
		/**
		 * Fill in the index information of the two connected nodes
		 */
		public void fillInConnectionState(int[] array, int insertionIndex) {
			final int l = bottomNode_.getIndex();
			final int r = topNode_.getIndex();
			if(l<r) {
				array[insertionIndex++] = l;
				array[insertionIndex] = r;
			} else {
				array[insertionIndex++] = r;
				array[insertionIndex] = l;
			}
		}
		/**
		 * Does nothing to fix up tree structure
		 */
		public void setNodes(UNode bottom, UNode top) {
			this.bottomNode_ = bottom;
			this.topNode_ = top;
		}
		/**
		 * @note does not change the nodes connection information. Leaves tree in an inconsitent state
		 */
		public void swapNode(UNode nodeToReplace, UNode replacement) {
			if(nodeToReplace==bottomNode_) {
				bottomNode_ = replacement;
			} else if(nodeToReplace==topNode_) {
				topNode_ = replacement;
			} else {
				throw new RuntimeException("Unknown node to replace");
			}
		}
		public final ConditionalProbabilityStore getBottomFlatConditionalProbabilities(SubstitutionModel model, boolean modelChanged) {
			return bottomNode_.getFlatConditionalProbabilities(model,modelChanged, this,0);
		}
		public final ConditionalProbabilityStore getTopFlatConditionalProbabilities(SubstitutionModel model, boolean modelChanged) {
			return topNode_.getFlatConditionalProbabilities(model,modelChanged, this,0);
		}

		//Branch Length stuff
		public final double getBranchLength() { return branchLength_; }
		public final void setBranchLength(double x) { this.branchLength_ = x; }

		public String toString(UNode caller) {
			if(caller==bottomNode_) {
				return topNode_.toString(this);
			}
			if(caller!=topNode_) {
				throw new RuntimeException("Unknown caller");
			}
			return bottomNode_.toString(this);
		}
		public void testLikelihood(SubstitutionModel model, ConstructionTool tool) {
			testLikelihood(null,model,tool);
		}
		public void testLikelihood(UNode caller, SubstitutionModel model, ConstructionTool tool) {
			System.out.printf("Likelihood: %f, %f, %f\n",
					calculateLogLikelihood(model, true,tool.allocateNewExternalCalculator(), tool),
					calculateLogLikelihood2(model, true,tool.allocateNewExternalCalculator(), tool),
					calculateLogLikelihood3(model, true,tool.allocateNewExternalCalculator(), tool)
					);
			if(caller!=bottomNode_) {
				bottomNode_.testLikelihood(this,model,tool);
			}
			if(caller!=topNode_){
				topNode_.testLikelihood(this,model,tool);
			}
		}

		public ConditionalProbabilityStore getExtendedConditionalProbabilities( SubstitutionModel model, boolean modelChanged, UNode caller, int depth) {
			UNode other = getOther(caller);
			return other.getExtendedConditionalProbabilities(branchLength_,model, modelChanged, this,depth);
		}
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( SubstitutionModel model, boolean modelChanged, UNode caller, LHCalculator.External externalCalculator, ConditionalProbabilityStore extendedStore) {
			UNode other = getOther(caller);
			return other.getExtendedConditionalProbabilities(branchLength_,model, modelChanged, this, externalCalculator, extendedStore);
		}
		public void setup( ConstructionTool tool, Connection[] allConnections) {
			//A call to the right node is not needed as the recursion will cut back that way
			if(tool.hasSequences()) {
				bottomNode_.rebuildPattern( tool );
			}
			for(int i = 0 ; i < allConnections.length ; i++) {
				allConnections[i].centerPatternValid_ = false;
			}
		}

		public void getAllConnections(ArrayList<Connection> store) {
			getAllConnections(store,null);
		}
		public void getAllConnections(ArrayList<Connection> store, UNode caller) {
			store.add(this);
			if(caller!=bottomNode_) {
				bottomNode_.getAllConnections(store,this);
			}
			if(caller!=topNode_) {
				topNode_.getAllConnections(store,this);
			}


		}
		public void getCenterPatternInfo(PatternInfo store, ConstructionTool tool) {
			PatternInfo bottom = bottomNode_.getPatternInfo(this);
			PatternInfo top    = topNode_.getPatternInfo(this);
			tool.build(store, bottom, top);
		}

		public UNode getOther(UNode caller) {
			if(bottomNode_==caller) {
				return topNode_;
			}
			if(topNode_==caller) {
				return bottomNode_;
			}
			throw new RuntimeException("Unknown caller!");
		}

		public final void doNNI(MersenneTwisterFast r) {
			doNNI(r.nextBoolean());
		}
		/**
		 * Does not reconstruct patterns.
		 * If arg is true, swap left child of bottom node with left child of top node,
		 * otherwise swap left child of bottom node with right child of top node.
		 * Returns false if this is a pendant edge (no NNI possible)
		 */
		public boolean doNNI(boolean bottomLeftSwapsWithTopLeft) {
			if (bottomNode_.isLeaf()) { return false; } // No NNI possible
			/*
			 * The bottom node must be the parent to the two other edges. To preserve the direction
			 * of this edge (so 'top' continues to point towards the root), we must swap with a top 
			 * edge which has the same direction - i.e. a child edge of topNode_, not the parent edge.
			 * (The complications are unnecessary but mostly harmless for unrooted trees.)
			 */
			boolean topSwapEdgeIsLeft;
			Connection topSwapEdge;
			if (topNode_.getLeft(this).isTopNode(topNode_)) {
				topSwapEdgeIsLeft = true;
				topSwapEdge = (Connection)topNode_.getLeft(this); // cannot be a Root.
			} else {
				topSwapEdgeIsLeft = false;
				topSwapEdge = (Connection)topNode_.getRight(this); // cannot be a Root.
			}
			Connection bottomSwapEdge;
			if (bottomLeftSwapsWithTopLeft ^ topSwapEdgeIsLeft) {
				bottomSwapEdge = (Connection)bottomNode_.getLeft(this); // cannot be a Root.
			} else {
				bottomSwapEdge = (Connection)bottomNode_.getRight(this); // cannot be a Root.
			} 
			bottomNode_.swapConnection(bottomSwapEdge,topNode_,topSwapEdge);
			return true;
		}

		public double calculateLogLikelihood2(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool) {
			PatternInfo pi = getCenterPatternInfo(tool);
			final ConditionalProbabilityStore bottomConditionalProbabilityProbabilties =
				bottomNode_.getFlatConditionalProbabilities(model, modelChanged, this,0);
			final ConditionalProbabilityStore topConditionalProbabilityProbabilties =
				topNode_.getExtendedConditionalProbabilities(branchLength_, model, modelChanged,  this,0);
			return calculator.calculateLogLikelihood(model, pi, bottomConditionalProbabilityProbabilties,topConditionalProbabilityProbabilties, this.edgeParams_);
		}
		public double calculateLogLikelihood3(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool) {
			PatternInfo pi = getCenterPatternInfo(tool);
			final ConditionalProbabilityStore bottomConditionalProbabilityProbabilties =
				bottomNode_.getExtendedConditionalProbabilities(branchLength_, model, modelChanged,  this,0);				
			final ConditionalProbabilityStore topConditionalProbabilityProbabilties =
				topNode_.getFlatConditionalProbabilities(model, modelChanged, this,0);
			return calculator.calculateLogLikelihood(model, pi, bottomConditionalProbabilityProbabilties,topConditionalProbabilityProbabilties, this.edgeParams_);
		}
		// MDW: changed to make this the function actually used, other one evaluated at other end of edge compared to calculator.
		// (which matters for non-time-reversable models, but this is still just a hack - proper evaluation of NTR models is more complex.)
		public double calculateLogLikelihood(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool) {
		 PatternInfo pi = getCenterPatternInfo(tool);
		 final ConditionalProbabilityStore bottom = bottomNode_.getFlatConditionalProbabilities(model, modelChanged, this,0);
		 final ConditionalProbabilityStore top    =    topNode_.getFlatConditionalProbabilities(model, modelChanged, this,0);
		 return calculator.calculateLogLikelihood(branchLength_, edgeParams_, model, pi, bottom,top,tool.newConditionalProbabilityStore(false));
	 }
		public SiteDetails calculateSiteDetails(SubstitutionModel model, boolean modelChanged, LHCalculator.External calculator, ConstructionTool tool) {
			PatternInfo pi = getCenterPatternInfo(tool);
			final ConditionalProbabilityStore bottom = bottomNode_.getFlatConditionalProbabilities(model, modelChanged, this,0);
			final ConditionalProbabilityStore top    =    topNode_.getFlatConditionalProbabilities(model, modelChanged, this,0);
		  return calculator.calculateSiteDetailsUnrooted(branchLength_, edgeParams_, model,pi,bottom,top,tool.newConditionalProbabilityStore(false));
		}
		// For debugging: return string of branch length plus stringification of ModelEdgeParameters, if any. 
		public String debugString() {
			if (edgeParams_ == null) {
				return Double.toString(branchLength_);
			} else {
				return Double.toString(branchLength_) + edgeParams_.toString(); 
			}
		}
	}

// -==--=-=-=-=-=-=-=-=-=-=-=-=-==--==-=-=--==--==-=--==-=--=-==-=--==--==-=-=-
// ==== OptimisationHandler ====
// -==--=-=-==--=-=-=-==--=-=-=-=-==--=-=-=-==-=-=-=-=-=-===--=-==--=-==--=-==-
	/*
	 * OptimisationHandler exists for the optimiseBranchLength method, which gets and edge (Connection) passed to it.
	 * There are two subclasses: UniOptimisationHandler and MultiOptimisationHandler.
	 * In almost all cases, the only relevant parameter to be optimised on the edge is the edge length.
	 * In this case, we use UniOptimisationHandler. If there are other edge parameters (stored in a ModelEdgeParameter object)
	 * which need to be optimised at the same time, we use MultiOptimisationHandler. We can ask the model which one
	 * is needed by the SubstitutionModel.mustOptimiseMultipleParametersOnEdge() method.
	 */
	
	/*
	 *  MDW: We have a function (likelihood) of many variables (edge lengths, models, patterns...)
	 *  but we want to make just one variable (an edge length) changeable to an exterior optimizer.
	 *  So we make a UnivariateFunction object and hide all of the remaining variables inside this
	 *  object. (I've added edgeParams to the store of hidden information.) 
	 */
	
	private static abstract class OptimisationHandler {
		protected ConditionalProbabilityStore bottomFlatConditionalProbabilities_;
		protected ConditionalProbabilityStore topFlatConditionalProbabilities_;

		protected final SubstitutionModel model_;
		protected PatternInfo currentPatternInfo_;
		protected ModelEdgeParameters edgeParams_ = null;
		protected final LHCalculator.External calculator_;
		protected final ConditionalProbabilityStore tempStore_;
		protected final ConstructionTool tool_;
		public OptimisationHandler(SubstitutionModel model,  ConstructionTool tool) {
			this.calculator_ = tool.allocateNewExternalCalculator();
			this.tempStore_ = tool.newConditionalProbabilityStore(false);
			this.tool_ = tool;
			this.model_ = model;
		}
		public abstract void setup(Connection c, boolean modelChanged);
		/**
		 * Optimise the branch length of a certain connection
		 *
		 * @param c
		 * @param modelChanged
		 * @return maximum
		 */
		public abstract double optimiseBranchLength(Connection c, boolean modelChanged); 

	}
	
	private static final class UniOptimisationHandler extends OptimisationHandler implements UnivariateFunction {

		private ModelEdgeParameters edgeParams_ = null;
		private final UnivariateMinimum um_;
		
		public UniOptimisationHandler(SubstitutionModel model,  ConstructionTool tool) {
			super(model, tool);
			this.um_ = new UnivariateMinimum();
		}
		/**
		 * Optimise the branch length of a certain connection
		 *
		 * @param c
		 * @param modelChanged
		 * @return maximum
		 */
		public double optimiseBranchLength(Connection c, boolean modelChanged) {
			setup(c,modelChanged);
			um_.findMinimum(c.getBranchLength(),this);
			
			c.setBranchLength(um_.minx);
			return -um_.fminx;
	}

		// MDW: Recalculates (if necessary) the conditional probabilities for the nodes at each end of the Connection.
		// Creates a combo PatternInfo from each node's subtree PatternInfo. 
		public void setup(Connection c, boolean modelChanged) {
			this.bottomFlatConditionalProbabilities_ = c.getBottomFlatConditionalProbabilities(model_, modelChanged);

			this.topFlatConditionalProbabilities_ = c.getTopFlatConditionalProbabilities(model_,modelChanged);
			this.currentPatternInfo_ = c.getCenterPatternInfo(tool_);
			this.edgeParams_ = c.edgeParams_;  // MDW: added. 
		}

		// MDW: Uses the info cached in 'setup' to very quickly recalculate the tree likelihood for different edge lengths.
		public double evaluate(double argument) {
			// MDW: added edgeParams_
			return -calculator_.calculateLogLikelihood(argument,edgeParams_,model_,currentPatternInfo_,bottomFlatConditionalProbabilities_,topFlatConditionalProbabilities_, tempStore_);
		}

		public double getLowerBound() {	return MINIMUM_BRANCH_LENGTH; }
		public double getUpperBound() { return MAXIMUM_BRANCH_LENGTH; }

	}
	
	
	private static final class MultiOptimisationHandler extends OptimisationHandler implements MultivariateFunction {

		private final MultivariateMinimum mm_;
		private int numArguments_=-1;
		
		public MultiOptimisationHandler(SubstitutionModel model,  ConstructionTool tool) {
			super(model, tool);
			this.mm_ = new OrthogonalSearch(); // TODO: Is this the best MultivariateMinimum to use?
		}
		/**
		 * Optimise the branch length of a certain connection
		 *
		 * @param c
		 * @param modelChanged
		 * @return maximum
		 */
		public double optimiseBranchLength(Connection c, boolean modelChanged) {
			setup(c,modelChanged);
			double[] xvec = new double[numArguments_];
			xvec[0] = c.getBranchLength();
			for (int i=1; i<numArguments_; i++) {
				xvec[i] = edgeParams_.getParameter(i-1);
			}
			double minimum = mm_.findMinimum(this,xvec);
			
			c.setBranchLength(xvec[0]);
			for (int i=1; i<numArguments_; i++) {
				edgeParams_.setParameter(xvec[i],i-1); // this.edgeParams_ is a pointer to c.edgeParams_.
			}

			return -minimum;
	    }

		// MDW: Recalculates (if necessary) the conditional probabilities for the nodes at each end of the Connection.
		// Creates a combo PatternInfo from each node's subtree PatternInfo. 
		public void setup(Connection c, boolean modelChanged) {
			this.bottomFlatConditionalProbabilities_ = c.getBottomFlatConditionalProbabilities(model_, modelChanged);

			this.topFlatConditionalProbabilities_ = c.getTopFlatConditionalProbabilities(model_,modelChanged);
			this.currentPatternInfo_ = c.getCenterPatternInfo(tool_);
			this.edgeParams_ = c.edgeParams_;  
			this.numArguments_ = edgeParams_.getNumParameters() + 1;
		}

		/*
		 * MultivariateFunction methods:
		 */
		
		// MDW: Uses the info cached in 'setup' to very quickly recalculate the tree likelihood for different edge lengths.
		public double evaluate(double argument[]) {
			double branchLength = argument[0];
			for (int i=1; i<numArguments_; i++) {
				edgeParams_.setParameter(argument[i], i-1);
			}
			// MDW: added edgeParams_. Split line for ease of debugging.
			double fx = -calculator_.calculateLogLikelihood(branchLength,edgeParams_,model_,currentPatternInfo_,bottomFlatConditionalProbabilities_,topFlatConditionalProbabilities_, tempStore_);
			return fx;
		}

		public int getNumArguments() { return numArguments_; }
		public double getLowerBound(int n) 
		{
			if (n==0) {
				return MINIMUM_BRANCH_LENGTH;
			} else {
				return edgeParams_.getLowerLimit(n-1);
			}
		}
		public double getUpperBound(int n) { 
			if (n==0) {
				return MAXIMUM_BRANCH_LENGTH;
			} else {
				return edgeParams_.getUpperLimit(n-1);
			}
		}
		public OrthogonalHints getOrthogonalHints() {
			// TODO: Currently the only existing ModelEdgeParameter object doesn't have
			// any orthogonal hints, but really should try to make use of any such hints if they exist.
			return null;
		}

	}
	
	
	
// -==--=-=-=-=-=-=-=-=-=-=-=-=-==--==-=-=--==--==-=--==-=--=-==-=--==--==-=-=-
// ==== OptimisationHandler ====
// -==--=-=-==--=-=-=-==--=-=-=-=-==--=-=-=-==-=-=-=-=-=-===--=-==--=-==--=-==-
	// MDW: Any use of ModelEdgeParameters with this is untested.
	private static final class NNIOptimisationHandler implements UnivariateFunction {
		private final double[][][] transitionProbabiltityStore_ ;
		private ConditionalProbabilityStore bottomFlatConditionalProbabilities_;
		private ConditionalProbabilityStore topFlatConditionalProbabilities_;

		private final ConditionalProbabilityStore bottomFlatStore_;
		private final ConditionalProbabilityStore topFlatStore_;
		private final ConditionalProbabilityStore bottomLeftExtendedStore_;
		private final ConditionalProbabilityStore bottomRightExtendedStore_;
		private final ConditionalProbabilityStore topLeftExtendedStore_;
		private final ConditionalProbabilityStore topRightExtendedStore_;

		private final ConditionalProbabilityStore tempStore_;

		private final SubstitutionModel model_;
		private ModelEdgeParameters edgeParams_; 
		private final int numberOfCategories_;
		private final int numberOfStates_;
		private PatternInfo currentPatternInfo_;

		private final LHCalculator.External calculator_;
		private final UnivariateMinimum um_;
		private final ConstructionTool tool_;
		private final PatternInfo bottomPatternStore_;
		private final PatternInfo topPatternStore_;
		private final PatternInfo centerPatternStore_;
		private final Connection[] allConnections_;
		public NNIOptimisationHandler(Connection[] allConnections, SubstitutionModel model, ConstructionTool tool) {
			numberOfStates_ = model.getDataType().getNumStates();
			this.allConnections_ = allConnections;
			this.calculator_ = tool.allocateNewExternalCalculator();

			numberOfCategories_ = model.getNumberOfTransitionCategories();
			this.transitionProbabiltityStore_ = new double[numberOfCategories_][numberOfStates_][numberOfStates_];
			this.model_ = model;
			this.tool_ = tool;
			this.bottomPatternStore_ = tool.constructFreshPatternInfo(true);
			this.topPatternStore_ = tool.constructFreshPatternInfo(true);
			this.centerPatternStore_ = tool.constructFreshPatternInfo(true);
			this.tempStore_ = tool.newConditionalProbabilityStore(false);
			this.bottomFlatStore_ = tool.newConditionalProbabilityStore(false);
			this.topFlatStore_ = tool.newConditionalProbabilityStore(false);
			this.bottomLeftExtendedStore_ = tool.newConditionalProbabilityStore(false);
			this.bottomRightExtendedStore_ = tool.newConditionalProbabilityStore(false);
			this.topLeftExtendedStore_ = tool.newConditionalProbabilityStore(false);
			this.topRightExtendedStore_ = tool.newConditionalProbabilityStore(false);

			this.um_ = new UnivariateMinimum();
		}
		private final ConditionalProbabilityStore calculateFlat(PatternInfo patternInfo, ConditionalProbabilityStore resultStore, ConditionalProbabilityStore leftExtend, ConditionalProbabilityStore rightExtend) {
			calculator_.calculateFlat(patternInfo,leftExtend,rightExtend,resultStore);
			return resultStore;
		}

		/**
		 * Optimise the branch length of a certain connection
		 *
		 * @param c
		 * @param modelChanged
		 * @return True likelihood (not negated_
		 */
		public double optimiseSimulataneousNNIBranchLength(Connection c, boolean modelChanged) {
			edgeParams_ = c.edgeParams_; 
			 UNode bottomNode = c.getBottom();
			 UNode topNode = c.getTop();
			PatternInfo bottomLeftPI   = bottomNode.getLeftPatternInfo(c);
			PatternInfo topRightPI =    topNode.getRightPatternInfo(c);
			if(bottomLeftPI==null||topRightPI==null) {
				//One of the nodes is a leaf or something incompatible for NNI
				setup(c,modelChanged);
				um_.findMinimum(c.getBranchLength(),this);
				c.setBranchLength(um_.minx);
				return -um_.fminx;
			}
// =--=-=-=-==-=-
			PatternInfo bottomRightPI = bottomNode.getRightPatternInfo(c);
			PatternInfo topLeftPI = topNode.getLeftPatternInfo(c);

			ConditionalProbabilityStore bottomLeftExtended  = bottomNode.getLeftExtendedConditionalProbabilities( model_,false,c,calculator_, bottomLeftExtendedStore_);
			ConditionalProbabilityStore bottomRightExtended = bottomNode.getRightExtendedConditionalProbabilities(model_,false,c,calculator_, bottomRightExtendedStore_);
			ConditionalProbabilityStore topLeftExtended     =    topNode.getLeftExtendedConditionalProbabilities( model_,false,c,calculator_, topLeftExtendedStore_);
			ConditionalProbabilityStore topRightExtended    =    topNode.getRightExtendedConditionalProbabilities(model_,false,c,calculator_, topRightExtendedStore_);

			tool_.build(topPatternStore_, topLeftPI,topRightPI);
			tool_.build(bottomPatternStore_, bottomLeftPI,bottomRightPI);
			tool_.build(centerPatternStore_, bottomPatternStore_,topPatternStore_);

			ConditionalProbabilityStore bottomFlat = calculateFlat(c.getBottomPatternInfo(),bottomFlatStore_,bottomLeftExtended,bottomRightExtended);
			ConditionalProbabilityStore topFlat = calculateFlat(c.getTopPatternInfo(),topFlatStore_, topLeftExtended,topRightExtended);
			setPatternInfo(bottomFlat,topFlat,c.getCenterPatternInfo(tool_));

			um_.findMinimum(c.getBranchLength(),this);
			final double standardMaxFX = -um_.fminx;
			final double standardMaxX = um_.minx;


			//Swap bottom left with top right
			tool_.build(bottomPatternStore_, topLeftPI,bottomRightPI);
			tool_.build(topPatternStore_, bottomLeftPI,topRightPI);
			tool_.build(centerPatternStore_, bottomPatternStore_,topPatternStore_);

			bottomFlat = calculateFlat(bottomPatternStore_,bottomFlatStore_,topLeftExtended,bottomRightExtended);
			topFlat = calculateFlat(topPatternStore_,topFlatStore_, bottomLeftExtended,topRightExtended);
			setPatternInfo(bottomFlat,topFlat,centerPatternStore_);

			um_.findMinimum(c.getBranchLength(),this);
			final double bottomLeftSwapMaxX = um_.findMinimum(c.getBranchLength(),this);
			final double bottomLeftSwapMaxFX = -um_.fminx;
// =--=-=-=-==-=-

			//Diag swap
			tool_.build(bottomPatternStore_, topRightPI,bottomRightPI);
			tool_.build(topPatternStore_, topLeftPI,bottomLeftPI);
			tool_.build(centerPatternStore_, bottomPatternStore_,topPatternStore_);

			bottomFlat = calculateFlat(bottomPatternStore_,bottomFlatStore_,topRightExtended,bottomRightExtended);
			topFlat = calculateFlat(topPatternStore_,topFlatStore_, topLeftExtended,bottomLeftExtended);



			setPatternInfo(bottomFlat,topFlat,centerPatternStore_);
			um_.findMinimum(c.getBranchLength(),this);

			final double diagSwapMaxX = um_.minx;
			final double diagSwapMaxFX = -um_.fminx;



			if(standardMaxFX>diagSwapMaxFX) {
				if(standardMaxFX>bottomLeftSwapMaxFX) {
					//Standard wins!
					c.setBranchLength(standardMaxX);
					return standardMaxFX;
				} else {
					//bottom left wins
					c.doNNI(true);
					c.setup(tool_,allConnections_);
					c.setBranchLength(bottomLeftSwapMaxX);
					return bottomLeftSwapMaxFX;
				}
			} else {
				if(diagSwapMaxFX>bottomLeftSwapMaxFX) {
					//Diag diag wins!
					c.setBranchLength(diagSwapMaxX);
					c.doNNI(false);
					c.setup(tool_,allConnections_);
					return diagSwapMaxFX;
				} else {
					//bottom left wins
					c.setBranchLength(bottomLeftSwapMaxX);
					c.doNNI(true);
					c.setup(tool_,allConnections_);
					return bottomLeftSwapMaxFX;
				}
			}
		}
		private void setup(Connection c, boolean modelChanged) {
			this.edgeParams_ = c.edgeParams_;
			this.bottomFlatConditionalProbabilities_ = c.getBottomFlatConditionalProbabilities(model_, modelChanged);
			this.topFlatConditionalProbabilities_ = c.getTopFlatConditionalProbabilities(model_,modelChanged);
			this.currentPatternInfo_ = c.getCenterPatternInfo(tool_);
		}
		public void setPatternInfo(
				ConditionalProbabilityStore bottomFlatConditionalProbabilities,
				ConditionalProbabilityStore topFlatConditionalProbabilities,
				PatternInfo pi
			) {
			this.bottomFlatConditionalProbabilities_ = bottomFlatConditionalProbabilities;
			this.topFlatConditionalProbabilities_ = topFlatConditionalProbabilities;
			this.currentPatternInfo_ = pi;
		}

		public double evaluate(double argument) {
			return -calculator_.calculateLogLikelihood(argument,edgeParams_,model_,currentPatternInfo_,bottomFlatConditionalProbabilities_,topFlatConditionalProbabilities_,tempStore_);
		}

		public double getLowerBound() {	return 0; }
		public double getUpperBound() { return 1; }

	}


	private static final class ConstructionTool {
		private final String[] names_;
		private final int[][] sequences_;
		private final int numberOfStates_;
		private final int numberOfCategories_;
		private final int numberOfSites_;
		private final DataType dt_;
		private int nextConnectionIndex_ = 0;
		private final ArrayList allUNodes_ = new ArrayList();
		private final LHCalculator.Generator calcGenerator_;
		public ConstructionTool(Alignment alignment, int numberOfStates, int numberOfCategories, LHCalculator.Factory calculatorFactory) {
		  if(alignment!=null) {
				DataType dt = alignment.getDataType();
				if( dt.isAmbiguous() ) {
					this.dt_ = dt.getAmbiguousVersion().getSpecificDataType();
				} else {
					this.dt_ = alignment.getDataType();
				}

				this.numberOfSites_ = alignment.getSiteCount();
			} else {
			  this.dt_ = null;
				this.numberOfSites_ = 0;
			}
			this.numberOfStates_ = numberOfStates;
			this.numberOfCategories_ = numberOfCategories;
			if(alignment!=null) {
		  	this.names_ = Identifier.getNames(alignment);
				this.sequences_ = pal.alignment.AlignmentUtils.getAlignedStates( alignment, dt_.getNumStates() );
				this.calcGenerator_ = calculatorFactory.createSeries(numberOfCategories,dt_);
			} else {
				this.names_ = null;
			  this.sequences_ = null;
				this.calcGenerator_ = null;
			}

		}
		public boolean hasSequences() { return sequences_!=null&&sequences_.length>0; }
		public PatternInfo constructFreshPatternInfo(boolean binaryPattern) {
			return new PatternInfo(numberOfSites_,binaryPattern);
		}
		public final ConditionalProbabilityStore newConditionalProbabilityStore(boolean isForLeaf) {
			return calcGenerator_.createAppropriateConditionalProbabilityStore(isForLeaf);
		}
		public final int allocateNextConnectionIndex() {
			return nextConnectionIndex_++;
		}
		public LHCalculator.Internal allocateNewInternalCalculator() {
			if(calcGenerator_!=null) {
				return calcGenerator_.createNewInternal();
			}
			return null;
		}
		public LHCalculator.External allocateNewExternalCalculator() {
			if(calcGenerator_!=null) {
				return calcGenerator_.createNewExternal();
			}
			return null;
		}
		public LHCalculator.Leaf createNewLeafCalculator(int[] patternStateMatchup, int numberOfPatterns) {
			if(calcGenerator_!=null) {
				return calcGenerator_.createNewLeaf( patternStateMatchup, numberOfPatterns );
			}
			return null;
		}
		public int build(PatternInfo beingBuilt, PatternInfo p1, PatternInfo p2) {
			return beingBuilt.build(p1,p2,numberOfSites_);
		}
		public final int allocateNextUNodeIndex(UNode node) {
			int index = allUNodes_.size();
			allUNodes_.add(node);
			return index;
		}
		public final UNode[] getOrderedNodes() {
			UNode[] result = new UNode[allUNodes_.size()];
			allUNodes_.toArray(result);
			return result;
		}
		public DataType getDataType() { return dt_; }
		public final int getNumberOfSites() { return numberOfSites_; }

		public int[] getSequence(String name) {
			if(sequences_==null) {
			  return null;
			}
			for(int i = 0 ; i < names_.length ; i++) {
				if(name.equals(names_[i])) {
					return sequences_[i];
				}
			}
			throw new IllegalArgumentException("Unknown sequence:"+name);
		}
		public int getNumberOfStates() { return numberOfStates_; }
		public int getNumberOfTransitionCategories() { return numberOfCategories_; }
	}
	
	
	/**
	 * Method 'staticTest' runs a simple test case, printing results to provided PrintWriter.
	 * 
	 * Note: see also 'possiblyRootedMLSearcher.testLikelihood' method, 
	 * which evaluates likelihood on each Connection/Root within the tree.
	 * 
	 * The test case is simple enough to be hand-calculatable. We have one site, four
	 * leaves and a two state model (actually implemented as a four state model where 
	 * two of the states are inaccessible). The same process acts on all branches. If 
	 * the root base frequencies are equal to the equilibrium frequencies of the process, 
	 * we have a time reversible model.
	 * 
	 * We test in three ways: 
	 * (1) as a time reversible model
	 * (2) as a non-time-reversible model where the root frequencies are the equilibrium frequencies
	 *  (3) as a non-time-reversible model where the root frequencies are not the equilibrium frequencies.
	 * Cases (1) and (2) must produce the same result. For any case, we must get the same result
	 * no matter where we evaluate the likelihood. (In principle, any edge can have the likelihood evaluated
	 * at either end.)
	 *  
	 * Here is the model:
	 * The The available states are A and G (probability vectors are in the order [A,G].)
	 * The Markov matrix (in row-sums-to-one orientation) on most edges is
	 * M=[[.8,.2],[.3,.7]]. A couple of edges are double length, so their Markov matrix is 
	 * M^2 = [[.7,.3],[.45,.55]]. The equilibrium frequencies of this process are [.6,.4]. 
	 * The (rooted) tree is ((u:1,v:1):2,(w:1,x:2):1).
	 * The observed pattern is (u,v,w,x)=(A,G,G,G) 
	 *  
	 * In the time reversible case (cases 1 and 2) the likelihood is 0.0346875 (exactly.) 
	 * In case 3, we use root frequencies [0.2,0.8] and then the likelihood is 0.0475 (exactly.) 
	 * 
	 * TODO: Extend the test to use two rate classes.
	 */
	/*
	 * Here is some Octave code which calculates all of the partial conditional probabilities
	 * all over the tree, using the calculations required for an NTR model. The interior nodes are 
	 * labeled 's' (parent of u, v), 't' (parent of w,x) and 'r' (root.) 
	 * Every vector gets a 3 letter name. The third letter is 't' or 'b' indicating whether it 
	 * is at the top or the bottom of an edge. The first two letters specify the edge: 'rs' is the edge
	 * between nodes r and s. The order of these two letters indicates the direction in which we are looking.
	 * So 'usb' is the at the bottom of edge 'us', and in direction 'u' to 's': i.e. likelihood of everything below 'u'
	 * (As u is a leaf in state A, usb = [1,0].) By contrast, 'sub' is the conditional probabilities of everything
	 * above node 'u' (including the branch from u to s). It is sub=[0.0345875,0.046187] in cases (1) and (2),
	 * [0.0475,0.0705] in case (3).
	 * 
	 *  The Octave code:
	 *  M=[[.8,.3];[.2,.7]]
	 *	M2=M*M
	 *	usb=[1,0]
	 *	vsb=[0,1]
	 *	wtb=[0,1]
	 *	xtb=[0,1]
	 *	ust=usb*M
	 *	vst=vsb*M
	 *	wtt=wtb*M
	 *	xtt=xtb*M2
	 *	srb=ust.*vst
	 *	trb=wtt.*xtt
	 *	srt=srb*M2
	 *	trt=trb*M
	 *	
	 *	rrb=[0.2,0.8]
	 *	rst=trt.*rrb
	 *	rtt=srt.*rrb
	 *	rsb=rst*M2'
	 *	rtb=rtt*M'
	 *	sut=rsb.*vst
	 *	svt=rsb.*ust
	 *	twt=rtb.*xtt
	 *	txt=rtb.*wtt
	 *	sub=sut*M'
	 *	svb=svt*M'
	 *	twb=twt*M'
	 *	txb=txt*M2'
	 *
	 *  The code above calculates for case 3. For case 2, use rbb=[0.6,0.4]. (For case 1, final answer is the
	 *  same as case 2, but calculations are slightly different.)
	 * 
	 */
	public static void staticTest(PrintWriter out) {
		
		// The expected answers:
		double target1 = 0.0346875;
		double target2 = target1;
		double target3 = 0.0475;
		
		// The two state model:
		NeoRateMatrix myRateMatrix = new GeneralREVRateMatrix(4,new int[] {0,1,0,0,1,0}, new double[] {0});
		double[] equilibriumFreqs = new double[] {0.6, 0.0, 0.4, 0.0}; 
		SubstitutionModel trModel = new SingleClassSubstitutionModel(myRateMatrix, Nucleotides.DEFAULT_INSTANCE, equilibriumFreqs, true);
		// and the non-equilibrium extension of this model:
		NonTimeReversibleSubstitutionModel ntrModel = new GeneralNonEquilibriumModel(trModel);
		
		// The test tree:
		Tree testTree = TreeUtils.stringToTree("((u:1,v:1):2,(w:1,x:2):1);");
		testTree = TreeUtils.getScaled(testTree,0.48*Math.log(2));
		// This scaling is required to give us Markov matrix [[.8,.2][.7,.3]] on edges of length 1.
		
		// The test alignment (single site):
		String alignmentString = new String(">u\nA\n>v\nG\n>w\nG\n>x\nG\n");
		StringReader sr = new StringReader(alignmentString);
		BufferedReader br = new BufferedReader(sr);
		Alignment alignment = null;
		try {
			alignment = AlignmentReaders.readFastaSequences(br,Nucleotides.DEFAULT_INSTANCE);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// For test 3
		double[] testRootFreqs = new double[] {0.2, 0.0, 0.8, 0.0};
		
		// The various calculator factories to test:
		LHCalculator.Factory[] calcFactories = new LHCalculator.Factory[] {
				SimpleLHCalculator.getFactory(),
				FastFourStateLHCalculator.getFactory(),
				SimpleModelFastFourStateLHCalculator.getFactory()
		};
		//public PossiblyRootedMLSearcher(Node root, Alignment alignment, SubstitutionModel model, LHCalculator.Factory calcFactory, long seed) {
		
		for (LHCalculator.Factory factory : calcFactories) {
			out.printf("Testing calculator factory %s:\n", factory.getClass().getCanonicalName() );
				
			// Test case 1:
			PossiblyRootedMLSearcher prmls = new PossiblyRootedMLSearcher(testTree,alignment,trModel);
			double answer1 = Math.exp(prmls.calculateLogLikelihood());		
			out.printf("PossiblyRootedMLSearcher test 1: likelihood = %10.8g, expected %10.8g\n", answer1, target1);
			for (int i=0; i<prmls.allConnections_.length; i++) {
				double ll = prmls.allConnections_[i].calculateLogLikelihood(trModel, false, prmls.tool_.allocateNewExternalCalculator(), prmls.tool_);
				if (Math.abs(Math.exp(ll)-target1)>1e-12) {
					out.printf("Test fails on edge %d\n", i);
				}
			}
			
			// Test case 2:
			ntrModel.setRootFrequencies(equilibriumFreqs);
			prmls = new PossiblyRootedMLSearcher(testTree,alignment,ntrModel);
			double answer2 = Math.exp(prmls.calculateLogLikelihood());		
			out.printf("PossiblyRootedMLSearcher test 2: likelihood = %10.8g, expected %10.8g\n", answer2, target2);
			for (int i=0; i<prmls.allConnections_.length; i++) {
				double ll = prmls.allConnections_[i].calculateLogLikelihood(ntrModel, false, prmls.tool_.allocateNewExternalCalculator(), prmls.tool_);
				if (Math.abs(Math.exp(ll)-target2)>1e-12) {
					out.printf("Test fails on edge %d\n", i);
				}
			}
			
			// Test case 3:
			ntrModel.setRootFrequencies(testRootFreqs);
			prmls = new PossiblyRootedMLSearcher(testTree,alignment,ntrModel);
			double answer3 = Math.exp(prmls.calculateLogLikelihood());		
			out.printf("PossiblyRootedMLSearcher test 3: likelihood = %10.8g, expected %10.8g\n", answer3, target3);
			for (int i=0; i<prmls.allConnections_.length; i++) {
				double ll = prmls.allConnections_[i].calculateLogLikelihood(ntrModel, false, prmls.tool_.allocateNewExternalCalculator(), prmls.tool_);
				if (Math.abs(Math.exp(ll)-target3)>1e-12) {
					out.printf("Test fails on edge %d\n", i);
				}
			}
		} // loop over factories
	}
}