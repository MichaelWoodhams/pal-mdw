// LHCalculator.java
//
// (c) 1999-2003 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package pal.eval;

/**
 * <p>Title: LHCalculator </p>
 * <p>Description: An LHCalculator object must be treated as a stateful, single threaded object that can be used
 * for calculating components in an overall likelihood calculation. </p>
 * <p>History<br>
 *  <ul>
 *    <li>25/10/2003 Added leaf handling interface </li>
 *    <li>30/3/2004 Changed certain methods to more intelligent ones (relating to posterior distribution of sites). Added abstract External class.
 *    <li>25/02/2009 MDW: change so every call which takes a distance also takes a ModelEdgeParameters argument.
 *  </ul>
 * </p>
 * @author Matthew Goode
 * @version 1.0
 * @note needs to have the use of the word likelihood altered in certain cases (to conditional probability)
 *
 */

 /*
  * Note on notation:
  * Wherever we have a 'left' and 'right' something as parameters, we are doing a calculation at a node, and the left/right are two incoming
  * edges (Connections). Wherever we have a 'top' and 'bottom' something as parameters, we are doing a calculation at an edge/Connection, and 
  * these are the incoming whatevers from each end of the edge. 
  * This is how things *should* be, but it used to be that 'left' and 'right' were used for both situations, and it may be 
  * that the process of refactoring all those variable names was imperfect. Violations should be corrected. 
  */


import pal.datatype.*;
import pal.substmodel.*;

public interface LHCalculator {
  /**
   * The External calculator does not maintain any state and is approapriate for
   * calculation where a store is provided
   */
  public static interface External extends java.io.Serializable {
    /**
     * This method is applied at a node to take component-wise product of conditional probabilities from two incoming edges
     * and present them as the likelihood at the base of the other edge.
     * @param centerPattern the pattern information
     * @param leftConditionalProbabilities Implementations must not overwrite or change
     * @param rightConditionalProbabilities Implementations must not overwrite or change
     * @param resultStore Where to stick the created categoryPatternState information
     * @note calls to getLastConditionalProbabilities() does not have to be valid after call this method
     */
    public void calculateFlat( PatternInfo centerPattern, ConditionalProbabilityStore leftConditionalProbabilities, ConditionalProbabilityStore rightConditionalProbabilities, ConditionalProbabilityStore resultStore );

    /**
     * Combines cond. prob.s from two incoming edges and then propagates them a distance along an outgoing edge.
     * @param distance
     * @param forwardsInTime: true if the node we are at is at the beginning (in time) of the outgoing edge.
     *            For calculation in time reversible models, we treat everything as being towards the root,
     *            so this arg should be 'false'. 
     * @param edgeParams: any extra parameters on the edge (usually null.)
     * @param model
     * @param centerPattern the pattern information
     * @param resultStore Where to stick the created categoryPatternState information
     * @param edgeParams any edge-dependent variables other than its length (distance)
     * @param leftConditionalProbabilities Implementations must not overwrite or change
     * @param rightConditionalProbabilities Implementations must not overwrite or change
     * @note calls to getLastConditionalProbabilities() does not have to be valid after call this method
     */
    public void calculateExtended( double distance,
    							   boolean forwardsInTime,
    							   ModelEdgeParameters edgeParams,
                                   SubstitutionModel model,
                                   PatternInfo centerPattern,
                                   ConditionalProbabilityStore leftConditionalProbabilities,
                                   ConditionalProbabilityStore rightConditionalProbabilities,
                                   ConditionalProbabilityStore resultStore );

    /**
     * Extend the conditionals back in time by some distance, with some model
     * @param distance The distance to extend by
     * @param model The model to use
     * @param conditionalProbabilities The probabilities to extend
     */
    public void calculateSingleExtendedDirect(
									double distance,
									ModelEdgeParameters edgeParams,
									SubstitutionModel model,
									int numberOfPatterns,
                                    ConditionalProbabilityStore conditionalProbabilities
                                  );
	/**
     * Extend the conditionals back in time by some distance, with some model
	 * @param distance The distance to extend by
     * @param model The model to use
     * @param baseConditionalProbabilities The probabilities to extend
     * @param resultConditionalProbabilities The probabilities to extend
     */
    public void calculateSingleExtendedIndirect(
									double distance,
									ModelEdgeParameters edgeParams,
									SubstitutionModel model,
									int numberOfPatterns,
                                    ConditionalProbabilityStore baseConditionalProbabilities,
                                    ConditionalProbabilityStore resultConditionalProbabilities
                                  );


    /**
     * Calculate the likelihood given two sub trees (left, right) and their flat (unextended) likelihood probabilities
     * (top/bottom ordering doesn't matter to time reversible models.) 
     * @param distance
     * @param edgeParams
     * @param model
     * @param centerPattern the pattern information
     * @param bottomFlatConditionalProbabilities (i.e. closer to the leaves)
     * @param topFlatConditionalProbabilities (i.e. closer to the root, because trees grow upside down.) 
     * @param tempStore may be used internally to calculate likelihood
	 * @return the log likelihood
     */
    public double calculateLogLikelihood( double distance, 
    								   ModelEdgeParameters edgeParams,
    								   SubstitutionModel model,
                                       PatternInfo centerPattern,
                                       ConditionalProbabilityStore bottomFlatConditionalProbabilities,
                                       ConditionalProbabilityStore topFlatConditionalProbabilities,
                                       ConditionalProbabilityStore tempStore
                                       );
    
    /**
     * Calculate the likelihood given two sub trees (bottom, top) and their extended likelihood probabilities
     * Basically does component-wise product of bottom probability vector, top probability vector and equilibrium frequencies.
     * In this case, the top/bottom labels are arbitrary: calculation is symmetric in these two arguments.
     * Used for calculating likelihood on an edge, where one of top/bottom is a flat CP, the other is an extended one.
     * @param model
     * @param centerPattern the pattern information
     * @param leftConditionalProbabilities
     * @param rightConditionalProbabilities
     * @param edgeParams normally null, see the class for more details.
     * @return the Log likelihood
     */
    public double calculateLogLikelihood(  SubstitutionModel model, 
    								   PatternInfo centerPattern,
                                       ConditionalProbabilityStore bottomConditionalProbabilities,
                                       ConditionalProbabilityStore topConditionalProbabilities,
                                       ModelEdgeParameters edgeParams);

    /**
     * Calculate the likelihood given the conditional probabilities at the root
     * (dot product of a single conditional probability with the root/equilibrium base frequencies.)
     * @param model The substitution model used
     * @param patternWeights the weights of each pattern
     * @param numberOfPatterns the number of patterns
     * @param conditionalProbabilities The conditionals
     * @return the Log likelihood
     */
    public double calculateLogLikelihoodSingle( SubstitutionModel model, int[] patternWeights, int numberOfPatterns,
                                       ConditionalProbabilityStore conditionalProbabilityStore);


    /**
     * Calculate the conditional probabilities of each pattern for each category
     * @param model
     * @param centerPattern the pattern information
     * @param leftConditionalProbabilities
     * @param rightConditionalProbabilities
     * @param categoryPatternLogLikelihoodStore after call will hold a matrix of values in the form [cat][pattern], where [cat][pattern] represents the site probability under a particular category/class, *not* multiplied by the category probability or pattern weights
     */
    public SiteDetails calculateSiteDetailsRooted( SubstitutionModel model,
      PatternInfo centerPattern,
      ConditionalProbabilityStore leftConditionalProbabilitiesStore,
      ConditionalProbabilityStore rightConditionalProbabilitiesStore
      );

    /**
     * Calculate the conditional probabilities of each pattern for each category
     * @param distance The distance between the two nodes
     * @param edgeParams any edge-dependent variables other than its length (distance)
     * @param model
     * @param centerPattern the pattern information
     * @param bottomConditionalProbabilities
     * @param topConditionalProbabilities
     * @param categoryPatternLogLikelihoodStore after call will hold a matrix of values in the form [cat][pattern], 
     * where [cat][pattern] represents the site probability under a particular category/class, *not* multiplied 
     * by the category probability or pattern weights
     */
	public SiteDetails calculateSiteDetailsUnrooted( 
				double distance, 
				ModelEdgeParameters edgeParams,
				SubstitutionModel model,
				PatternInfo centerPattern,
				ConditionalProbabilityStore bottomConditionalProbabilitiesStore,
				ConditionalProbabilityStore topConditionalProbabilitiesStore,
				ConditionalProbabilityStore tempStore
    );
  } //End of class External
// =================================================================================================
// ================= Internal ======================================================================
// =================================================================================================
  /**
   * The Internal calculator may maintain state and is an appropriate permanent attachment
   * to internal nodes of the tree structure
   */
  public static interface Internal {
    /**
     * calculate flat probability information (not extended over a branch).
     * @param centerPattern the pattern information
     * @param leftConditionalProbabilities Implementations should be allowed to overwrite in certain cases
     * @param rightConditionalProbabilities Implementations should be allowed to overwrite in certain cases
     * @return true if results built from cached information
     * @note An assumption may be made that after a call to this method the leftConditionals and rightConditionals are not used again!
     */
    public ConditionalProbabilityStore calculateFlat( PatternInfo centerPattern, ConditionalProbabilityStore leftConditionalProbabilities, ConditionalProbabilityStore rightConditionalProbabilities );
    
    /**
    * Combines cond. prob.s from two incoming edges and then propagates them a distance along an outgoing edge.
    * @param distance
     * @param edgeParameters Optional model parameters specific to this edge
     * @param forwardsInTime: true if the node we are at is at the beginning (in time) of the outgoing edge.
     *            For calculation in time reversible models, we treat everything as being towards the root,
     *            so this arg should be 'false'. 
     * @param model
     * @param centerPattern the pattern information
     * @param bottomConditionalProbabilities
     * @param topConditionalProbabilities
     * @param modelChangedSinceLastCall this should be true if the substituion model has altered since the last call to this method on this particular object, false otherwise
    * @return resulting conditional probabilities
    * @note An assumption may be made that after a call to this method the leftConditionals and rightConditionals are not used again!
    */
    public ConditionalProbabilityStore calculateExtended( 
    		double distance, 
    		boolean forwardsInTime,
    		ModelEdgeParameters edgeParameters, 
    		SubstitutionModel model, 
    		PatternInfo centerPattern,
    	    final ConditionalProbabilityStore bottomConditionalProbabilities, 
    	    final ConditionalProbabilityStore topConditionalProbabilities, 
    	    boolean modelChangedSinceLastCall );

    public ConditionalProbabilityStore calculatePostExtendedFlat( 
    		double distance,
    		ModelEdgeParameters edgeParameters, 
    		SubstitutionModel model, 
    		PatternInfo centerPattern, 
    		final ConditionalProbabilityStore bottomConditionalProbabilities,
    		final ConditionalProbabilityStore topConditionalProbabilities, 
    		boolean modelChangedSinceLastCall );

  } //End of Internal

// =================================================================================================
// ================= Leaf ==========================================================================
// =================================================================================================
	/**
	 * A LHCalculator.Leaf object is attached to each leaf node and can be used to calculated conditional probabilities across the related branch.
	 * Allows for quick implementations as well as implementations that cope correctly with ambiguous characters
	 * @note Should not be made serializable!
	 */
	public static interface Leaf {
		public ConditionalProbabilityStore getFlatConditionalProbabilities();
		public ConditionalProbabilityStore getExtendedConditionalProbabilities( 
				double distance, 
				ModelEdgeParameters edgeParam, 
				SubstitutionModel model, 
				boolean modelChanged);
		/**
		 * Create a new Leaf calculator that has exactly the same properties as this one (but is different such that it may be used independently)
		 * @return a copy of this leaf calculator
		 */
		public Leaf getCopy();
	}

// =================================================================================================
// ================= Factory ==========================================================================
// =================================================================================================	
	public static interface Factory extends java.io.Serializable  {
		public Generator createSeries( int numberOfCategories, DataType dt );
	}
	
// =================================================================================================
// ================= Generator ==========================================================================
// =================================================================================================
	public static interface Generator extends java.io.Serializable  {
		/**
		 * Create anew leaf calculator
		 * @param patternStateMatchup The sequence as reduced to patterns. This should just be one state per pattern.
		 * For example given a sequence [ 0, 1,0,1,3,0] a patternMatchup may be [0,1,3] (the first element is the first
		 * pattern, which is state 0, the second element is the second pattern which is 1, and the third element is the
		 * third pattern (novel pattern) which is state 3)
		 * @param numberOfPatterns The number of patterns in the patternStateMatchup array
		 * @return a leaf calculator object
		 */
		public Leaf createNewLeaf(int[] patternStateMatchup, int numberOfPatterns);

		public Leaf createNewLeaf(int[] patternStateMatchup, int numberOfPatterns, Generator parentGenerator );

		public External createNewExternal();

		public Internal createNewInternal();

		public boolean isAllowCaching();

	    /**
	     * An obscure method, primarily used by the High Accuracy calculator
	     * @param parentGenerator A reference to an encompassing generator (that may for example
	     * wish to impose it's own choice on the creation of ConditionalProbabilityStores)
	     * @throws IllegalArgumentException Generator does not allow being a subserviant generator
	     * @return
	     */
	    public External createNewExternal( Generator parentGenerator ) throws IllegalArgumentException;
	
	    /**
	     * An obscure method, primarily used by the High Accuracy calculator
	     * @param parentGenerator A reference to an encompasing generator (that may for example
	     * wish to impose it's own choice on the creation of ConditionalProbabilityStores)
	     * @throws IllegalArgumentException Generator does not allow being a subserviant generator
	     * @return
	     */
	    public Internal createNewInternal( Generator patentGenerator ) throws IllegalArgumentException;
	
	    public ConditionalProbabilityStore createAppropriateConditionalProbabilityStore( boolean isForLeaf );

	}
// ======================================================================================
	public abstract class AbstractExternal {
		public final SiteDetails calculateSiteDetailsUnrooted(
	 				double distance,
	 				ModelEdgeParameters edgeParams,
	 				SubstitutionModel model,
	 				PatternInfo centerPattern,
	 				ConditionalProbabilityStore bottomFlatConditionalProbabilities,
	 				ConditionalProbabilityStore topFlatConditionalProbabilities,
	 				ConditionalProbabilityStore tempStore
		) {
			double[][] store = new double[model.getNumberOfTransitionCategories()][centerPattern.getNumberOfPatterns()];
	 		calculateCategoryPatternProbabilities(distance, edgeParams, model, centerPattern, bottomFlatConditionalProbabilities,topFlatConditionalProbabilities,tempStore, store);
	 		double[] siteLikelihoods = calculateSiteLikelihoods(
	 					store,
	 					model.getTransitionCategoryProbabilities(),
	 					model.getNumberOfTransitionCategories(),
	 					centerPattern.getSitePatternMatchup(),
	 					centerPattern.getNumberOfSites()
	 		);
	 		return SiteDetails.Utils.create(
	 					store,false,model,
	 					centerPattern.getNumberOfPatterns(),
	 					centerPattern.getSitePatternMatchup(),
	 					centerPattern.getNumberOfSites(),
	 					siteLikelihoods
	 		);
		}
		private final double[] calculateSiteLikelihoods( double[][] conditionals, final double[] catProbabilities, int numberOfCategories, int[] sitePatternMatchup,  int numberOfSites) {
		    final double[] siteLikeihoods = new double[numberOfSites];
			for(int site= 0 ; site < numberOfSites ; site++) {
				double total = 0;
				int pattern = sitePatternMatchup[site];
				for(int cat = 0 ; cat < numberOfCategories ; cat++) {
				  total+=catProbabilities[cat]*conditionals[cat][pattern];
				}
				siteLikeihoods[site] = total;
			}
			return siteLikeihoods;
		}
		public final SiteDetails calculateSiteDetailsRooted(SubstitutionModel model,
				PatternInfo centerPattern,
				ConditionalProbabilityStore leftConditionalProbabilitiesStore,
				ConditionalProbabilityStore rightConditionalProbabilitiesStore
		) {
			double[][] store = new double[model.getNumberOfTransitionCategories()][centerPattern.getNumberOfPatterns()];
			calculateCategoryPatternProbabilities(model,centerPattern,leftConditionalProbabilitiesStore,rightConditionalProbabilitiesStore,store);
			final  double[] siteLikelihoods = calculateSiteLikelihoods(
				    store,
					model.getTransitionCategoryProbabilities(),
					model.getNumberOfTransitionCategories(),
					centerPattern.getSitePatternMatchup(),centerPattern.getNumberOfSites()
			);
			return SiteDetails.Utils.create(store,false,model,centerPattern.getNumberOfPatterns(),centerPattern.getSitePatternMatchup(),centerPattern.getNumberOfSites(),siteLikelihoods);
		}
		protected abstract void calculateCategoryPatternProbabilities(
				double distance,
				ModelEdgeParameters edgeParams,
				SubstitutionModel model, 
				PatternInfo centerPattern,
				ConditionalProbabilityStore bottomFlatConditionalProbabilities,
				ConditionalProbabilityStore topFlatConditionalProbabilities,
				ConditionalProbabilityStore tempStore,
				double[][] categoryPatternLogLikelihoodStore
		);

		protected abstract void calculateCategoryPatternProbabilities(
				SubstitutionModel model, PatternInfo centerPattern,
				ConditionalProbabilityStore leftConditionalProbabilities,
				ConditionalProbabilityStore rightConditionalProbabilities,
				double[][] categoryPatternLikelihoodStore
		);
	}
}