package pal.substmodel;

import pal.math.SetOfSmallIntegers;


public class MosaicModelNodeParameters implements ModelNodeParameters {
	// One set for each adjacent edge: specifies which submodels are represented in the subtree
	// at the far end of that edge.
	private SetOfSmallIntegers[] submodelsInSubtree_;
	private int[] subtreeEdgeParameterCounts_;
	private MosaicModelEdgeParameters[] adjacentEdgeParameters_;
	
	private SetOfSmallIntegers[] markSubmodelsInSubtree_;
	private int[] markSubtreeEdgeParameterCounts_;
	private MosaicModelEdgeParameters[] markAdjacentEdgeParameters_;
	
	// Constructor for leaf nodes - specify the set of submodels which apply to this taxon.
	public MosaicModelNodeParameters(SetOfSmallIntegers submodels) {
		submodelsInSubtree_ = new SetOfSmallIntegers[] {submodels};
		subtreeEdgeParameterCounts_ = new int[] {0}; // we have no subtree, so count is zero.
		adjacentEdgeParameters_ = null;
		
		markSubmodelsInSubtree_ =null;
		markSubtreeEdgeParameterCounts_ = null;
		markAdjacentEdgeParameters_ = null;
	}
	
	public MosaicModelNodeParameters(int adjacentEdges) {
		this.setNumberOfAdjacentEdges(adjacentEdges);
		markSubmodelsInSubtree_ =null;
		markSubtreeEdgeParameterCounts_ = null;
		markAdjacentEdgeParameters_ = null;
	}
	
	public void setNumberOfAdjacentEdges(int adjacentEdges) {
		submodelsInSubtree_ = new SetOfSmallIntegers[adjacentEdges];
		subtreeEdgeParameterCounts_ = new int[adjacentEdges];
		adjacentEdgeParameters_ = new MosaicModelEdgeParameters[adjacentEdges];
		for (int i=0; i<adjacentEdges; i++) {
			submodelsInSubtree_[i] = new SetOfSmallIntegers();
			subtreeEdgeParameterCounts_[i] = -1;
		}
	}
	
	public SetOfSmallIntegers getSubmodelsInSubtree(int subtree) {
		return submodelsInSubtree_[subtree];
	}
	
	public MosaicModelNodeParameters getInstance(int adjacentEdges) {
		return new MosaicModelNodeParameters(adjacentEdges);
	}
	
	public void resetForTopologyChanged() {
		for (SetOfSmallIntegers set : submodelsInSubtree_) {
			set.setToEmpty();
		}
	}

	// For the subtree rooted by this node as seen from Connection number 'direction',
	// find the set of submodels represented in this subtree. Associate that subset
	// with the 'direction' Connection.
	// (If this is a leaf node, nothing will happen: 'i!=direction' will never be true.)
	public void setUnionOfSubmodels(
			ModelNodeParameters[] adjacentMNPs, 
			int[] nodeBackPointers, 
			ModelEdgeParameters[] adjacentMEPs, 
			int[] edgeBackPointers, 
			int direction) 
	{
		SetOfSmallIntegers setToBuild = submodelsInSubtree_[direction];
		for (int i=0; i<submodelsInSubtree_.length; i++) {
			// nodeBackPointers[i] == -1 if the adjacent ConnectionOrRoot is the Root
			if (i!=direction && nodeBackPointers[i]>=0) {
				MosaicModelNodeParameters adjacentNode = (MosaicModelNodeParameters)adjacentMNPs[i];
				setToBuild.unionEquals(adjacentNode.submodelsInSubtree_[nodeBackPointers[i]]);
			}
		}
	}
	
	/*
	 * Get count of free edge weight parameters on the subtree in the given direction.
	 * From this nodes point of view, the subtree on (and including) the edge in the given direction.
	 */
	private int getSingleSubtreeEdgeParametersCount(int direction) {
		if (subtreeEdgeParameterCounts_[direction] == -1) {
			// We have not got a cached value, so need to ask neighbouring node to calculate it for us
			if (adjacentEdgeParameters_[direction] != null) {
				subtreeEdgeParameterCounts_[direction] = adjacentEdgeParameters_[direction].getSubtreeEdgeParametersCount(this);
			} else {
				subtreeEdgeParameterCounts_[direction] = 0;
			}
		}
		return subtreeEdgeParameterCounts_[direction];
	}
	/*
	 * Return number of free edge weight parameters on the subtree rooted at this node
	 * as seen from the calling edge.
	 */
	public int getSubtreeEdgeParametersCount(MosaicModelEdgeParameters caller) {
		
		int callerIndex = -1;
		for (int i=0; i<adjacentEdgeParameters_.length; i++) {
			if (caller == adjacentEdgeParameters_[i]) {
				callerIndex = i;
				break;
			}
		}
		if (callerIndex == -1) { throw new RuntimeException("Bug in Mosaic model edge parameter counting - couldn't find caller");}
		// could combine the two loops in this routine, but at risk of an infinite recursion if pointers are wrong.
		int sum = 0;
		for (int i=0; i<adjacentEdgeParameters_.length; i++) {
			if (i!=callerIndex) {
				sum += getSingleSubtreeEdgeParametersCount(i);
			}
		}
		return sum;
	}
	
	/**
	 * Returns orphan generation of lowest orphan gen of adjacent edge except for 'caller'
	 * and the set of submodels applying to the lowest orphan generation edge(s)
	 * @param caller
	 * @param subset used to return set of submodels.
	 * @return
	 */
	public int getOrphanGeneration(MosaicModelEdgeParameters caller, SetOfSmallIntegers subset) {
		SetOfSmallIntegers edgeSubset = new SetOfSmallIntegers();
		int bestOrphanGeneration = Integer.MAX_VALUE-1; // -1 to be smaller than value from Root edge, saves a little work
		for (MosaicModelEdgeParameters edge : adjacentEdgeParameters_) {
			if (edge != caller) {
				int edgeOrphanGen = edge.getOrphanGeneration(this, edgeSubset);
				if (edgeOrphanGen < bestOrphanGeneration) {
					// Found new lowest neighbouring orphan generation
					bestOrphanGeneration = edgeOrphanGen;
					subset.set(edgeSubset);
				} else if (edgeOrphanGen == bestOrphanGeneration) {
					// Found an equal to current best orphan generation
					subset.unionEquals(edgeSubset);
				}
			}
		}
		return bestOrphanGeneration; // and value of 'subset' returned by reference.
	}
	
	public void update(
			ModelNodeParameters[] adjacentMNPs, 
			int[] nodeBackPointers, 
			ModelEdgeParameters[] adjacentMEPs, 
			int[] edgeBackPointers, 
			int indexTowardsRoot,
			boolean firstPass) 
	{
		int n = submodelsInSubtree_.length;
		if ( adjacentMNPs.length != n || adjacentMEPs.length != n || indexTowardsRoot >= n || indexTowardsRoot <-1) {
			// This code will never be executed
			throw new IllegalArgumentException();
		}
		if (firstPass) {
			adjacentEdgeParameters_ = new MosaicModelEdgeParameters[adjacentMEPs.length];
			for (int i=0; i<adjacentMEPs.length; i++) {
				adjacentEdgeParameters_[i] = (MosaicModelEdgeParameters)adjacentMEPs[i];
				subtreeEdgeParameterCounts_[i] = -1;
			}
			if (indexTowardsRoot == -1) {
				// Special case: we're at the root. All subtrees have their submodels sets pointing towards us up to date.
				// We need to update this nodes submodel sets in all directions.
				for (int i=0; i<n; i++) {
					this.setUnionOfSubmodels(adjacentMNPs, nodeBackPointers, adjacentMEPs, edgeBackPointers, i);
				}
			} else {
				// Normal case: only in the towards-root direction 
				this.setUnionOfSubmodels(adjacentMNPs, nodeBackPointers, adjacentMEPs, edgeBackPointers, indexTowardsRoot);
			}
		} else {
			// Second pass: the towards-root one has been set, but now we have the information to be able to se the rest.
			 // We need to update this nodes submodel sets in all directions.
			if (indexTowardsRoot != -1) {
				// If we are at root (indexTowardsRoot == -1) then first pass completely updated us - do nothing.
				for (int i=0; i<n; i++) {
					if (i != indexTowardsRoot) {
						this.setUnionOfSubmodels(adjacentMNPs, nodeBackPointers, adjacentMEPs, edgeBackPointers, i);
					}
				}
			}
		}
	}
	
	public void mark() {
		markSubmodelsInSubtree_ = submodelsInSubtree_.clone();
		markSubtreeEdgeParameterCounts_ = subtreeEdgeParameterCounts_.clone();
		markAdjacentEdgeParameters_ = adjacentEdgeParameters_.clone();
	}
	public void undoToMark() {
		submodelsInSubtree_ = markSubmodelsInSubtree_;
		subtreeEdgeParameterCounts_ = markSubtreeEdgeParameterCounts_;
		adjacentEdgeParameters_ = markAdjacentEdgeParameters_;
	}
	
	public String toString() {
		StringBuffer out = new StringBuffer("MosaicModelNodeParameters(");
		for (SetOfSmallIntegers set : submodelsInSubtree_) {
			out.append(set.toString());
		}
		out.append(")");
		return out.toString();
	}

}
