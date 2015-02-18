package pal.substmodel;

import pal.algorithmics.Markable;


/**
 * Unofficial modification to PAL by Michael Woodhams (School of IT, University of Sydney)
 * 
 * This is a do-nothing base class. Any unusual substitution models which require extra parameters on an
 * node should subclass this one. In TreeSearch.UnrootedMLSearcher, the UNode will have one of these
 * objects, and will call updates on it from Connection.setup(), during tree topology changes.
 * 
 * Used by the MosaicSubstitutionModel to track which models apply within which subtrees.
 * 
 * To attach one of these to an AttributeNode, store it as an attribute using 
 * ModelNodeParameters.attributeLabel as the label.
 */

public interface ModelNodeParameters extends Markable {
	public static final String ATTRIBUTE_LABEL = "ModelNodeParameters"; 
	
	// Should just have a constructor which takes and 'int', but Java won't let me insist on that
	// whether this is an interface or an abstract class.
	// Calls to this method are allowed to delete all data.
	public void setNumberOfAdjacentEdges(int adjacentEdges);
	// Gets called during tree traversal after topology change in UnrooteedMLSearcher.InternalNode.rebuildPattern
	// For the i-th adjacent node, the current node is that node's connection number nodeBackpointers[i]. 
	// Similarly for edges: 0 for left, 1 for right stored in edgeBackPointers.
	// 'indexTowardsRoot' is the connection number pointing towards the root. '-1' if this node is the root. 
	public void update(ModelNodeParameters[] adjacentMNPs, 
			int[] nodeBackPointers, 
			ModelEdgeParameters[] adjacentMEPs, 
			int[] edgeBackPointers, 
			int indexTowardsRoot,
			boolean firstPass);
	// Return a new instance of whatever class this object is.
	public ModelNodeParameters getInstance(int adjacentEdges);
	// reset to a known 'blank' state. Gets called on each node by UnrooteedMLSearcher.InternalNode.rebuildPattern
	// when tree topology has changed
	public void resetForTopologyChanged();
}


