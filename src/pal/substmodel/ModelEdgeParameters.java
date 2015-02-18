/**
 * Unofficial modification to PAL by Michael Woodhams (School of IT, University of Sydney)
 * 
 * This is a do-nothing base class. Any unusual substitution models which require more parameters on an
 * edge that simply an edge length should store those parameters in a subclass of this one.
 * 
 * For example: we are mixing models on each edge, the weights given to each model being allowed to
 * vary arbitrarily on each edge. 
 * 
 * Example 2: We are allowing a different rate matrix on every edge.
 * 
 * The toString method is used by NodeUtils.printNH(...), so implement something suitable
 */

package pal.substmodel;
import pal.algorithmics.Markable;
import pal.misc.*;

public interface ModelEdgeParameters extends Parameterized, Markable {
	public static final String ATTRIBUTE_LABEL = "ModelEdgeParameters";
	
	// Called from PossiblyRootedMLSearcher during adjustment of tree following a topology change.
	// This edge (Connection) has a 'near' node which has ModelNodeParameters MNP, and it
	// from the near node, this is Connection number 'nearBackPointer'. Ditto for 'far'.
	public void update(
			ModelNodeParameters nearMNP, 
			int nearBackPointer,
			ModelNodeParameters farMNP,
			int farBackPointer,
			boolean firstPass);
	// Set parameters to some default value: used when we want to optimize from a 'fresh start'
	// rather than starting with the current parameters.
	public void initializeParameters();
}
