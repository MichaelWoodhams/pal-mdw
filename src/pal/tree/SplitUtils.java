// SplitUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package pal.tree;

import pal.misc.*;

/**
 * utilities for split systems
 *
 * @version $Id: SplitUtils.java,v 1.6 2002/06/05 23:23:14 matt Exp $
 *
 * @author Korbinian Strimmer
 * 
 * TODO: Can be made more efficient by reimplementing using SetOfSmallIntegers
 * (which has not yet been shifted into PAL at the time of writing.) 
 */
public class SplitUtils
{
	//
	// Public stuff
	//

	/**
	 * creates a split system from a tree
	 * (using a pre-specified order of sequences)
	 *
	 * @param idGroup  sequence order for the matrix
	 * @param tree
	 */
	
	/*
	 * Old buggy code: if root of Tree was beside a leaf, the leaf split would be included.
	 * Otherwise, the split between two sides of the root would get included twice.
	 * My solution is clunky and not too efficient (wasted creation of boolean[]'s 
	 * to start with.) There must be a better solution where we just chose one
	 * of the internal nodes to not send to 'getSplit' - but also worry about
	 * non-binary trees.

	public static SplitSystem oldgetSplits(IdGroup idGroup, Tree tree)
	{
		int size = tree.getInternalNodeCount()-1;
		SplitSystem splitSystem = new SplitSystem(idGroup, size);

		boolean[][] splits = splitSystem.getSplitVector();

		for (int i = 0; i < size; i++)
		{
			getSplit(idGroup, tree.getInternalNode(i), splits[i]);
		}

		return splitSystem;
	}
	 */
	public static SplitSystem getSplits(IdGroup idGroup, Tree tree)
	{	
		int numInternalNodes = tree.getInternalNodeCount()-1;
		int numLeaves = idGroup.getIdCount();
		boolean[][] rawSplits = new boolean[numInternalNodes][];

		int numGoodSplits = 0;
		for (int i = 0; i < numInternalNodes; i++)
		{
			boolean[] split = new boolean[numLeaves];
			getSplit(idGroup, tree.getInternalNode(i), split);
			
			// Add split to list if:
			// * it isn't a leaf
			// * it isn't already there
			int count = 0;
			for (int j=0; j<split.length; j++) {
				count += (split[j] ? 1 : 0);
			}
			if (count == 1 || count == numLeaves-1) {
				break;
			}
			boolean foundMatch = false; // found a matching split?
			for (int j=0; j<numGoodSplits && !foundMatch; j++) {
				boolean foundMismatch = false; // found a point of disagreement within this particular split?
				for (int k=0; k<numLeaves && !foundMismatch; k++) {
					foundMismatch = (rawSplits[j][k] != split[k]);
				}
				foundMatch = !foundMismatch;
			}
			if (!foundMatch) {
				rawSplits[numGoodSplits++] = split;
			}
		}
		
		
		SplitSystem splitSystem = new SplitSystem(idGroup, numGoodSplits);
		boolean[][] splits = splitSystem.getSplitVector();

		for (int i = 0; i < numGoodSplits; i++)
		{
			splits[i] = rawSplits[i];
		}


		return splitSystem;
	}



	/**
	 * creates a split system from a tree
	 * (using tree-induced order of sequences)
	 *
	 * @param tree
	 */
	public static SplitSystem getSplits(Tree tree)
	{
		IdGroup idGroup = TreeUtils.getLeafIdGroup(tree);

		return getSplits(idGroup, tree);
	}



	/**
	 * get split for branch associated with internal node
	 *
	 * @param idGroup order of labels
	 * @param internalNode Node
	 * @param boolean[] split
	 * 
	 * This method has a known bug - see source code for details.
	 */
	/*
	 * TODO: Fix this bug.
	 * Following gives distance of 0.5, not 0. First tree gives a split system with
	 * one split duplicated (inefficient but perhaps not incorrect), the second
	 * tree gives a 'leaf split' (one single 'false') which should not be there.
	System.out.printf("Bug: R-F distance between identical trees = %f\n",
			TreeUtils.getRobinsonFouldsDistance(
					TreeUtils.stringToTree("(((A,B),(C,D)),(((E,F),(G,H)),Z));"),
					TreeUtils.stringToTree("((((A,B),(C,D)),((E,F),(G,H))),Z);")
		    )
	);
	*/

	public static void getSplit(IdGroup idGroup, Node internalNode, boolean[] split)
	{
		if (internalNode.isLeaf() || internalNode.isRoot())
		{
			throw new IllegalArgumentException("Only internal nodes (and no root) nodes allowed");
		}

		// make sure split is reset
		for (int i = 0; i < split.length; i++)
		{
			split[i] = false;
		}

		// mark all leafs downstream of the node

		for (int i = 0; i < internalNode.getChildCount(); i++)
		{
			markNode(idGroup, internalNode, split);
		}

		// standardize split (i.e. first index is alway true)
		if (split[0] == false)
		{
			for (int i = 0; i < split.length; i++)
			{
				if (split[i] == false)
					split[i] = true;
				else
					split[i] = false;
			}
		}
	}

	/**
	 * checks whether two splits are identical
	 * (assuming they are of the same length
	 * and use the same leaf order)
	 *
	 * @param s1 split 1
	 * @param s2 split 2
	 */
	public static boolean isSame(boolean[] s1, boolean[] s2)
	{
		boolean reverse;
		if (s1[0] == s2[0]) reverse = false;
		else reverse = true;

		if (s1.length != s2.length)
			throw new IllegalArgumentException("Splits must be of the same length!");

		for (int i = 0; i < s1.length; i++)
		{
			if (reverse)
			{
				// splits not identical
				if (s1[i] == s2[i]) return false;
			}
			else
			{
				// splits not identical
				if (s1[i] != s2[i]) return false;
			}
		}

		return true;
	}

	//
	// Package stuff
	//

	static void markNode(IdGroup idGroup, Node node, boolean[] split)
	{
		if (node.isLeaf())
		{
			String name = node.getIdentifier().getName();
			int index = idGroup.whichIdNumber(name);

			if (index < 0)
			{
				throw new IllegalArgumentException("INCOMPATIBLE IDENTIFIER (" + name + ")");
			}

			split[index] = true;
		}
		else
		{
			for (int i = 0; i < node.getChildCount(); i++)
			{
				markNode(idGroup, node.getChild(i), split);
			}
		}
	}

}
