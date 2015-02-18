package pal.tree;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.NoSuchElementException;
import pal.util.AlgorithmCallback;

/**
 * Given a list of taxa names, produce all possible binary rooted tree topologies on those taxa. 
 * Create an object via the static newIterator method - there is no constructor.
 * 
 * Interfaces Iterator<SimpleTree> and TreeIterator are pretty much synonymous, but have different
 * method names.
 * 
 * @author Michael Woodhams
 */

public abstract class AllRootedTreeIterator implements Iterator<SimpleTree>, TreeIterator {
	
	private static final double BRANCH_LENGTH = 1; // TODO: (possibly) allow user to specify branch length to use.
	
	public static AllRootedTreeIterator newIterator(String[] taxonNames) {
		return newIterator(taxonNames, taxonNames.length);
	}
	public static AllRootedTreeIterator newIterator(String[] taxonNames, int nTaxa) {
		if (nTaxa == 1) {
			return new TrivialIterator(taxonNames[0]);
		} else if (nTaxa == 2) {
			return new TrivialIterator(taxonNames[0], taxonNames[1]);
		} else if (nTaxa <=0) {
			throw new IllegalArgumentException(); 
		} else {
			return new NonTriviaIterator(taxonNames, nTaxa);
		}
	}
	
	public SimpleTree getNextTree(AlgorithmCallback callback) {
		// Note we clone the iterator's output - otherwise 
		// subsequent trees from iterator will share nodes, and any
		// topology manipulations on one will make the other invalid.
		return new SimpleTree(new SimpleNode(getNextNode()));
	}
	
	public abstract void resetIterator();
	/** returns next tree, but just as the root node (none of the overhead of a full Tree object.)
	 * WARNING: you must clone this node if you want to manipulate tree topology OR if you want
	 * to hold more than one output of the iterator in memory at once.
	 * Best to use 'getNextTree' instead which is safe.
	 * (I should probably make this method 'private' but I'd need to figure out how to let
	 * AllUnrootedTreeIterator have access to it.)
	 * @return
	 */
	public abstract SimpleNode getNextNode();
	
	// alternate method names to implement java.util.Iterator
	/**
	 * Redirects to isMoreTrees()
	 */
	public boolean hasNext() { return isMoreTrees(); }
	/**
	 * redirects to getNextTree()
	 */
	public SimpleTree next() { return getNextTree(null); }
	/**
	 * Unimplemented
	 */
	public void remove() { throw new UnsupportedOperationException(); }
	

	
	private static class TrivialIterator extends AllRootedTreeIterator {
		// For one or two taxa: there is only one tree.
		private SimpleNode theTree_;
		private boolean hasNext_;
		
		public TrivialIterator(String name) {
			theTree_ = new SimpleNode(name,BRANCH_LENGTH);
			hasNext_ = true;
		}
		
		public TrivialIterator(String name1, String name2) {
			theTree_ = new SimpleNode(new SimpleNode[] {new SimpleNode(name1,BRANCH_LENGTH),new SimpleNode(name2,BRANCH_LENGTH)},BRANCH_LENGTH);
			hasNext_ = true;
		}
		
		public boolean isMoreTrees() { return hasNext_; }
		public SimpleNode getNextNode() {
			if (hasNext_) {
				hasNext_ = false;
				return theTree_;
			} else {
				throw new NoSuchElementException();
			}
		}
		public void resetIterator() { hasNext_ = true; }
	} // class TrivialIterator

	// NonTrivialIterator: for where there is more than one tree to iterate through (>= 3 taxa.)
	private static class NonTriviaIterator extends AllRootedTreeIterator {
		private String[] taxa_;
		private int nTaxa_;
		private int split_; //The most recent split tried, encoded as bits. An int has plenty of bits for this purpose.
		private AllRootedTreeIterator leftIterator_;
		private AllRootedTreeIterator rightIterator_;
		// We need to return all possible left tree/right tree pairs. This keeps the current left tree while
		// we cycle through all the possible right trees.
		private Node leftTree_; 
		private String[] leftSubset_;
		private String[] rightSubset_; // non-persistent scratch storage space.
		private int leftSubsetCount_, rightSubsetCount_;
		
		/**
		 * An iterator which produces all possible binary tree topologies labeled by the first
		 * nTaxa identifiers of taxaIds
		 * @param taxaIDs
		 * @param nTaxa
		 */
		public NonTriviaIterator(String[] taxaNames, int nTaxa) {
			nTaxa_ = nTaxa;
			taxa_ = new String[nTaxa];
			System.arraycopy(taxaNames, 0, taxa_, 0, nTaxa);
			leftSubset_ = new String[nTaxa-1];
			rightSubset_ = new String[nTaxa-1];
			resetIterator();
		}

		public void resetIterator() {
			split_ = (1<<(nTaxa_-1))-1; // e.g. nTaxa_ = 4, split_ <- 111(binary).
			newSubsets();
		}
		
		public boolean isMoreTrees() {
			return rightIterator_.isMoreTrees() || leftIterator_.isMoreTrees() || split_ > 0;
		}

		public SimpleNode getNextNode() {
			if (rightIterator_.isMoreTrees()) {
				SimpleNode root = new SimpleNode();
				root.addChild(leftTree_);
				root.addChild(rightIterator_.getNextNode());
				root.setBranchLength(BRANCH_LENGTH);
				return root;
			} else if (leftIterator_.isMoreTrees()) {
				leftTree_ = leftIterator_.getNextNode();
				rightIterator_.resetIterator();
				return getNextNode();
			} else {
				// left and right iterators are exhausted, so we need a new split
				if (split_ == 0) {
					throw new NoSuchElementException();
				}
				newSubsets();

				return getNextNode();
			}
		} // getNextnode()
		
		/*
		 * Generate next subsets from 'split_'.
		 * Sets left/rightIterator_, left/rightSubset_, left/rightSubsetCount_.
		 */
		private void newSubsets() {
			int bits = split_--;
			leftSubsetCount_ = 0;
			rightSubsetCount_ = 0;
			for (int i=0; i<nTaxa_-1;i++) {
				if ((bits & 1) == 0) {
					leftSubset_[leftSubsetCount_++] = taxa_[i];
				} else {
					rightSubset_[rightSubsetCount_++] = taxa_[i];
				}
				bits >>= 1;
			}
			// final taxon is always put in leftSubset_
			leftSubset_[leftSubsetCount_++] = taxa_[nTaxa_-1];
			leftIterator_  = newIterator(leftSubset_,  leftSubsetCount_);
			rightIterator_ = newIterator(rightSubset_, rightSubsetCount_);
			leftTree_ = leftIterator_.getNextNode();
		} // newSubsets()
	} // class NonTrivialIterator
	
	public static void test(PrintWriter out) {
		out.println("Test should produce 15 distinct rooted trees");
		AllRootedTreeIterator testIterator = newIterator(new String[] {"A","B","C","D"});
		while (testIterator.hasNext()) {
			TreeUtils.printNH(testIterator.next(), out, false, false);
			out.flush();
		}
	}
} // class AllRootedTreeIterator
