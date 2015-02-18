package pal.tree;

import java.io.PrintWriter;
import java.util.Iterator;

import pal.util.AlgorithmCallback;
/**
 * Given a list of taxa names, produce all possible binary unrooted tree topologies on those taxa. 
 * 
 * @author Michael Woodhams
 */
public class AllUnrootedTreeIterator implements Iterator<SimpleTree>, TreeIterator {

	private static final double BRANCH_LENGTH = 1; // TODO: (possibly) allow user to specify branch length to use.
	private AllRootedTreeIterator rootedIterator_;
	private SimpleNode finalNode_;
	private boolean forceBinary_;
	
	/**
	 * An iterator over binary unrooted trees (or more accurately, arbitrarily rooted.)
	 * You can get output trees in two formats: root node is ternary (all the rest are binary,
	 * i.e. we've taken a binary unrooted tree and called one of the internal nodes the root)
	 * or root node is binary (i.e. we've taken a binary unrooted tree and inserted a root into
	 * an arbitrary edge.) 
	 * The tree topology search methods require the root node to be binary.
	 * @param taxonNames
	 * @param forceBinary: if true, add an extra root node so that all nodes are binary.
	 */
	public AllUnrootedTreeIterator(String taxonNames[], boolean forceBinary) {
		// make a node of the final taxon
		finalNode_ = new SimpleNode(taxonNames[taxonNames.length-1],BRANCH_LENGTH);
		// and a rooted iterator over all but the final taxon
		rootedIterator_ = AllRootedTreeIterator.newIterator(taxonNames, taxonNames.length-1);
		forceBinary_ = forceBinary;
	}
	
	/**
	 * An iterator over binary unrooted trees (or more accurately, arbitrarily rooted.)
	 * Will return a tree where root node is ternary, all others binary.
	 * @param taxonNames
	 */
	public AllUnrootedTreeIterator(String taxonNames[]) {
		this(taxonNames,false);
	}

	public boolean isMoreTrees() {
		return rootedIterator_.isMoreTrees();
	}

	public SimpleTree getNextTree(AlgorithmCallback callback) {
		// Note restOfTree is a clone of iterator's output - otherwise 
		// subsequent trees from iterator will share nodes, and any
		// topology manipulations on one will make the other invalid.
		SimpleNode restOfTree = new SimpleNode(rootedIterator_.getNextNode());
		SimpleNode root;
		if (forceBinary_) {
			root = new SimpleNode();
			root.addChild(restOfTree);
			root.addChild(new SimpleNode(finalNode_));
		} else {
			root = restOfTree;
			root.addChild(finalNode_);
		}
		return new SimpleTree(root);
	}
	
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
	
	public static void test(PrintWriter out) {
		out.println("Test should produce 15 distinct unrooted trees");
		AllUnrootedTreeIterator testIterator = new AllUnrootedTreeIterator(new String[] {"A","B","C","D","E"},true);
		while (testIterator.hasNext()) {
			Tree tree = testIterator.next();
			TreeUtils.printNH(tree, out, false, false);
			if (!TreeUtils.isTreeValid(tree)) { out.println("Invalid tree - bad links!"); }
			out.flush();
		}
		testIterator = new AllUnrootedTreeIterator(new String[] {"A","B","C","D","E","G"},false);
		int count = 0;
		while (testIterator.hasNext()) { 
			Tree tree = testIterator.next(); 
			if (!TreeUtils.isTreeValid(tree)) { out.println("Invalid tree - bad links!"); }
			count++; 
		}
		out.printf("For 6 taxa, generated %d tree topologies (105 expected)\n",count);
		testIterator = new AllUnrootedTreeIterator(new String[] {"A","B","C","D","E","G","H"},false);
		count = 0;
		while (testIterator.hasNext()) { 
			Tree tree = testIterator.next(); 
			if (!TreeUtils.isTreeValid(tree)) { out.println("Invalid tree - bad links!"); }
			count++; 
		}
		out.printf("For 7 taxa, generated %d tree topologies (945 expected)\n",count);
	}
}
