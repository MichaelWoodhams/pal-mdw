package pal.tree;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import pal.math.SetOfSmallishIntegers;
import pal.misc.IdGroup;
import pal.misc.SimpleIdGroup;

/**
 * Given a collection of trees on the same taxa, produce counts of how often each split occurs.
 * 
 * This is a very simple class, written to the task immediately to hand, and could do with 
 * much elaboration.
 * 
 * @author woodhams
 */

public class BootstrapAnalysis {
	public final static String BOOTSTRAP_ATTRIBUTE = "Bootstrap";
	private IdGroup identifiers_;
	private Map<SetOfSmallishIntegers,Integer> map_; // keys = splits (as SetOfSmallishIntegers), values = count.

	/**
	 * Constructor with a single tree, which both supplies the taxon set and the first set of splits to count.
	 * @param tree
	 */
	public BootstrapAnalysis(Tree tree) {
		identifiers_ = new SimpleIdGroup(tree);
		// TODO: specify initial capacity, load factor for greater efficiency.
		map_ = new HashMap<SetOfSmallishIntegers, Integer>();
		addTree(tree);
	}
	public BootstrapAnalysis(String treeString) {
		this(TreeUtils.stringToTree(treeString));
	}
	// TODO: Add some more constructors, e.g. array of trees
	
	public void addTree(Tree tree) {
		SplitSystem splits = SplitUtils.getSplits(identifiers_, tree);
		int n = splits.getSplitCount();
		for (int i=0; i<n; i++) {
			SetOfSmallishIntegers split = new SetOfSmallishIntegers(splits.getSplit(i));
			if (map_.containsKey(split)) {
				// What a complicated way to add one to something. I hope Java optimizes this down to something sensible.
				map_.put(split, new Integer(1+((Integer)map_.get(split)).intValue()));
			} else {
				map_.put(split, new Integer(1));
			}
		}
	}
	public void addTree(String treeString) {
		addTree(TreeUtils.stringToTree(treeString));
	}

	// TODO: some form of output other than printing.
	public void printSplitCounts(PrintWriter out) {
		Iterator<SetOfSmallishIntegers> splitIterator = map_.keySet().iterator();
		while (splitIterator.hasNext()) {
			SetOfSmallishIntegers split = splitIterator.next();
			out.printf("%3d: %s\n", map_.get(split),setToString(split));
		}
	}
	
	private static final String EMPTY = "";
	private static final String COMMA = ",";
	private String setToString(SetOfSmallishIntegers split) {
		StringBuffer buffer = new StringBuffer("(");
		String separator = EMPTY;
		for (int i=0; i<identifiers_.getIdCount(); i++) {
			if (split.isMember(i)) {
				buffer.append(separator).append(identifiers_.getIdentifier(i).toString());
				separator = COMMA;
			}
		}
		buffer.append(") | (");
		separator = EMPTY;
		for (int i=0; i<identifiers_.getIdCount(); i++) {
			if (!split.isMember(i)) {
				buffer.append(separator).append(identifiers_.getIdentifier(i).toString());
				separator = COMMA;
			}
		}
		buffer.append(")");
		return buffer.toString();
	}
	
	/**
	 * For each internal edge in tree, give it a Bootstrap attribute
	 * containing the % bootstrap support for that edge.
	 * Tree must contain AttributeNodes, not plain nodes.
	 * (At time of writing, PAL has no Nodes which are not AttributeNodes
	 * so this is unlikely to cause problems.)
	 * @param tree
	 * @param nBootStraps - the number of bootstraps performed.
	 */
	public void applyToTree(Tree tree, int nBootStraps) {
		int nNodes = tree.getInternalNodeCount();
		int nLeaves = tree.getExternalNodeCount();
		boolean[] booleanSplit = new boolean[nLeaves];
		for (int i=0; i<nNodes; i++) {
			AttributeNode node = (AttributeNode)tree.getInternalNode(i);
			if (node != tree.getRoot()) {
				SplitUtils.getSplit(identifiers_, node, booleanSplit);
				SetOfSmallishIntegers split = new SetOfSmallishIntegers(booleanSplit);
				int count=0;
				if (map_.containsKey(split)) {
					count = ((Integer)map_.get(split)).intValue();
				}
				// I'm not sure the following ever finds more splits,
				// but good to include it for future-proofing.
				split.complement();
				if (map_.containsKey(split)) {
					count += ((Integer)map_.get(split)).intValue();
				}
				if (count > nBootStraps) {
					throw new IllegalArgumentException("Found over 100% bootstrap support");
				}
				double support = (double)count*100.0/(double)nBootStraps;
				node.setAttribute(BOOTSTRAP_ATTRIBUTE, new Double(support));
			}
		}
	}
	
	/**
	 * Print a tree with bootstrap supports.
	 * If you want more control over how it is done, use TreeUtils.printNH
	 * @param tree
	 * @param out
	 * @param printBootstrapsAsInternalLabels
	 */
	public static void printNHWithBootstraps(Tree tree, PrintWriter out, boolean printBootstrapsAsInternalLabels) {
		TreeUtils.printNH(tree, out, true, false, false, true, printBootstrapsAsInternalLabels, 0, false);
	}
}

