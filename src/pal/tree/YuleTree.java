package pal.tree;

import pal.math.MersenneTwisterFast;
import pal.misc.IdGroup;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.util.AlgorithmCallback;

/**
 * A very quick random tree generator, based on the Yule process.
 *
 * I ignore the AlgorithmCallback. Nothing in here should take very long anyway.
 * 
 * @author Michael Woodhams
 *
 *
 * TODO:
 * Alternative constructors taking array of Strings, IdGroup, or just number of leaves.
 */

public class YuleTree implements TreeGenerator {
	private Identifier[] leaves_;
	private MersenneTwisterFast rng_;
	private int nLeaves_;
	private SimpleNode[] leafList_; // working space: nonpersistent.
	private double speciationRate_;
	
	public YuleTree(IdGroup leafIDs, double rate, long seed) {
		this(leafIDs,rate);
		rng_.setSeed(seed);
	}
	
	public YuleTree(IdGroup leafIDs, double rate) {
		nLeaves_ = leafIDs.getIdCount();
		leaves_ = new Identifier[nLeaves_];
		for (int i=0; i<nLeaves_; i++) {
			leaves_[i] = leafIDs.getIdentifier(i);
		}
		rng_ = new MersenneTwisterFast(); // seeds from system time.
		leafList_ = new SimpleNode[nLeaves_]; // preallocate storage.
		speciationRate_ = rate;
	}
	
	public YuleTree(IdGroup leafIDs) {
		this(leafIDs,1.0);
	}
	
	public YuleTree(String[] leafNames, double rate, long seed) {
		this(new SimpleIdGroup(leafNames), rate, seed);
	}
	
	public Tree getNextTree(AlgorithmCallback callback) { 
		SimpleNode root = new SimpleNode();
		leafList_[0] = root;
		int size = 1; // current number of leaves in leafList.
		while (size < nLeaves_) {
			int toBranch = rng_.nextInt(size);
			SimpleNode parent = leafList_[toBranch];
			SimpleNode child1 = new SimpleNode();
			SimpleNode child2 = new SimpleNode();
			parent.addChild(child1);
			parent.addChild(child2);
			leafList_[toBranch] = child1; // replaces parent, which is no longer a leaf.
			leafList_[size++]   = child2;
			double x = rng_.nextExponential(1/(speciationRate_*size));
			for (int i=0; i<size; i++) {
				leafList_[i].setBranchLength(leafList_[i].getBranchLength()+x);
			}
		}
		rng_.shuffle(leaves_);
		for (int i=0; i<nLeaves_; i++) {
			leafList_[i].setIdentifier(leaves_[i]);
		}
		return new SimpleTree(root);
	}

}
