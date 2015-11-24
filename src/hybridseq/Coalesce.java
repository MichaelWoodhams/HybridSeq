package hybridseq;
/**
 * A static methods only class.
 */

import java.util.LinkedHashSet;
import java.util.Iterator;
import java.util.Set;

import mdwUtils.StepwiseFunction;

import palExtensions.ExtRandom;
import pal.tree.Node;
import pal.tree.SimpleNode;
import pal.tree.SimpleTree;
import pal.tree.Tree;

public class Coalesce {
	
	public static Tree coalescentTree(Tree speciesTree, StepwiseFunction<Double> coalescenceRate, ExtRandom rng) {
		Set<Node> uncoalesced = uncoalescedFromNode(speciesTree.getRoot(), coalescenceRate, 0, rng);
		coalesceAlongBranch(uncoalesced, coalescenceRate, 0, Double.NEGATIVE_INFINITY, rng);
		// assert: uncoalesced now has only one member. I exercise the Axiom of Choice:
		Node root = uncoalesced.iterator().next();
		// root will have infinite branch length, so truncate this:
		root.setBranchLength(0);
		return new SimpleTree(root); 
	}
	
	
	/**
	 * Returns a set of uncoalesced nodes, representing the state of the gene
	 * tree at the top of the branch leading to this node.
	 * 
	 * @param speciesTreeNode
	 * @param topOfBranchTime  the time (distance from root) of this node's parent
	 * @param rng
	 * @return
	 */

	private static Set<Node> uncoalescedFromNode(Node speciesTreeNode, StepwiseFunction<Double> coalescenceRate, double topOfBranchTime, ExtRandom rng) {
		/*
		 * Use of LinkedHashSet here is required for replicability.
		 * Even with identical input (including RNG seeds) HashSet can vary
		 * its output order from run to run. 
		 */
		Set<Node> uncoalesced = new LinkedHashSet<Node>();
		if (speciesTreeNode.isLeaf()) {
			// TODO: worry about whether branch length needs to be set to zero
			Node newLeaf = new SimpleNode(speciesTreeNode);
			uncoalesced.add(newLeaf);
		} else {
			double nodeTime = topOfBranchTime + speciesTreeNode.getBranchLength(); // the time at this node
			for (int i=0; i<speciesTreeNode.getChildCount(); i++) {
				Node child = speciesTreeNode.getChild(i);
				// The recursion:
				uncoalesced.addAll(uncoalescedFromNode(child,coalescenceRate,nodeTime,rng));
			}
			coalesceAlongBranch(uncoalesced, coalescenceRate, nodeTime, topOfBranchTime, rng);
		}
		return uncoalesced;
	}
	
	

	/**
	 * Evolves a set of nodes backwards in time along a branch, coalescing them as it goes.
	 * 
	 * @param uncoalesced Tree nodes: initially those present at start of branch. Returns with those at end of branch
	 * @param coalescenceRate: coalescenceRate as a function of time
	 * @param bottomTime: time at bottom of branch/coalescence process
	 * @param topTime: time at top of branch
	 * 
	 * As we're going backwards in time, bottomTime >= topTime
	 */
	private static void coalesceAlongBranch(Set<Node> uncoalesced, StepwiseFunction<Double> coalescenceRate, double bottomTime, double topTime, ExtRandom rng) {
		int nNodes = uncoalesced.size();
		if (nNodes==0) throw new IllegalArgumentException("Must provide at least one node");
		if (bottomTime < topTime) throw new IllegalArgumentException("Time travel is forbidden");
		if (nNodes==1) {
			// Nothing to coalesce, so advance to top of branch and return
			advanceTime(uncoalesced, bottomTime-topTime);
			return;
		}
		/*
		 * Which event will happen 'next' (backwards in time)?
		 * + Nodes coalesce
		 * + coalescenceRate size changes
		 * + end of branch
		 */
		double coalescenceTime = bottomTime - rng.nextExponential(1/(coalescenceRate.f(bottomTime)*nNodes*(nNodes-1)/2));
		double popChangeTime = coalescenceRate.lastStep(bottomTime);
		if (topTime>coalescenceTime && topTime>popChangeTime) {
			// nothing happens until top of branch
			advanceTime(uncoalesced, bottomTime-topTime);
			return;
		}
		double currentTime = Double.NaN;
		if (popChangeTime>coalescenceTime) {
			// coalescenceRate size changes before coalescence
			advanceTime(uncoalesced, bottomTime-popChangeTime);
			currentTime = popChangeTime;
		} else {
			// Next thing to happen is a coalescence.
			// Advance to coalescence time:
			advanceTime(uncoalesced, bottomTime-coalescenceTime);
			currentTime = coalescenceTime;
			// Choose two distinct random index numbers for the two nodes to coalesce
			int a = rng.nextInt(nNodes); // a in range 0..nNodes-1
			int b = rng.nextInt(nNodes-1); // b in range 0..nNodes-2
			if (b>=a) b++; // now b is in range 0..nNodes-1 and !=a.
			Node nodeA = null; 
			Node nodeB = null;
			Iterator<Node> iter = uncoalesced.iterator();
			do {
				Node node = iter.next();
				if (a==0) nodeA=node;
				if (b==0) nodeB=node;
				a--; b--;
				if (a<0&&b<0) break;
			} while (iter.hasNext());
			// remove the nodes to coalesce
			uncoalesced.remove(nodeA);
			uncoalesced.remove(nodeB);
			// create a new parent node to replace them
			Node newNode = new SimpleNode();
			newNode.addChild(nodeA);
			newNode.addChild(nodeB);
			newNode.setBranchLength(0);
			uncoalesced.add(newNode);
		}
		// and use tail-recursion to evolve the rest of the way up the branch
		coalesceAlongBranch(uncoalesced, coalescenceRate, currentTime, topTime, rng);
	}
	
	private static void advanceTime(Set<Node> nodes, double timeInc) {
		for (Node node : nodes) node.setBranchLength(node.getBranchLength()+timeInc);
	}
}
