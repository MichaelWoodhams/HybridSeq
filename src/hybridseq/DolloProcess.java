package hybridseq;

import pal.alignment.Alignment;
import pal.alignment.SimpleAlignment;
import pal.datatype.TwoStates;
import pal.math.MersenneTwisterFast;
import pal.misc.IdGroup;
import pal.tree.AttributeNode;
import pal.tree.Tree;
import palExtensions.IdGroupUtils;

public class DolloProcess {
	protected MersenneTwisterFast rng;
	protected double lossRate;
	protected IdGroup leafIds;

	public DolloProcess(IdGroup ids, double rate, MersenneTwisterFast rng) {
		this.rng = rng;
		lossRate = rate;
		leafIds = IdGroupUtils.copySorted(ids);
	}
	
	public Alignment getDolloAlignment(Tree[] trees, int nSitesPerTree) {
		int nTrees = trees.length;
		if (nSitesPerTree*nTrees <= 0) return null;
		int nLeafs = leafIds.getIdCount();
		char[][] characters = new char[nLeafs][nTrees*nSitesPerTree];
		char[] sitePattern;
		for (int treeIndex=0; treeIndex<nTrees; treeIndex++) {
			Tree tree = trees[treeIndex];
			for (int site=0; site<nSitesPerTree; site++) {
				sitePattern = simulateCharacter(tree, leafIds);
				for (int leaf=0; leaf<nLeafs; leaf++) {
					characters[leaf][treeIndex*nSitesPerTree+site] = sitePattern[leaf];
				}
			}
		}
		SimpleAlignment align = new SimpleAlignment(leafIds, characters, TwoStates.DEFAULT_INSTANCE);
		return align;
	}
	
	protected static final String PRESENT = "Dollo character present";
	protected void characterDescent(AttributeNode node) {
		// Character was present at this node's parent. Did it reach this node?
		if (rng.nextDouble() > Math.exp(-lossRate*node.getBranchLength())) {
			// Was lost, so nothing more to do
			return;
		}
		// Else was not lost
		if (node.isLeaf()) {
			node.setAttribute(PRESENT, true);
		} else {
			for (int i=0; i<node.getChildCount(); i++) {
				characterDescent((AttributeNode)node.getChild(i));
			}
		}
	}
	
	public char[] simulateCharacter(Tree tree, IdGroup sortedIds) {
		// set PRESENT to 'false' on all leaves
		int nLeaves = tree.getExternalNodeCount();
		for (int i=0; i<nLeaves; i++) {
			AttributeNode node = (AttributeNode)(tree.getExternalNode(i));
			node.setAttribute(PRESENT, false);
		}
		// Select the highest node at which the character is present
		DiscreteDistribution<AttributeNode> nodeDistribution = calculateNodeDistribution(tree);
		AttributeNode birthNode;
		boolean birthAccepted = false;
		do {
			birthNode = nodeDistribution.draw(rng);
			if (birthNode.isRoot()) {
				birthAccepted = true;
			} else {
				/*
				 * Given that a birth event occurred on the branch above this node (uniform
				 * distribution of where on the branch the event occurred) and edge above node
				 * has branch length l, the probability that the character will not have been
				 * lost again before reaching the node is (1-exp(-r*b))/(r*b) where r = loss rate
				 */
				double x = lossRate*birthNode.getBranchLength();
				birthAccepted = (rng.nextDouble() < (1-Math.exp(-x))/x);
			}
		} while (!birthAccepted);
		// recurse through the tree
		characterDescent(birthNode);
		// and examine which nodes the character managed to reach. For consistency, we must
		// search the nodes in sorted order.
		char[] charPresence = new char[nLeaves];
		for (int i=0; i<nLeaves; i++) {
			int nodeNumber = tree.whichIdNumber(sortedIds.getIdentifier(i).toString());
			Object present = ((AttributeNode)tree.getExternalNode(nodeNumber)).getAttribute(PRESENT);
			charPresence[i] = ((Boolean)present) ? '1' : '0';
		}
		return charPresence;
	}
	
	protected DiscreteDistribution<AttributeNode> calculateNodeDistribution(Tree tree) {
		int nLeaves = tree.getExternalNodeCount();
		int nInternal = tree.getInternalNodeCount();
		AttributeNode[] nodes = new AttributeNode[nLeaves+nInternal];
		double[] weights = new double[nLeaves+nInternal];
		for (int i=0; i<nInternal; i++) {
			AttributeNode node = (AttributeNode) tree.getInternalNode(i);
			nodes[i] = node;
			if (node.isRoot()) {
				/*
				 * In combination with always accepting a birth event on the root node,
				 * this value gives the correct ratio of birth events prior to root to
				 * birth events within the tree.
				 */
				weights[i] = 1/lossRate; 
			} else {
				weights[i] = node.getBranchLength();
			}
		}
		for (int i=0; i<nLeaves; i++) {
			AttributeNode node = (AttributeNode) tree.getExternalNode(i);
			nodes[nInternal+i] =  node;
			weights[nInternal+i] = node.getBranchLength();
		}
		return new DiscreteDistribution<AttributeNode>(nodes,weights);
	}
}
