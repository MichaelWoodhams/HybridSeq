package hybridseq;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashSet;
import java.util.LinkedList;

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.bio.phylo.io.nexus.NexusBlock;

import pal.distance.DistanceMatrix;
import pal.math.MersenneTwisterFast;
import pal.misc.IdGroup;
import pal.misc.Identifier;
import pal.misc.SimpleIdGroup;
import pal.tree.AttributeNode;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.SimpleTree;
import pal.tree.Tree;
import pal.tree.TreeDistanceMatrix;
import pal.tree.TreeUtils;
import palExtensions.IdGroupUtils;
import biojavaExtensions.GenericBlockBuilder;
import biojavaExtensions.NexusUtils;

/**
 * HybridNetwork is to NetworkNode as SimpleTree is to SimpleNode.
 * 
 * This class was created by copying SimpleTree and then editing.
 * 
 * @author woodhams
 *
 */

/*
 *  TODO:
 *  I'm not happy with storing filoSitesPerTree, dolloSitesPerTree, numTrees as data members.
 *  They should instead be passed as parameters to getDolloAlignment etc.
 *  but currently those numbers are not easily available to where getDolloAlignment gets called.
 *  
 *  Also: consider making random number generator a data member, rather than keep passing one
 *  as a parameter to many methods.
 *  
 *  Alternatively: consider adding as data members the data HybridForests constructor gets from HybridNetworkParameters
 */


@SuppressWarnings("serial")
public class HybridNetwork implements Tree {
	// Static stuff:
	// Following limit only applies if we invoke toTrees or treeWeights, which produce output exponential in num hybr. nodes
	public static int MAX_HYBRID_NODES = 12; 
	// Attribute name used in eNewick format output routines
	private static final String ATTR_NET_NODE_NAME = "netNodeName";
	
	// Private stuff
	
	/** root node */
	private NetworkNode root;
	
	/** Total tree depth (should equal root-to-leaf sum of branch lengths for any leaf, up to rounding error) */
	private double treeDepth = Double.NEGATIVE_INFINITY; // default value to make it obvious if hasn't been set yet.

	/** list of internal nodes (including root) */
	private NetworkNode[] internalNode = null;
	/** number of internal nodes (including root) */
	private int numInternalNodes;

	/** list of hybrid nodes - subset of internalNode */
	private HybridNetworkNode[] hybridNode = null;
	/** Number of hybrid nodes */
	private int numHybridNodes;
	
	/** list of external nodes */
	private NetworkNode[] externalNode = null;
	/** number of external nodes */
	private int numExternalNodes;

	/** attributes attached to this tree. */
	private Hashtable<String,Object>[] attributes = null;

	/** holds the units of the trees branches. */
	private int units = pal.misc.Units.EXPECTED_SUBSTITUTIONS;
	
	/** Can be set by outside process (i.e. HybridNetworkGenerator.) */
	// In principle, can be calculated from the network, but I've no need for this so haven't coded it. 
	private DistanceMatrix distMatrix = null;
	
	/** constructor tree consisting solely of root node */
	public HybridNetwork() {
		// Default configuration
		root = new SimpleNetworkNode();
		root.setBranchLength(0.0);
		root.setBranchLengthSE(0.0);
	}

	/** constructor taking a root node */
	public HybridNetwork(NetworkNode r) {
		root = r;
		createNodeList();
	}
	
	public void setTreeDepth(double d) {
		treeDepth = d;
	}
	
	public double getTreeDepth() {
		return treeDepth;
	}
	
	public DistanceMatrix getDistanceMatrix() {
		return distMatrix;
	}
	
	/** Intended for use by HybridNetworkGenerator only */
	public void setDistanceMatrix(DistanceMatrix dm) {
		distMatrix = dm;
	}

	/**
	 * Resolve all the hybrid nodes at random (using the node's probabilities)
	 * and return the resulting tree.
	 * @return
	 */
	public Tree randomTree(MersenneTwisterFast rng) {
		for (int i=0; i<numHybridNodes; i++) {
			HybridNetworkNode node = hybridNode[i];
			int parent = (rng.nextDouble() < node.getParentWeight(0)) ? 0 : 1;
			node.setCurrentParent(parent);
		}
		/*
		 * If the first event was a reticulation event, then the randomly resolved tree's
		 * root will have only one child. This may be interpreted by other programs as
		 * not being a well formed tree, so in this instance we go down levels until 
		 * we find a node with multiple children.
		 */
		Node rootCandidate = rootOfResolvedTree();
		while (rootCandidate.getChildCount()==1) {
			Node newCandidate = rootCandidate.getChild(0);
			newCandidate.setBranchLength(newCandidate.getBranchLength()+rootCandidate.getBranchLength());
			newCandidate.setParent(null);
			rootCandidate = newCandidate;
		}
		return new SimpleTree(rootCandidate);
	}
		
	public void removeRandomHybridization(MersenneTwisterFast rng) {
		int nHybridNodes = getHybridNodeCount();
		if (nHybridNodes == 0) throw new RuntimeException("No hybrid node to remove");
		HybridNetworkNode oldNode = (HybridNetworkNode)getHybridNode(rng.nextInt(nHybridNodes));
		oldNode.remove();
		createNodeList(); // node lists need recalculating
	}
	
	
	/**
	 * The "most likely" tree, i.e. the one generated by the most likely parent assignments
	 * of hybrid nodes. Does not account for the fact that multiple assignments can result
	 * in the same tree.
	 * @return
	 */
	public Tree dominantTree() {
		for (int i=0; i<numHybridNodes; i++) {
			HybridNetworkNode node = hybridNode[i];
			int parent = (node.getParentWeight(0) > 0.5) ? 0 : 1;
			node.setCurrentParent(parent);
		}
		return new SimpleTree(rootOfResolvedTree());
	}
	
	/*
	 * This is new to HybridNetwork.
	 * Using the current parent allocations of all HybridNetworkNodes, return
	 * the resulting tree.
	 */
	public AttributeNode rootOfResolvedTree() {
		AttributeNode treeRoot = root.convertToAttributeNode();
		cleanConvertedTree(treeRoot,-1);
		return treeRoot;
	}
	
	/*
	 * Recursively applying NetworkNode.convertToSimpleNode() to a network can 
	 * result in various problems for the output 'tree': 
	 * We can have nodes with indegree 1 and outdegree 1 which should be removed.
	 * We can have subtrees where all the leaves are orphaned internal nodes,
	 * rather than taxa. These also should be removed.
	 * 
	 * This node's child number of its parent is passed, as we will be
	 * eliminating this node if it has indegree and outdegree == 1,
	 * and thus may need to replace the parent's child pointer.
	 */
	private static void cleanConvertedTree(AttributeNode n, int childNum) {
		// Clean all child nodes:
		for (int i=0; i<n.getChildCount(); i++) cleanConvertedTree((AttributeNode)n.getChild(i),i);
		// Check for and remove dead branches.
		for (int i=n.getChildCount()-1; i>=0; i-- ) {
			AttributeNode child = (AttributeNode)n.getChild(i);
			if (child.getAttribute(NetworkNode.LEAF_STATUS_PROPERTY)==NetworkNode.LEAF_DEAD) {
				n.removeChild(i);
			}
		}
		// If we have no live children, and are not specifically ourselves marked as alive
		// then we're dead.
		if (n.getChildCount()==0 && !(n.getAttribute(NetworkNode.LEAF_STATUS_PROPERTY)==NetworkNode.LEAF_ALIVE)) {
			n.setAttribute(NetworkNode.LEAF_STATUS_PROPERTY, NetworkNode.LEAF_DEAD);
		}
		// If we have indegree==1 and outdegree==1, we are a redundant node, so remove ourselves.
		if (n.getChildCount() == 1 && n.getParent() != null) {
			AttributeNode child = (AttributeNode)n.getChild(0);
			Node parent = n.getParent();
			child.setBranchLength(  child.getBranchLength()  +n.getBranchLength());
			child.setBranchLengthSE(child.getBranchLengthSE()+n.getBranchLengthSE());
			parent.setChild(childNum, child);
			child.setParent(parent);
		}
	}
	

	
	/**
	 * Decompose into all the trees embedded in this network. Does not eliminate duplicates.
	 * Public access should be to getWeightedForest().
	 */
	private Tree[] toTrees() {
		if (numHybridNodes > MAX_HYBRID_NODES) 
			throw new RuntimeException("Too many hybrid nodes - number of trees would be excessive");
		int nTrees = 1 << numHybridNodes; // 2^numHybridNodes
		Tree[] tree = new Tree[nTrees];
		for (int i=0; i<nTrees; i++) {
			// Set hybrid nodes according to the bitpattern of i
			int iCopy = i;
			for (int bitNum=0; bitNum<numHybridNodes; bitNum++) {
				int bit = iCopy & 1;
				iCopy >>= 1;
				hybridNode[bitNum].setCurrentParent(bit);
			}
			tree[i] = new SimpleTree(rootOfResolvedTree());
		}
		return tree;
	}
	
	/**
	 * Calculate weights of all the trees that toTrees produces.
	 * Use getWeightedForest() for public access.
	 */
	private double[] treeWeights() {
		if (numHybridNodes > MAX_HYBRID_NODES) 
			throw new RuntimeException("Too many hybrid nodes - number of trees would be excessive");
		int nTrees = 1 << numHybridNodes; // 2^numHybridNodes
		double[] weights = new double[nTrees];
		for (int i=0; i<nTrees; i++) {
			double weight = 1;
			int iCopy = i;
			for (int bitNum=0; bitNum<numHybridNodes; bitNum++) {
				int bit = iCopy & 1;
				iCopy >>= 1;
				weight *= hybridNode[bitNum].getParentWeight(bit);
			}
			weights[i] = weight;
		}
		return weights;
	}
	
	
	/*
	 * Returns WeightedForest of the unique trees which the network decomposes into.
	 * 
	 */
	
	public WeightedForest getWeightedForest() {
		Tree[] trees = this.toTrees();
		double[] weights = this.treeWeights();
		HashMap<String,WeightedTree> hash = new HashMap<String,WeightedTree>(trees.length);
		for (int i=0; i<trees.length; i++) {
			StringWriter sw = new StringWriter();
			PrintWriter pw = new PrintWriter(sw);
			NodeUtils.printNH(pw, trees[i].getRoot(), true, false, 0, false);
			String nhTree = sw.toString();
			if (hash.containsKey(nhTree)) {
				double wt = hash.get(nhTree).getWeight() + weights[i];
				hash.get(nhTree).setWeight(wt);
			} else {
				hash.put(nhTree, new SimpleWeightedTree(trees[i],weights[i]));
			}
		}
		return new WeightedForest(hash.values());
	}

	// TODO: Change how GenericBlock checks for valid field names so that
	// the try/catch below is not needed.
	public NexusBlock makeENewickBlock(String blockComment) {
		GenericBlockBuilder builder = new GenericBlockBuilder();
		builder.startBlock("eNewick");
		NexusUtils.addComment(builder, blockComment);
		String[][] eNewick = this.toENewick(true, false, false);
		try {
			for (String[] subtree : eNewick) {
				builder.addFieldPlusWhitespace(subtree[0], "=", subtree[1]);
			}
		} catch (ParseException e) {
			// this could only be triggered if the block had a valid keys list.
			e.printStackTrace();
			throw new RuntimeException("Can't happen");
		}
		builder.endBlock();
		return builder.getNexusBlock();
	}
	
	/*
	 * Primarily for debugging purposes:
	 */
	public void printSummary(PrintWriter out) {
		
		double height = 0;
		Node node = externalNode[0];
		do {
			height += node.getBranchLength();
			node = node.getParent();
		} while (!node.isRoot());
		out.printf("Tree depth is %f with %d reticulation events and %d taxa\n", height, numHybridNodes, numExternalNodes);
		out.print("\nNetwork in eNewick format:\n");
		this.printENewick(out, true, false, true);
		out.println();
		if (numHybridNodes > 4) {
			out.println("Most likely tree:");
			NodeUtils.printNH(out, dominantTree().getRoot(), true, false, 0, false);
			out.println("\nDisplay of remaining trees suppressed due to number");
		} else {
			out.println("Exact list of lineage trees:");
			WeightedForest forest = getWeightedForest();
			for (WeightedTree tree : forest) {
				out.printf("Weight %5f:  ", tree.getWeight());
				NodeUtils.printNH(out, tree.getRoot(), true, false, 0, false);
				out.println();
			}
			out.println();
		}
		
		// And the harder bit: the overall distance matrix.
		if (numExternalNodes > 15) {
			//out.println("Display of distance matrix suppressed due to size");
		} else if (numHybridNodes > 12) {
			//out.println("Display of distance matrix suppressed due to number of trees");
		} else {
			/*
			WeightedForest forest = getWeightedForest();
			int nLeaves = this.numExternalNodes;
			IdGroup leafIds = new SimpleIdGroup(nLeaves);
			for (int i=0; i<nLeaves; i++) leafIds.setIdentifier(i, externalNode[i].getIdentifier());
			leafIds = IdGroupUtils.copySorted(leafIds); 
			double[][] distFromTrees = new double[nLeaves][nLeaves];
			
			for (WeightedTree tree : forest) {
				TreeDistanceMatrix treeMat = new TreeDistanceMatrix(tree,leafIds);
				for (int j=0; j<nLeaves; j++)
					for (int k=0; k<nLeaves; k++)
						distFromTrees[j][k] += tree.getWeight()*treeMat.getDistance(j, k);
			}
			
			out.print("Distance matrix:");
			for (int i=0; i<nLeaves; i++) {
				out.printf("%8s ", leafIds.getIdentifier(i).getName());
			}
			out.println();
			for (int i=0; i<nLeaves; i++) {
				out.printf("%18s ", leafIds.getIdentifier(i).getName());
				for (int j=0; j<nLeaves; j++) {
					out.printf("%8.5f ", distFromTrees[i][j]);
				}
				out.println();
			}
			*/
		}
	}
	
	/*
	 * Primarily for debugging purposes:
	 */
	public void oldprint(PrintWriter out) {
		Tree[] trees = this.toTrees();
		double[] weights = this.treeWeights();
		int nLeaves = this.numExternalNodes;
		IdGroup leafIds = new SimpleIdGroup(nLeaves);
		for (int i=0; i<nLeaves; i++) leafIds.setIdentifier(i, externalNode[i].getIdentifier());
		leafIds = IdGroupUtils.copySorted(leafIds); 
		double[][] distFromTrees = new double[nLeaves][nLeaves];

		for (int i=0; i<trees.length; i++) {
			out.printf("Weight %5f:  ", weights[i]);
			NodeUtils.printNH(new PrintWriter(out), trees[i].getRoot(), true, false, 0, false);
			out.write(";\n");
			TreeDistanceMatrix treeMat = new TreeDistanceMatrix(trees[i],leafIds);
			for (int j=0; j<nLeaves; j++)
				for (int k=0; k<nLeaves; k++)
					distFromTrees[j][k] += weights[i]*treeMat.getDistance(j, k);
		}
		
		out.print("Distance matrix:");
		for (int i=0; i<nLeaves; i++) {
			out.printf("%8s ", leafIds.getIdentifier(i).getName());
		}
		out.println();
		for (int i=0; i<nLeaves; i++) {
			out.printf("%18s ", leafIds.getIdentifier(i).getName());
			for (int j=0; j<nLeaves; j++) {
				out.printf("%8.5f ", distFromTrees[i][j]);
			}
			out.println();
		}
	}

	/**
	 * Return the units that this tree is expressed in.
	 */
	public final int getUnits() {
		return units;
	}

	/**
	 * Sets the units that this tree is expressed in.
	 */
	public final void setUnits(int units) {
		this.units = units;
	}
		
	/**
	 * Returns the number of external nodes.
	 */
	public final int getExternalNodeCount() {
		if(externalNode==null) {
			createNodeList();
		}
		return numExternalNodes;
	}

	/**
	 * Returns the ith external node.
	 */
	public final Node getExternalNode(int i) {
		if(externalNode==null) {
			createNodeList();
		}
		return externalNode[i];
	}

	/**
	 * Returns the number of internal nodes.
	 */
	public final int getInternalNodeCount() {
		if(internalNode==null) {
			createNodeList();
		}
		return numInternalNodes;
	}

	/**
	 * Returns the ith internal node.
	 */
	public final Node getInternalNode(int i) {
		if(internalNode==null) {
			createNodeList();
		}
		return internalNode[i];
	}
	
	/**
	 * Returns the number of hybrid nodes.
	 */
	public final int getHybridNodeCount() {
		if(hybridNode==null) {
			createNodeList();
		}
		return numHybridNodes;
	}

	/**
	 * Returns the ith hybrid node.
	 */
	public final Node getHybridNode(int i) {
		if(hybridNode==null) {
			createNodeList();
		}
		return hybridNode[i];
	}

	/**
	 * Returns the root node of this tree.
	 */
	public final Node getRoot() {
		return root;
	}

	/**
	 * Set a new node as root node.
	 */
	public final void setRoot(NetworkNode r) {
		root = r;
		createNodeList();
	}
	

	/** count and list external, internal and hybrid nodes and
		compute heights of each node */
	public void createNodeList()
	{
		/*
		 * This needed a complete rewrite compared to SimpleTree, as now
		 * we have a network and may be able to get to a given node by
		 * several roots. We do a depth-first-search, putting found nodes into
		 * a set, and skipping any nodes we reach which have been seen before.
		 */
		numInternalNodes = 0;
		numExternalNodes = 0;
		numHybridNodes = 0;
		/*
		 *  LinkedHashSet vs HashSet is required, otherwise order of nodes (and hence
		 *  results of downstream calculations) can vary from run to run, with identical
		 *  start conditions (input parameters, RNG seeds.)
		 */
		LinkedHashSet<NetworkNode> processed = new LinkedHashSet<NetworkNode>();
		LinkedList<NetworkNode> unchecked = new LinkedList<NetworkNode>();
		
		unchecked.add(root);

		while (!unchecked.isEmpty()) {
			NetworkNode node = unchecked.pop();
			if (!processed.contains(node)) {
				if (node.isLeaf()) {
					numExternalNodes++;
				} else {
					numInternalNodes++;
				}
				if (node instanceof HybridNetworkNode) {
					numHybridNodes++;
				}
				for (int i=0; i<node.getChildCount(); i++) {
					unchecked.push((NetworkNode)node.getChild(i));
				}
				processed.add(node);
			}
		} 
		
		

		internalNode = new NetworkNode[numInternalNodes];
		externalNode = new NetworkNode[numExternalNodes];
		hybridNode   = new HybridNetworkNode[numHybridNodes];
		int iIntern=0;
		int iExtern=0;
		int iHybrid=0;
		for(NetworkNode node : processed) {
			if (node.isLeaf()) {
				node.setNumber(iExtern);
				externalNode[iExtern++]=node;
			} else {
				node.setNumber(iIntern);
				internalNode[iIntern++]=node;
			}
			if (node instanceof HybridNetworkNode) {
				hybridNode[iHybrid++]=(HybridNetworkNode)node;
			}
		}

		// NOT DONE: Calculating node heights.
	}

	// Will convert to NH tree format on the basis of current parent assignments of HybridNetworkNodes
	public String toString() {
		StringWriter sw = new StringWriter();
		NodeUtils.printNH(new PrintWriter(sw), rootOfResolvedTree(), true, false, 0, false);
		sw.write(";");
		
		return sw.toString();
	}


	/**
	 * return node with number num (as stored in the 'number' field)
	 *
	 * @param num number of node
	 *
	 * @return node
	 */
	public Node findNode(int num)
	{
		createNodeList();

		if (num <= numExternalNodes)
		{
			return externalNode[num-1];
		}
		else
		{
			return internalNode[num-1-numExternalNodes];
		}
	}

	private int getIndex(Node node) {
		if (node.isLeaf()) return node.getNumber();
		return getExternalNodeCount() + node.getNumber();
	}

	/**
	 * Sets a named attribute for a given node.
	 * @param node the node whose attribute is being set.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	@SuppressWarnings("unchecked")
	public void setAttribute(Node node, String name, Object value) {
		if (node instanceof AttributeNode) {
			((AttributeNode)node).setAttribute(name, value);
		} else {
			int index = getIndex(node);
			if (attributes == null) {
				/* This line requires the @SuppressWarnings
				 * http://stackoverflow.com/questions/14917375/cannot-create-generic-array-of-how-to-create-an-array-of-mapstring-obje
				 * gives an alternative method, but it still requires @SuppressWarnings
				 */
				attributes = new Hashtable[getExternalNodeCount() + getInternalNodeCount()];
			}
			if (attributes[index] == null) {
				attributes[index] = new Hashtable<String,Object>();
			}
			attributes[index].put(name, value);
		}
	}

// ========= IdGroup stuff ===============================
	public int getIdCount() {
		return getExternalNodeCount();
	}
	public Identifier getIdentifier(int i) {
		return getExternalNode(i).getIdentifier();
	}
	public void setIdentifier(int i, Identifier id) {
		getExternalNode(i).setIdentifier(id);
	}
	public int whichIdNumber(String s) {
		return IdGroup.Utils.whichIdNumber(this,s);
	}

//========================================================
	/**
	 * @return an object representing the named attributed for the numbered node.
	 * @param node the node being interrogated.
	 * @param name the name of the attribute of interest.
	 */
	public Object getAttribute(Node node, String name) {
		if (node instanceof AttributeNode) {
			return ((AttributeNode)node).getAttribute(name);
		} else {
			int index = getIndex(node);
			if (attributes == null || attributes[index] == null) {
				return null;
			}
			return attributes[index].get(name);
		}
	}


	// interface Report

	public void report(PrintWriter out)
	{
		TreeUtils.report(this, out);
	}

	public Tree getCopy() {
		return new SimpleTree(this);
	}
	
	/**
	 * Set a new node as root node.
	 * Required by Tree interface, but liable to result in time travel
	 * if implemented on a HybridNetwork, so diabled:
	 */
	public final void setRoot(Node r) {
		throw new UnsupportedOperationException("Rerooting a network is a bad idea.");
	}
	
	// Build a test hybrid network, and decompose it into trees.
	public static void test() {
		// First, create all the nodes:
		// Leaves at time 6
		SimpleNetworkNode leafA = new SimpleNetworkNode("A",1.0);
		SimpleNetworkNode leafB = new SimpleNetworkNode("B",2.0);
		SimpleNetworkNode leafC = new SimpleNetworkNode("C",2.0);
		SimpleNetworkNode leafD = new SimpleNetworkNode("D",3.0);
		SimpleNetworkNode leafE = new SimpleNetworkNode("E",1.0);
		// Introgression from E ancestor to A ancestor at time 5
		SimpleNetworkNode introgressorE = new SimpleNetworkNode("int to A",3.0);
		HybridNetworkNode introgressionAE = new HybridNetworkNode("int A from E");
		// Introgression from C ancestor to B ancestor at time 4
		SimpleNetworkNode introgressorC = new SimpleNetworkNode("int to B",1.0);
		HybridNetworkNode introgressionBC = new HybridNetworkNode("int B from C");
		// Hybridization from B, D ancestors to form C ancestor at time 3
		SimpleNetworkNode hybridizorB = new SimpleNetworkNode("B hyb to C",2.0);
		SimpleNetworkNode hybridizorD = new SimpleNetworkNode("D hyb to C",1.0);
		HybridNetworkNode hybridBD = new HybridNetworkNode("hybrid BD");
		// Speciation of D and E at time 2
		SimpleNetworkNode specDE = new SimpleNetworkNode("DE",2.0);
		// speciation of A and B at time 1
		SimpleNetworkNode specAB = new SimpleNetworkNode("AB",1.0);
		// Root at time 0
		SimpleNetworkNode root = new SimpleNetworkNode("root",0.0);
		
		// Now join them all together
		// Introgression from E ancestor to A ancestor at time 5
		introgressionAE.setCurrentParent(1);
		introgressionAE.setBranchLength(0);
		introgressionAE.setParent0Weight(0.9);
		introgressionAE.addChild(leafA);
		introgressorE.addChild(leafE);		
		introgressorE.addChild(introgressionAE);
		introgressionAE.setCurrentParent(0);
		introgressionAE.setBranchLength(4.0); // to AB parent node
		// Introgression from C ancestor to B ancestor at time 4
		introgressionBC.setCurrentParent(1);
		introgressionBC.setBranchLength(0.0); 
		introgressionBC.setParent0Weight(0.8);
		introgressionBC.addChild(leafB);
		introgressorC.addChild(leafC);
		introgressorC.addChild(introgressionBC);
		introgressionBC.setCurrentParent(0);
		introgressionBC.setBranchLength(1.0);
		
		// Hybridization from B, D ancestors to form C ancestor at time 3
		hybridBD.setCurrentParent(0);
		hybridBD.setBranchLength(0.0);
		hybridBD.setParent0Weight(0.6);
		hybridBD.addChild(introgressorC);
		hybridizorB.addChild(introgressionBC);
		hybridizorB.addChild(hybridBD);
		hybridBD.setCurrentParent(1);
		hybridBD.setBranchLength(0.0);
		hybridizorD.addChild(leafD);	
		hybridizorD.addChild(hybridBD);
		// Speciation of D and E at time 2
		specDE.addChild(hybridizorD);
		specDE.addChild(introgressorE);
		// speciation of A and B at time 1
		specAB.addChild(introgressionAE);
		specAB.addChild(hybridizorB);
		// Root at time 0
		root.addChild(specAB);
		root.addChild(specDE);

		HybridNetwork testNetwork = new HybridNetwork(root);
		
		Tree[] trees = testNetwork.toTrees();
		double[] weights = testNetwork.treeWeights();
		
		PrintWriter stdout = new PrintWriter(System.out);
		for (int i=0; i<trees.length; i++) {
			stdout.printf("Weight %5f:  ", weights[i]);
			NodeUtils.printNH(new PrintWriter(stdout), trees[i].getRoot(), true, false, 0, false);
			stdout.write(";\n");
		}
		stdout.close();
		
	}
	
	public void printENewick(PrintWriter out, boolean printLengths, boolean printInternalLabels, boolean breakLines) {
		String[][] eNewick = this.toENewick(printLengths, printInternalLabels, breakLines);
		for (String[] subtree : eNewick) {
			out.printf("\t%s=%s;\n", subtree[0], subtree[1]);
		}
	}
	
	public String[][] toENewick(boolean printLengths, boolean printInternalLabels, boolean breakLines) {
		int n=this.getHybridNodeCount();
		String[][] subtrees = new String[n+1][2];
		try {
			StringWriter sw = new StringWriter();
			// Ensure each hybrid node has a unique name. (If the node actually has an identifier, that will take precedence.)
			for (int i=0; i<n; i++) {
				String label = String.format("hybr%d", i);
				((AttributeNode)this.getHybridNode(i)).setAttribute(ATTR_NET_NODE_NAME, label);
			}
			subtrees[0][0] = "root";
			printENewickNode(sw, (AttributeNode)this.getRoot(), null, printLengths, printInternalLabels, breakLines, 5);
			subtrees[0][1] = sw.toString();
			for (int i=0; i<n; i++) {
				AttributeNode aNode = (AttributeNode)this.getHybridNode(i);
				String label = aNode.getIdentifier().toString();
				if (label.length()==0) label = (String)aNode.getAttribute(ATTR_NET_NODE_NAME);
				subtrees[i+1][0] = label;
				sw.getBuffer().setLength(0);
				printENewickNode(sw, aNode, null, printLengths, printInternalLabels, breakLines, label.length()+1);
				subtrees[i+1][1] = sw.toString();
			}
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException("Can't happen");
		}
		return subtrees;
	}
	
	/*
	 * Mildly edited version of pal.tree.NodeUtils.printNH
	 * 
	 * There are four main cases:
	 * 1 'node' is a leaf. Print its identifier and possibly branch length.
	 * 2 'node' is the root, or is a hybrid node acting as the root of its own subtree. This is signaled by parentNode being
	 *   null. In this case, print its children but not its name or branch length.
	 * 3 'node' is an internal non-hybrid node. Recursively print its children, then optionally name and branch length
	 *   (this has subcases of leaf and non-leaf node.)
	 * 4 'node' is a hybrid node which has been reached by recursion from this method. In this case, 
	 *   print its name plus weighting but descend no further.
	 */
	private static int printENewickNode(Writer out, AttributeNode node, Node parentNode,
			boolean printLengths, boolean printInternalLabels, 
			boolean breakLines, int column) throws IOException {

		if (breakLines) column = breakLine(out, column);
		int nodeType=0;
		if (node.isLeaf()) nodeType = 1;
		else if (parentNode == null) nodeType = 2;
		else if (node instanceof HybridNetworkNode) nodeType = 4;
		else nodeType = 3;
		
		if (nodeType==2 || nodeType==3) {
			// Recurse into children if we have ordinary internal node, or hybrid node being treated as root
			out.write("(");
			column++;
			for (int i = 0; i < node.getChildCount(); i++) {
				if (i != 0) {
					out.write(",");
					column++;
				}
				column = printENewickNode(out, (AttributeNode)(node.getChild(i)), node,  printLengths, printInternalLabels, breakLines, column);
			}
			out.write(")");
			column++;
		}
		if (nodeType==4) {
			// Hybrid node reached by recursion. Print name, branch length, weight but do not recurse.
			// (we need special code to find the correct name and branch length.)
			HybridNetworkNode hNode = (HybridNetworkNode)node;
			if (breakLines) column = breakLine(out, column);
			String id = hNode.getIdentifier().toString();
			if (id.length()==0) id = (String)hNode.getAttribute(ATTR_NET_NODE_NAME);
			out.write(id);
			column += id.length();
			hNode.setCurrentParentTo(parentNode);
			if (printLengths) {
				out.write(":");
				column++;
				if (breakLines) column = breakLine(out, column);
				String decimal = String.format("%7f", hNode.getBranchLength());
				out.write(decimal);
				column += decimal.length();
			}
			if (breakLines) column = breakLine(out, column);
			String wt = String.format("[wt=%f]", hNode.getParentWeight());
			out.write(wt);
			column += wt.length();
		} 	
		if (nodeType==1 || (nodeType==3 && printInternalLabels)) {
			// Print label for leaf or ordinary internal node 
			if (breakLines) column = breakLine(out, column);
			String id = node.getIdentifier().toString();
			out.write(id);
			column += id.length();
		}
		if (printLengths && (nodeType==1 || nodeType==3)) {
			// branch lengths for leaves and ordinary internal nodes
			out.write(":");
			column++;
			if (breakLines) column = breakLine(out, column);
			String decimal = String.format("%7f", node.getBranchLength());
			out.write(decimal);
			column += decimal.length();
		}
		return column;
	}
	
	// Direct copy from pal.tree.NodeUtils
	private static int breakLine(Writer out, int column) throws IOException {
		if (column > 70) {
			out.write('\n');
			column = 0;
		}
		return column;
	}
}
