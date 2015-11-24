package hybridseq;
import pal.misc.Identifier;
import pal.tree.Node;

/**
 * 
 * A network node representing a hybridization event: either introgression 
 * (gene flow from another node into this one) or the genesis of a hybrid 
 * species (two other nodes contribute the genes for the child of this one.)
 * 
 * @author woodhams
 *
 */
/*
 * WARNING! The arrays altParent, altLength and altLengthSE are not always up to date.
 * This is due to PAL's SimpleNode.setParent() etc being final. As such,
 * use the getParent(), getParent(n) methods to find parents, not access to altParent[]
 * and similar for the other arrays.
 * 
 * TODO: (possibly) make NetworkNode not extend SimpleNode to avoid this.
 */
@SuppressWarnings("serial")
public class HybridNetworkNode extends NetworkNode {
	private int currentParent = 0; // zero or one.
	/* The natural thing to do would be to make parent, length and lengthSE into
	 * length 2 arrays, but some of SimpleNode's methods are 'final', which 
	 * precludes this approach. So we keep length 2 arrays and copy
	 * values to 'parent', 'length', 'lengthSE' when we change
	 * the current parent.
	 */
	private Node[] altParent;
	private double[] altLength;
	private double[] altLengthSE;
	private double parent0Weight; // how much of genome comes from parent 0 (between 0 to 1)
	
	public HybridNetworkNode() {
		super();
		currentParent = 0;
		altParent = new Node[]{null,null};
		altLength = new double[]{0.0,0.0};
		altLengthSE = new double[]{0.0,0.0};
	}
	
	public HybridNetworkNode(String name) {
		this();
		setIdentifier(new Identifier(name));
	}
	
	/**
	 * Remove this hybridization node from a network, patching the network
	 * to look tree-like in it absence. (Essentially the same effect as
	 * setting parent0Weight = 1) 
	 */
	public void remove() {
		assert (this.getChildCount()==1);
		NetworkNode child = (NetworkNode)this.getChild(0);
		// can't access altparent[] directly as contents may be stale (due to workaround of 'final' methods in PAL.)
		NetworkNode p0 = (NetworkNode)this.getParent(0); // cast to NetworkNode for the 'removeChild(Node)' method.
		NetworkNode p1 = (NetworkNode)this.getParent(1);
		
		// need to ensure we have the appropriate parent set as 'current', in case child is a HybridNetworkNode
		child.setCurrentParentTo(this);
		child.setBranchLength(child.getBranchLength()+altLength[0]);
		p0.removeChild(this);
		p0.addChild(child);
		
		// now set child to be p1's child which is not this node
		child = null;
		for (int i=0; i<p1.getChildCount(); i++) {
			if (p1.getChild(i)!=this) child = (NetworkNode)p1.getChild(i);
		}
		child.setCurrentParentTo(p1);
		NetworkNode p1p = (NetworkNode)p1.getParent(); // parent 1's parent
		child.setBranchLength(child.getBranchLength()+p1.getBranchLength());
		p1p.removeChild(p1);
		p1p.addChild(child);
	}
	
	// get/set methods for which parent is current
	
	public void setCurrentParent(int p) {
		// Nothing to do if we're not changing currentParent
		if (p == currentParent) return;
		if (p==0 || p==1) {
			// We are changing current parent: Store current values
			// into the arrays, restore the alternate values
			int o = currentParent; // old
			altParent[o] = getParent();
			altLength[o] = getBranchLength();
			altLengthSE[o] = getBranchLengthSE();
			setParent(altParent[p]);
			setBranchLength(altLength[p]);
			setBranchLengthSE(altLengthSE[p]);
			currentParent = p;			
		} else {
			throw new IllegalArgumentException("Current parent must be 0 or 1");
		}
	}
	public int getCurrentParent() {
		return currentParent;
	}
	
	/**
	 * 'parent' must be one of the already set two alternate parents.
	 * Sets the node to consider that node to be the current parent.
	 * @param parent
	 */
	@Override
	public void setCurrentParentTo(Node parent) {
		if (this.getParent()==parent) return; // nothing to do
		if (altParent[1-currentParent]==parent) {
			this.setCurrentParent(1-currentParent);
			return;
		}
		throw new IllegalArgumentException("Supposed parent node was not a parent node");
	}
	
	public void setParent0Weight(double wt) {
		if (wt < 0 || wt > 1) throw new IllegalArgumentException("Must have 0 <= wt <= 1");
		parent0Weight = wt;
	}
	
	public double getParentWeight() {
		return getParentWeight(currentParent);
	}
	public double getParentWeight(int p) {
		return (p == 0) ? parent0Weight : 1.0 - parent0Weight;
	}
	
	// 'set' methods which are aware there are two parents
	
	public void setParent(int p, Node n) {
		altParent[p] = n;
		if (p==currentParent) setParent(n);
	}
	
	public void removeParent(int p) {
		altParent[p] = null;
		if (p==currentParent) removeParent();
	}

	public void setBranchLength(int p, double len) {
		altLength[p] = len;
		if (p==currentParent) setBranchLength(len);
	}
	
	public void setBranchLengthSE(int p, double len) {
		altLengthSE[p] = len;
		if (p==currentParent) setBranchLengthSE(len);
	}

	public Node getParent(int p) {
		if (p==currentParent) altParent[p] = getParent();
		return altParent[p];
	}
	
}
