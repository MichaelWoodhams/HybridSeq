package hybridseq;
import pal.tree.AttributeNode;
import pal.tree.Node;
import pal.tree.SimpleNode;
/**
 * A node in a hybrid network. 
 * 
 * Subclass HybridNetworkNode accommodates hybridization events.
 * Such nodes can be set as to which of two parents is the 'current'
 * one. Then this class's convertToSimpleNode returns a SimpleNode
 * which takes just the current parental relationships, and has
 * recursively converted the child nodes of this one also.
 * 
 * In addition: a SimpleNode returned by convertToSimpleNode might
 * be a leaf because the NetworkNode was a leaf, or it might be a 
 * a leaf because all its children are HybridNetworkNodes which
 * currently do not consider it to be their parent. In the later
 * case, we set a property "leaf status" to be "dead", else
 * we set it to "alive".
 * 
 * @author woodhams
 *
 */

@SuppressWarnings("serial")
public abstract class NetworkNode extends SimpleNode {
	protected static String LEAF_STATUS_PROPERTY = "leaf status";
	protected static String LEAF_ALIVE = "alive";
	protected static String LEAF_DEAD = "dead";

	/**
	 * Returns an AttributeNode (actually a SimpleNode) 
	 * @return
	 */
	AttributeNode convertToAttributeNode() {
		SimpleNode newNode = new SimpleNode();
		newNode.setIdentifier(getIdentifier());
		newNode.setBranchLength(getBranchLength());
		newNode.setBranchLengthSE(getBranchLengthSE());
		newNode.setNodeHeight(getNodeHeight());
		if (this.isLeaf()) {
			newNode.setAttribute(LEAF_STATUS_PROPERTY, LEAF_ALIVE);
		} else {
			int legitimateChildCount = 0;
			int n = getChildCount();
			for (int i=0; i<n; i++) {
				if (getChild(i).getParent() == this) {
					// Child believes we are its parent.
					NetworkNode childNode = (NetworkNode) getChild(i); // cast should never fail.
					AttributeNode simpleChild = childNode.convertToAttributeNode();
					legitimateChildCount++;
					newNode.addChild(simpleChild);
				}
			}
			if (legitimateChildCount == 0) {
				newNode.setAttribute(LEAF_STATUS_PROPERTY, LEAF_DEAD);
			}
		}
		return newNode;
	}
	
	void removeChild(Node child) {
		int n = getChildCount();
		int i;
		for (i=0; i<n && getChild(i)!=child; i++) {}
		if (i==n) throw new IllegalArgumentException("Child to remove is not a child of this node");
		removeChild(i);
	}
	
	/**
	 * A do-nothing method for now, but is overridden in HybridNetworkNode
	 * @param parent
	 */
	void setCurrentParentTo(Node parent) {
		if (this.getParent()!=parent) throw new IllegalArgumentException ("Tried to set current parent to non-parent node");
	}
}
