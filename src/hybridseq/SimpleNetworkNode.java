package hybridseq;
import pal.misc.Identifier;



@SuppressWarnings("serial")
public class SimpleNetworkNode extends NetworkNode {
	public SimpleNetworkNode(Identifier id, double branchLength) {
		super();
		setIdentifier(id);
		setBranchLength(branchLength);
	}
	public SimpleNetworkNode(String name, double branchLength) {
		this(new Identifier(name),branchLength);
	}
	public SimpleNetworkNode() {
		super();
	}

}
