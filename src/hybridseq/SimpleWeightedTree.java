package hybridseq;

import pal.tree.SimpleTree;
import pal.tree.Tree;

@SuppressWarnings("serial")
public class SimpleWeightedTree extends SimpleTree implements WeightedTree {
	private double weight;

	public SimpleWeightedTree(Tree tree, double weight) {
		super(tree);
		this.weight = weight;
	}
	public void setWeight(double wt) {
		weight = wt;
	}

	public double getWeight() {
		return weight;
	}

}
