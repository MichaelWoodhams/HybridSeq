package hybridseq;

import pal.tree.Tree;

/**
 * A tree with a weight.
 * 
 * @author woodhams
 */

public interface WeightedTree extends Tree {
	void setWeight(double wt);
	double getWeight();
}
