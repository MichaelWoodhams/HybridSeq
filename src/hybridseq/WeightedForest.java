package hybridseq;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.RandomAccess;

import org.biojavax.bio.phylo.io.nexus.TreesBlock;

import biojavaExtensions.NexusUtils;

/**
 * For now, this is just a wrapper around ArrayList<WeightedTree>. 
 * Defined as a class for ease of later enhancements.
 * @author woodhams
 *
 */

/*
 * Alternative: 
 * 
 * public class WeightedForest extends ArrayList<WeightedTree>{}
 * 
 * Much shorter, but I'm not sure I want WeightedForest to *be* an ArrayList,
 * although I could just promise myself never to count on it being one.
 * 
 * Also:
 * Consider implementing List<WeightedTree>.
 */

public class WeightedForest implements Iterable<WeightedTree>, Collection<WeightedTree>, RandomAccess {
	private ArrayList<WeightedTree> list;
	
	public WeightedForest() {
		list = new ArrayList<WeightedTree>();
	}
	
	public WeightedForest(Collection<? extends WeightedTree> collection) {
		list = new ArrayList<WeightedTree>(collection);
	}
	
	public TreesBlock toTreesBlock(String blockComment) {
		return toTreesBlock(blockComment,"T");
	}
	public TreesBlock toTreesBlock(String blockComment, String treePrefix) {
		int nTrees = list.size();
		String[] treeStrings = new String[nTrees];
		String[] weightStrings = new String[nTrees];
		for (int i=0; i<nTrees; i++) {
			treeStrings[i] = list.get(i).toString();
			weightStrings[i] = String.format("weight=%06f", list.get(i).getWeight()); 
		}			
		return NexusUtils.makeTreesBlock(treeStrings, weightStrings, blockComment, treePrefix);
	}
	
	/*
	 * From here on, we do nothing except pass methods through to 'list'.
	 */
	
	public void add(int n, WeightedTree tree) {
		list.add(n,tree);
	}

	public WeightedTree get(int n) {
		return list.get(n);
	}


	public boolean add(WeightedTree tree) {
		return list.add(tree);
	}

	public boolean addAll(Collection<? extends WeightedTree> trees) {
		return list.addAll(trees);
	}

	public void clear() {
		list.clear();
	}

	public boolean contains(Object tree) {
		return list.contains(tree);
	}

	public boolean containsAll(Collection<?> trees) {
		return list.containsAll(trees);
	}

	public boolean isEmpty() {
		return list.isEmpty();
	}

	public boolean remove(Object tree) {
		return list.remove(tree);
	}

	public boolean removeAll(Collection<?> trees) {
		return list.removeAll(trees);
	}

	public boolean retainAll(Collection<?> trees) {
		return list.retainAll(trees);
	}

	public int size() {
		return list.size();
	}

	public Object[] toArray() {
		return list.toArray();
	}

	public <T> T[] toArray(T[] trees) {
		return list.toArray(trees);
	}

	public Iterator<WeightedTree> iterator() {
		return list.iterator();
	}
}
