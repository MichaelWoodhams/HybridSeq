package hybridseq;

import java.util.HashMap;
import java.util.Iterator;

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.bio.phylo.io.nexus.CharactersBlockBuilder;
import org.biojavax.bio.phylo.io.nexus.NexusBlock;
import org.biojavax.bio.phylo.io.nexus.NexusComment;
import org.biojavax.bio.phylo.io.nexus.NexusFile;
import org.biojavax.bio.phylo.io.nexus.NexusFileBuilder;
import org.biojavax.bio.phylo.io.nexus.NexusObject;

import biojavaExtensions.GenericBlock;
import biojavaExtensions.GenericBlockBuilder;
import biojavaExtensions.GenericNexusStuff;
import biojavaExtensions.NexusUtils;

import pal.alignment.Alignment;
import pal.tree.Tree;
import palExtensions.ExtRandom;

public class HybridForests {
	private HybridNetwork net;
	private ExtRandom rng;
	private final int nTree;
	private Tree[] lineageTrees;
	private Tree[] coalescentTrees;
	private String[] lineageTreesAsStrings;
	private String[] coalescentTreesAsStrings;
	private Alignment dolloAlignment;
	private HybridNetworkParameters params;
	
	public HybridForests(HybridNetwork net, HybridNetworkParameters params) {
		this(net,params,params.seed);
	}
	public HybridForests(HybridNetwork net, HybridNetworkParameters params, long seed) {
		rng = new ExtRandom(seed+1);
		this.net = net;
		this.params = params;
		nTree = params.numRandomTrees; // redundant, but used often enough to cache it
		lineageTrees = new Tree[nTree];
		lineageTreesAsStrings = new String[nTree];
		makeRandomLineageTrees();
		if (params.coalesce) {
			coalescentTrees = new Tree[nTree];
			coalescentTreesAsStrings = new String[nTree];
			for (int i=0; i<nTree; i++) {
				coalescentTrees[i] = Coalesce.coalescentTree(lineageTrees[i], params.coalescenceRate, rng);
				coalescentTreesAsStrings[i] = coalescentTrees[i].toString();
			}
		}
		if (params.dolloSitesPerTree>0) {
			DolloProcess dollo = new DolloProcess(net, params.dolloRate, rng);
			if (params.coalesce) {
				dolloAlignment = dollo.getDolloAlignment(coalescentTrees,params.dolloSitesPerTree);
			} else {
				dolloAlignment = dollo.getDolloAlignment(lineageTrees,params.dolloSitesPerTree);
			}
		}
	}

    /**
     * 	
     * @param rng
     */
	private void makeRandomLineageTrees() {
		// This method is complicated by wanting to collect together cases where lineage tree is the same.
		HashMap<String,Tree> stringToTree = new HashMap<String,Tree>();
		for (int i=0; i<nTree; i++) {
			Tree tree = net.randomTree(rng);
			String str = tree.toString();
			stringToTree.put(str, tree);
			lineageTreesAsStrings[i] = str;
		}
		java.util.Arrays.sort(lineageTreesAsStrings);
		// TODO: (possibly) scan through array and make identical strings into the same string object.
		for (int i=0; i<nTree; i++) {
			lineageTrees[i] = stringToTree.get(lineageTreesAsStrings[i]);
		}
	}
	
	public Tree[] getLineageTrees() {
		return lineageTrees;
	}
	public Tree[] getCoalescentTrees() {
		return coalescentTrees;
	}
	
	public NexusBlock getTaxaBlock() {
		int nTaxa = net.getIdCount();
		String[] taxa = new String[nTaxa];
		for (int i=0; i<nTaxa; i++) taxa[i] = net.getIdentifier(i).toString();
		return NexusUtils.makeTaxaBlock(taxa);
	}
	
	// parameters nexus block is params.getNexusParametersBlock()
	// eNewick block is net.makeENewickBlock(null);
	
	public NexusBlock getLineageTreesBlock() {
		return NexusUtils.makeTreesBlock(
			lineageTreesAsStrings,
			null,
			"Randomly selected lineage (i.e. no coalescent) trees",
			"RLT"
		);
	}
	
	public NexusBlock getCoalescentTreesBlock() {
		return NexusUtils.makeTreesBlock(
			coalescentTreesAsStrings, 
			lineageTreesAsStrings, 
			"Randomly selected coalescent trees (with generating lineage trees as comments)",
			"CT"
		);
	}
	
	public NexusBlock getExactLineageTreesBlock() {
		return net.getWeightedForest().toTreesBlock("Lineage trees with exact probabilities","ELT");
	}
	
	/**
	 * Return Filo block generated from template Filo block found in the input Nexus file
	 * @param inputNexusFile
	 * @return
	 */
	public NexusBlock getFiloBlock(NexusFile inputNexusFile, String[] trees, int sitesPerTree) {
		NexusBlock filoBlock=null;
		@SuppressWarnings("unchecked")
		Iterator<NexusBlock> blockIter = inputNexusFile.blockIterator();
		while (blockIter.hasNext()) {
			NexusBlock block = blockIter.next();
			if (block instanceof GenericBlock && block.getBlockName().equalsIgnoreCase("filo")) {
				filoBlock = getFiloBlock((GenericBlock)block, trees, sitesPerTree);
				break; // Should there be more than one Filo block, ignore the rest.
			}
		}
		return filoBlock;
	}
	
	public NexusBlock getFiloBlock(GenericBlock filoTemplate, String[] trees, int sitesPerTree) {
		GenericBlock filoBlock=null;
		try {
			filoBlock = modifyFiloBlock(
					filoTemplate, 
					trees, 
					sitesPerTree
			);
		} catch (ParseException e) {
			System.err.println("Parse error modifying Filo block - can't happen");
			throw new RuntimeException(e);
		}
		return filoBlock;
	}
	
	/*
	 * Some of this code may be misplaced: why do we write the parameters block here? That 
	 * doesn't belong to the forest. 
	 * TODO: consider this. Maybe make this method only add some Trees blocks to an existing
	 * NexusFile.
	 */
	public NexusFile makeNexusFile(NexusFile inputNexusFile, NexusComment[] prefixComments) {
		NexusFileBuilder outputBuilder = new NexusFileBuilder();
		outputBuilder.startFile();
		NexusFile outputNexusFile = outputBuilder.getNexusFile();

		for (NexusComment comment : prefixComments) {
			outputNexusFile.addObject(comment);
		}
		
		/* 
		 * Filo has a bug: it can generate parse errors on Nexus blocks
		 * which are NOT Filo blocks but occur before the Filo block.
		 * In particular, the HybridSim block causes a "Number format failure" error message. 
		 * To avoid this, we need to ensure that the Filo block is the first block.
		 */
		if (params.filoSitesPerTree>0) {
			String[] trees = (coalescentTreesAsStrings==null) ? lineageTreesAsStrings : coalescentTreesAsStrings;
			outputNexusFile.addObject(getFiloBlock(inputNexusFile, trees, params.filoSitesPerTree));
		}
		outputNexusFile.addObject(params.getNexusParametersBlock()); // HybridSim parameters block
		outputNexusFile.addObject(getTaxaBlock());
		if (params.dolloSitesPerTree>0) {
			if (params.coalesce) {
				outputNexusFile.addObject(alignmentToBlock(dolloAlignment,"Characters evolved on coalescent trees"));
			} else {
				outputNexusFile.addObject(alignmentToBlock(dolloAlignment,"Characters evolved on lineage trees"));
			}
		}

		if (net.getHybridNodeCount() <= HybridNetwork.MAX_HYBRID_NODES) { 
			outputNexusFile.addObject(getExactLineageTreesBlock());
		} else {
			outputNexusFile.addObject(NexusUtils.newNexusComment("Exact lineage trees omitted - there are too many of them."));
		}
		outputNexusFile.addObject(getLineageTreesBlock());
		if (params.coalesce) {
			outputNexusFile.addObject(getCoalescentTreesBlock());
		}
		outputNexusFile.addObject(net.makeENewickBlock(null));
		
		// Now any (non-auto-generated) comments and unrecognized Nexus blocks in input Nexus file
		// will pass through unmodified to the output Nexus file.
		int deleteComments = 2; // delete two comments if they are first two items and they contain 'HybridSim'.
		// This is to account for the auto-added HybridSim comments when using an output file as new input.
		@SuppressWarnings("unchecked")
		Iterator<NexusObject>objIter = inputNexusFile.objectIterator();
		while (objIter.hasNext()) {
			NexusObject obj = objIter.next();
			if (obj instanceof NexusComment) {
				NexusComment comment = (NexusComment)obj;
				// Pass comment through unless it contains HybridSim
				// comment.toString does not do what it should, hence extra complication
				if (deleteComments>0 && NexusUtils.toString(comment).contains("HybridSim")) {
					// do not pass comment, decrement count
					deleteComments--;
				} else {
					outputNexusFile.addObject(obj);
					deleteComments=0;
				}
			} else {
				deleteComments=0; // only delete comments prior to first block
				// object is a block
				NexusBlock block = (NexusBlock)obj;
				String name = block.getBlockName();
				if (!name.equalsIgnoreCase("filo") &&
					!name.equalsIgnoreCase("hybridSim") && 
				    !name.equalsIgnoreCase("characters") &&
				    !name.equalsIgnoreCase("taxa") &&
				    !name.equalsIgnoreCase("trees") &&
				    !name.equalsIgnoreCase("eNewick") &&
				    !name.equalsIgnoreCase("distance")
				) {
					/*
					 *  The block is not one which HybridSim needs to replace,
					 *  so just replicate it unchanged
					 *  EXCEPT: BioJava bug: read in an 'unknown' block and write it out
					 *  inserts a newline after the 'begin', so need to fix that,
					 *  via UseableNexusBlock class
					 */
					outputNexusFile.addObject(block);
				}
			}
		}
		return outputNexusFile;
	}
	
	private static NexusBlock alignmentToBlock(Alignment align, String blockComment) {
		CharactersBlockBuilder builder = new CharactersBlockBuilder();
		builder.startBlock("characters");
		NexusUtils.addComment(builder, "Characters generated by a Dollo parsimonious process");
		NexusUtils.addComment(builder, blockComment);
		builder.setDimensionsNTax(align.getIdCount());
		builder.setDimensionsNChar(align.getSiteCount());
		builder.setDataType("standard");
		builder.addSymbol("0");
		builder.addSymbol("1");
		builder.setGap("-"); // Biojava bug: gap must be initialized.
		
		for (int i=0; i<align.getIdCount(); i++) {
			builder.addMatrixEntry(align.getIdentifier(i).toString());
			builder.appendMatrixData(align.getIdentifier(i).toString(), align.getAlignedSequenceString(i));
		}
		builder.endBlock();
		return builder.getNexusBlock();
	}

	
	private static GenericBlock modifyFiloBlock(GenericBlock inputBlock, String[] trees, int sitesPerTree) throws ParseException {
		/*
		 * Take some stuff from the input filo block and copy it through to the output,
		 * other stuff needs to be replaced:
		 * 'tree' fields in input are ignored and new ones generated.
		 * 'treeparams' fields in input are ignored and new ones generated.
		 * 'run' command in input is ignored and new 'run' added at end of the block
		 * 'params' field in input gets modified to change the overall sequence length.
		 * Anything else passes through unchanged. 
		 */
		GenericBlockBuilder builder = new GenericBlockBuilder();
		builder.startBlock("Filo");
		Iterator<GenericNexusStuff> iter = inputBlock.stuffIterator();
		while (iter.hasNext()) {
			GenericNexusStuff entry = iter.next();
			// write entry to builder unless it is one that needs overwriting: params, tree, treeparams or run
			if (entry.isComment()) {
				builder.addStuff(entry);
			} else if (entry.isWhite()) {
				// do nothing: we'll reformat stuff to add required whitespace
			} else if (entry.getKey().equalsIgnoreCase("params")) {
				// Delete length specification (if any) and add our own
				String newValue = String.format("\t\tl %d\n", trees.length*sitesPerTree)+entry.getValue().replaceAll("\\s+l +\\d+\\n", "");
				builder.addFieldPlusWhitespace("params", "\n", newValue);
			} else if (entry.getKey().equalsIgnoreCase("run") || entry.getKey().matches("tree.*")) {
				// do nothing: delete these entries
			} else {
				// pass through, after reformatting
				entry.reformat();
				builder.addWhiteSpace("\t");
				builder.addStuff(entry);
				builder.addWhiteSpace("\n");
			}
		}
		/*
		 *  Follow up with writing new values: It is complicated a bit because we 
		 *  want to amalgamate identical trees (which, because of the 'sort' in
		 *  HybridNetwork.makeRandomTrees(), will be adjacent to each other.)
		 */
		String lastTree = null;
		int siteCount = 0;
		int treeCount = 0;
		for (String tree : trees) {
			if (tree.equals(lastTree)) {
				siteCount+=sitesPerTree;
			} else {
				addTreeToFiloBlock(builder,lastTree,siteCount,treeCount); // checks for lastTree==null
				lastTree = tree;
				siteCount=sitesPerTree;
				treeCount++;
			}
		}
		addTreeToFiloBlock(builder,lastTree,siteCount,treeCount);
		builder.addFieldPlusWhitespace("run", "", "");
		builder.endBlock();
		return (GenericBlock)builder.getNexusBlock();
	}
	
	private static void addTreeToFiloBlock(GenericBlockBuilder builder, String tree, int siteCount, int treeCount) throws ParseException {
		if (tree != null) {
			String treeName = String.format("t%d", treeCount);
			String trimmedTree = tree.replace(";", ""); // else get two ';': one from tree string, one from Nexus field termination
			builder.addFieldPlusWhitespace("tree "+treeName, "=", trimmedTree);
			builder.addFieldPlusWhitespace("treeparams "+treeName, "\n", String.format("\t\tl %d\n\t",siteCount));
		}
	}

}
