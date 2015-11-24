package hybridseq;

import hybridseq.HybridNetworkGenerator.HybridFunction;
import hybridseq.HybridNetworkGenerator.RateFunction;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import mdwUtils.StepwiseFunction;

import org.biojavax.bio.phylo.io.nexus.NexusBlock;
import org.biojavax.bio.phylo.io.nexus.NexusFile;

import biojavaExtensions.GenericBlock;
import biojavaExtensions.GenericBlockBuilder;
import biojavaExtensions.GenericNexusStuff;

// TODO: some sort of small data structure per parameter. Will greatly compress this code.
// Possibly: for each datatype, have a HashMap which maps param name to its value.
public class HybridNetworkParameters {
	
	// The parameters:
	public double[] epochs; // EPOCHS
	public StepwiseFunction<Double> speciationRateOverTime; // SPEC_RATE
	public RateFunction speciationLeafFunction; // SPEC_LEAF (linear)
	public StepwiseFunction<Double> hybridizationRateOverTime; // HYBR_RATE 
	public RateFunction hybridizationLeafFunction; // HYBR_LEAF (quadratic)
	public StepwiseFunction<Double> introgressionRateOverTime; // INTR_RATE
	public RateFunction introgressionLeafFunction;  // INTR_LEAF (quadratic)
	// HYBR_THRESH, scale over which hybridization probability falls with distance.
	// Except for exponential decay, hybridizationThreshold = minimum distance at which hybridization probability = 0. 
	public StepwiseFunction<Double> reticulatationThresholdOverTime;  // HYBR_THRESH
	public HybridFunction reticFunction; // shape of hybridization probability decline with distance.
	public StepwiseFunction<Double> coalescenceRate; // COAL_RATE
	public boolean coalesce; // COAL. if true, use coalescent trees for tree output and filo/dollo sequences
	// list of mix/probability pairs
	public DiscreteDistribution<Double> hybridizationMixPDF = 
			new DiscreteDistribution<Double>(new Double[]{0.5}, new double[]{1}); // H_MIX_PDF
	// list of mix/probability pairs
	public DiscreteDistribution<Double> introgressionMixPDF =
			new DiscreteDistribution<Double>(new Double[]{0.5}, new double[]{1}); // I_MIX_PDF
	public int numRandomTrees; // NUM_TREE
	public int filoSitesPerTree; // FILO_PER_T: only affects output to Filo file.
	public int dolloSitesPerTree; // DOLLO_PER_T: Number of Dollo characters to simulate = dolloSitesPerTree * numRandomTrees
	public double dolloRate; // DOLLO_RATE: rate of loss of characters in Dollo simulation
	public long seed; // RNG_SEED
	public int minReticEvents; // if successful hybridizations fewer than this, rerun simulation. No effect if negative.
	public int reduceReticEventsTo; // if more than this many hybrid events, randomly remove some. No effect if negative.
	// Halting condition parameters. Simulation will halt when any one of these is reached.
	public double haltTime; // HALT_TIME
	public int maxHybridEvents; // HALT_HYBR
	public int maxTaxa; // HALT_TAXA
	
	// Field names for parameters
	public final static String BLOCK_NAME  = "HybridSim"; 
	public final static String EPOCHS      = "epochs";
	public final static String SPEC_RATE   = "speciation rate";
	public final static String HYBR_RATE   = "hybridization rate";
	public final static String INTR_RATE   = "introgression rate";
	public final static String COAL_RATE   = "coalescence rate";
	public final static String COAL_TIME   = "coalescence time"; // inverse of coalescence rate.
	public final static String COAL        = "coalesce";
	public final static String DOLLO_RATE  = "dollo rate";
	public final static String RETIC_THRESH= "reticulation threshold";
	public final static String RETIC_FUNC  = "reticulation function"; 
	public final static String HALT_TIME   = "halt time";
	public final static String SPEC_LEAF   = "speciation leaf function";
	public final static String HYBR_LEAF   = "hybridization leaf function";
	public final static String INTR_LEAF   = "introgression leaf function";
	public final static String HALT_TAXA   = "halt taxa";
	public final static String HALT_RETIC  = "halt reticulations";
	public final static String NUM_TREE    = "number random trees";
	public final static String FILO_PER_T  = "filo sites per tree";
	public final static String DOLLO_PER_T = "dollo sites per tree";
	public final static String RNG_SEED    = "seed";
	public final static String H_MIX_PDF   = "hybridization distribution";
	public final static String I_MIX_PDF   = "introgression distribution";
	public final static String MIN_RETIC   = "minimum reticulations";
	public final static String RED_RETIC   = "reduce reticulations to";
	/*
	 * HYBRIDSIM_VALID_KEYS is used by MainConfiguration, to detect
	 * invalid key names early in parsing process. This could be skipped
	 * with little harm (errors would just occur later in setParameter.)
	 */
	public final static Set<String> HYBRIDSIM_VALID_KEYS = new HashSet<String>(); 
	static {
		HYBRIDSIM_VALID_KEYS.add(EPOCHS);
		HYBRIDSIM_VALID_KEYS.add(SPEC_RATE);
		HYBRIDSIM_VALID_KEYS.add(HYBR_RATE);
		HYBRIDSIM_VALID_KEYS.add(INTR_RATE);
		HYBRIDSIM_VALID_KEYS.add(COAL_RATE);
		HYBRIDSIM_VALID_KEYS.add(COAL_TIME);
		HYBRIDSIM_VALID_KEYS.add(COAL);
		HYBRIDSIM_VALID_KEYS.add(DOLLO_RATE);
		HYBRIDSIM_VALID_KEYS.add(RETIC_THRESH);
		HYBRIDSIM_VALID_KEYS.add(RETIC_FUNC);
		HYBRIDSIM_VALID_KEYS.add(HALT_TIME);
		HYBRIDSIM_VALID_KEYS.add(SPEC_LEAF);
		HYBRIDSIM_VALID_KEYS.add(HYBR_LEAF);
		HYBRIDSIM_VALID_KEYS.add(INTR_LEAF);
		HYBRIDSIM_VALID_KEYS.add(HALT_TAXA);
		HYBRIDSIM_VALID_KEYS.add(HALT_RETIC);
		HYBRIDSIM_VALID_KEYS.add(NUM_TREE);
		HYBRIDSIM_VALID_KEYS.add(FILO_PER_T);
		HYBRIDSIM_VALID_KEYS.add(DOLLO_PER_T);
		HYBRIDSIM_VALID_KEYS.add(RNG_SEED);
		HYBRIDSIM_VALID_KEYS.add(H_MIX_PDF);
		HYBRIDSIM_VALID_KEYS.add(I_MIX_PDF);
		HYBRIDSIM_VALID_KEYS.add(MIN_RETIC);
		HYBRIDSIM_VALID_KEYS.add(RED_RETIC);
	}
	
	
	public HybridNetworkParameters(NexusFile nexusFile) {
		this(); // set defaults
		// Find the HybridSim nexus block:
		@SuppressWarnings("unchecked")
		Iterator<NexusBlock> iter = (Iterator<NexusBlock>)nexusFile.blockIterator();
		GenericBlock block = null;
		while (block == null && iter.hasNext()) {
			NexusBlock nextBlock = iter.next();
			if (nextBlock.getBlockName().equalsIgnoreCase(BLOCK_NAME)) {
				block = (GenericBlock)nextBlock;
				break;
			}
		}
		if (block == null) throw new RuntimeException("HybridSim block missing from Nexus file");
		// Find and set epochs first, to avoid error messages
		setParameter(EPOCHS,block.getValueTrimmed(EPOCHS));
		// Then loop over the rest of the assignments
		Iterator<GenericNexusStuff> stuffIter = block.stuffIterator();
		while (stuffIter.hasNext()) {
			GenericNexusStuff stuff = stuffIter.next();
			if (stuff.isField() && !stuff.getKey().equals(EPOCHS)) {
				setParameter(stuff.getKey(),stuff.getValue());
			}
		}
	}
	
	/**
	 * Set all default values
	 */
	private static final double[] DEFAULT_EPOCHS = new double[]{};
	private static final StepwiseFunction<Double> DEFAULT_STEPWISE_FUNCTION = new StepwiseFunction<Double>(DEFAULT_EPOCHS, new Double[]{1.0});
	private static final DiscreteDistribution<Double> DEFAULT_HYBRID_MIX = new DiscreteDistribution<Double>(new Double[]{0.5}, new double[]{1});
	private static final DiscreteDistribution<Double> DEFAULT_INTROG_MIX = new DiscreteDistribution<Double>(new Double[]{0.1}, new double[]{1});
	public HybridNetworkParameters() {
		maxTaxa = 20;  
		maxHybridEvents = 20;
		numRandomTrees = 1;
		filoSitesPerTree = 0;
		dolloSitesPerTree = 0;
		minReticEvents = 0;
		reduceReticEventsTo = -1; // -1 => no reduction in hybrid nodes
		seed = 4; // http://xkcd.com/221/
		dolloRate = 1;
		haltTime = 10;
		coalesce = true;
		speciationRateOverTime         = DEFAULT_STEPWISE_FUNCTION;
		hybridizationRateOverTime      = DEFAULT_STEPWISE_FUNCTION;
		introgressionRateOverTime      = DEFAULT_STEPWISE_FUNCTION;
		reticulatationThresholdOverTime= DEFAULT_STEPWISE_FUNCTION;
		coalescenceRate                = DEFAULT_STEPWISE_FUNCTION;
		speciationLeafFunction    =  RateFunction.LINEAR;    // Yule process
		hybridizationLeafFunction =  RateFunction.QUADRATIC; // (quadratic = proportional to number of pairs of taxa)
		introgressionLeafFunction =  RateFunction.QUADRATIC;
		epochs = DEFAULT_EPOCHS;
		reticFunction = HybridFunction.LINEAR;
		hybridizationMixPDF = DEFAULT_HYBRID_MIX;
		introgressionMixPDF = DEFAULT_INTROG_MIX;
	}
	
	public void setParameter(String name, String value) {
		// perhaps convert 'name' to lowercase?
		value = value.trim();
		switch (name) {
		case HALT_TAXA   : maxTaxa              = Integer.valueOf(value); break;
		case HALT_RETIC  : maxHybridEvents      = Integer.valueOf(value); break;
		case NUM_TREE    : numRandomTrees       = Integer.valueOf(value); break;
		case FILO_PER_T  : filoSitesPerTree     = Integer.valueOf(value); break;
		case DOLLO_PER_T : dolloSitesPerTree    = Integer.valueOf(value); break;
		case MIN_RETIC   : minReticEvents       = Integer.valueOf(value); break;
		case RED_RETIC   : reduceReticEventsTo  = Integer.valueOf(value); break;
		case RNG_SEED    : seed                 =    Long.valueOf(value); break;
		case DOLLO_RATE  : dolloRate            =  Double.valueOf(value); break;
		case HALT_TIME   : haltTime             =  Double.valueOf(value); break;
		case COAL        : coalesce             = Boolean.valueOf(value); break;
		case SPEC_RATE   : speciationRateOverTime         = new StepwiseFunction<Double>(epochs, parseList_Double(value)); break;
		case HYBR_RATE   : hybridizationRateOverTime      = new StepwiseFunction<Double>(epochs, parseList_Double(value)); break;
		case INTR_RATE   : introgressionRateOverTime      = new StepwiseFunction<Double>(epochs, parseList_Double(value)); break;
		case RETIC_THRESH: reticulatationThresholdOverTime= new StepwiseFunction<Double>(epochs, parseList_Double(value)); break;
		case COAL_RATE   : coalescenceRate                = new StepwiseFunction<Double>(epochs, parseList_Double(value)); break;
		case COAL_TIME   : coalescenceRate                =     invertedStepwiseFunction(epochs, parseList_Double(value)); break;
		case H_MIX_PDF   : hybridizationMixPDF = DiscreteDistribution.valueOf(value); 
			if (!hybridizationMixPDF.valuesInRange(0.0,1.0)) throw new IllegalArgumentException(H_MIX_PDF+": values must be in [0,1]");
			break;
		case I_MIX_PDF   : introgressionMixPDF = DiscreteDistribution.valueOf(value); 
			if (!introgressionMixPDF.valuesInRange(0.0,1.0)) throw new IllegalArgumentException(I_MIX_PDF+": values must be in [0,1]");
			break;
		case RETIC_FUNC  : reticFunction            = HybridFunction.insensitiveValueOf(value); break;
		case SPEC_LEAF   : speciationLeafFunction    =  RateFunction.insensitiveValueOf(value); break;
		case HYBR_LEAF   : hybridizationLeafFunction =  RateFunction.insensitiveValueOf(value); break;
		case INTR_LEAF   : introgressionLeafFunction =  RateFunction.insensitiveValueOf(value); break;
		case EPOCHS      :
			// Special case: it is an error to set epoch-dependent parameters, and then change the epochs.
			epochs = parseList_double(value); 
			if (   speciationRateOverTime         != DEFAULT_STEPWISE_FUNCTION
				|| hybridizationRateOverTime      != DEFAULT_STEPWISE_FUNCTION
				|| introgressionRateOverTime      != DEFAULT_STEPWISE_FUNCTION
				|| reticulatationThresholdOverTime!= DEFAULT_STEPWISE_FUNCTION
				|| coalescenceRate                != DEFAULT_STEPWISE_FUNCTION) {
				// Not really an IllegalArgument, but can't be bothered making new exception class just for this
				throw new IllegalArgumentException("Can't change epochs after defining epoch-dependent parameters");
			}
			break;
		default     : throw new IllegalArgumentException("Could not parse value '"+value+"' for parameter '"+name+"'");
		}
	}
	
	private static StepwiseFunction<Double> invertedStepwiseFunction(double[] epochs, Double[] inverseValues) {
		Double[] revertedValues = new Double[inverseValues.length];
		for (int i=0; i<inverseValues.length; i++) {
			revertedValues[i]=1/inverseValues[i];
		}
		return new StepwiseFunction<Double>(epochs, revertedValues);
	}
				
	/**
	 * Write all of the parameters out to a Nexus block.
	 * EXCEPTION:
	 * the constant/linear/quadratic rate function dependencies are unlikely to change, so only print 
	 * them if they differ from default.
	 * @return
	 */
	public NexusBlock getNexusParametersBlock() {
		GenericBlockBuilder builder = new GenericBlockBuilder();
		builder.startBlock(BLOCK_NAME); 
		boolean printFunctionalForms = (speciationLeafFunction != RateFunction.LINEAR
				                  || hybridizationLeafFunction != RateFunction.QUADRATIC
				                  || introgressionLeafFunction != RateFunction.QUADRATIC);
		boolean printHybridLimits = (reduceReticEventsTo>=0 || minReticEvents>=0);
		try {
			builder.addFieldPlusWhitespace(EPOCHS, "=", mdwUtils.Strings.arrayToString(epochs));
			builder.addFieldPlusWhitespace(SPEC_RATE,  "=", speciationRateOverTime.getValuesAsString());
			if (printFunctionalForms) builder.addFieldPlusWhitespace(SPEC_LEAF,  "=", speciationLeafFunction.toString());
			builder.addFieldPlusWhitespace(RETIC_THRESH,"=", reticulatationThresholdOverTime.getValuesAsString());
			builder.addFieldPlusWhitespace(RETIC_FUNC,  "=", reticFunction.toString());
			builder.addFieldPlusWhitespace(HYBR_RATE,  "=", hybridizationRateOverTime.getValuesAsString());
			if (printFunctionalForms) builder.addFieldPlusWhitespace(HYBR_LEAF,  "=", hybridizationLeafFunction.toString());
			builder.addFieldPlusWhitespace(H_MIX_PDF,  "=", hybridizationMixPDF.toString());
			builder.addFieldPlusWhitespace(INTR_RATE,  "=", introgressionRateOverTime.getValuesAsString());
			if (printFunctionalForms) builder.addFieldPlusWhitespace(INTR_LEAF,  "=", introgressionLeafFunction.toString());
			builder.addFieldPlusWhitespace(INTR_RATE,  "=", introgressionRateOverTime.getValuesAsString());
			builder.addFieldPlusWhitespace(I_MIX_PDF,  "=", introgressionMixPDF.toString());
			builder.addFieldPlusWhitespace(COAL,       "=", Boolean.toString(coalesce));
			if (coalesce) builder.addFieldPlusWhitespace(COAL_RATE,  "=", coalescenceRate.getValuesAsString());
			if (printHybridLimits) builder.addFieldPlusWhitespace(MIN_RETIC, "=", Integer.toString(minReticEvents));
			if (printHybridLimits) builder.addFieldPlusWhitespace(RED_RETIC, "=", Integer.toString(reduceReticEventsTo));
			builder.addFieldPlusWhitespace(HALT_TIME,  "=", Double.toString(haltTime));
			builder.addFieldPlusWhitespace(HALT_TAXA,  "=", Integer.toString(maxTaxa));
			builder.addFieldPlusWhitespace(HALT_RETIC,  "=", Integer.toString(maxHybridEvents));
			builder.addFieldPlusWhitespace(NUM_TREE,   "=", Integer.toString(numRandomTrees));
			builder.addFieldPlusWhitespace(FILO_PER_T, "=", Integer.toString(filoSitesPerTree));
			builder.addFieldPlusWhitespace(DOLLO_PER_T,"=", Integer.toString(dolloSitesPerTree));
			builder.addFieldPlusWhitespace(DOLLO_RATE, "=", Double.toString(dolloRate));
			builder.addFieldPlusWhitespace(RNG_SEED,   "=", Long.toString(seed));
		} catch (org.biojava.bio.seq.io.ParseException e) {
			// Should never happen
			throw new RuntimeException("Field names need updating in writeNexusParametersBlock");
		}
		builder.endBlock();
		return builder.getNexusBlock();
	}
	
	public static String getParameterDatatype(String paramName) {
		switch (paramName) {
		case HALT_TAXA   :
		case HALT_RETIC   :
		case NUM_TREE    :
		case FILO_PER_T  :
		case DOLLO_PER_T :
		case MIN_RETIC   :
		case RED_RETIC   : return "int";
		case RNG_SEED    : return "long";
		case DOLLO_RATE  : 
		case HALT_TIME   : return "double";
		case COAL        : return "boolean";
		case SPEC_RATE   : 
		case HYBR_RATE   : 
		case INTR_RATE   : 
		case RETIC_THRESH: 
		case COAL_RATE   : 
		case COAL_TIME   : return "StepwiseFunction"; // all StepwiseFunctions are StepwiseFunction<Double>
		case H_MIX_PDF   : 
		case I_MIX_PDF   : return "DiscreteDistribution";
		case RETIC_FUNC  : return "HybridFunction";
		case SPEC_LEAF   : 
		case HYBR_LEAF   :
		case INTR_LEAF   : return "RateFunction";
		default : throw new IllegalArgumentException("'"+paramName+"' is not a legal parameter name");
		}
	}
	
	private static double[] parseList_double(String line) {
		// Remove any surrounding parentheses (which are optional)
		line = line.replaceAll("[()]", "");
		if (line.matches("^\\s*$")) return new double[0]; // number list contained no numbers, just whitespace.
		String delimiter = ",";
		String[] tokens = line.split(delimiter);
		int n = tokens.length;
		double[] list = new double[n];
		for (int i=0; i<n; i++) {
			list[i] = Double.valueOf(tokens[i]);
		}
		return list;
	}

	private static Double[] parseList_Double(String line) {
		// Remove any surrounding parentheses (which are optional)
		line = line.replaceAll("[()]", "");
		if (line.matches("^\\s*$")) return new Double[0]; // number list contained no numbers, just whitespace.
		String delimiter = ",";
		String[] tokens = line.split(delimiter);
		int n = tokens.length;
		Double[] list = new Double[n];
		for (int i=0; i<n; i++) {
			list[i] = Double.valueOf(tokens[i]);
		}
		return list;
	}
}
