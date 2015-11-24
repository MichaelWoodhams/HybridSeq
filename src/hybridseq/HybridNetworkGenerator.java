package hybridseq;
import java.io.PrintWriter;
import java.util.ArrayList;

import pal.distance.DistanceMatrix;
import pal.math.MersenneTwisterFast;
import pal.misc.SimpleIdGroup;
import pal.tree.Node;
import palExtensions.ExtRandom;

/**
 * @author woodhams
 *
 */
public class HybridNetworkGenerator {
	/** Helper class RateFunction: 
	 * can be constant, linear or quadratic.
	 * constant => rate independent of number of taxa n
	 * linear => rate proportional to n (e.g. speciation rate for a Yule process)
	 * quadratic = rate proportional to pairs of taxa, n(n-1). 
	 * @author woodhams
	 *
	 */
	public enum RateFunction {
		CONST, LINEAR, QUADRATIC;
		
		public String toString() {
			return super.name().toLowerCase();
		}
		
		public static RateFunction insensitiveValueOf(String name) {
			return valueOf(name.toUpperCase());
		}
	};
	
	/**
	 * Helper class HybridFunction:
	 * the probability as a function of genetic distance that a hybridization will succeed.
	 * 
	 */
	public enum HybridFunction {
		LINEAR, STEP, QUADRATIC, EXPONENTIAL, SNOWBALL;
		
		public double probability(double threshold, double dist) {
			switch (this) {
				case LINEAR : return Math.max(0, 1-dist/threshold);
				case QUADRATIC: Math.max(0, 1-dist*dist/(threshold*threshold));
				case STEP: return (dist < threshold) ? 1 : 0;
				case EXPONENTIAL: return Math.exp(-dist/threshold);
				case SNOWBALL: return Math.exp(-dist*dist/threshold);
				default : throw new RuntimeException("Can't happen");
			}
		}
			
		public String toString() {
			return super.name().toLowerCase();
		}
			
		public static HybridFunction insensitiveValueOf(String name) {
			// strip white space and convert to upper case
			
			return valueOf(name.toUpperCase().trim());
		}
	}
	
	/**
	 * Helper class EventType
	 * 
	 * @author woodhams
	 *
	 */
	
	private enum EventType {SPECIATION, HYBRIDIZATION, INTROGRESSION, HALT, EPOCH_CHANGE};
	private class Event {
		public EventType type;
		public double interval;
		public int primaryNode = -1;
		public int secondaryNode = -1;
		
		
		public Event() {}
		
		public Event(EventType t, double d, int p, int s) {
			type = t;
			interval = d;
			primaryNode = p;
			secondaryNode = s;
		}
		
		/*
		 * Choose primaryNode and secondaryNode as distinct random leaves 
		 */
		private void randomNodes() {
			primaryNode = rng.nextInt(nLeaves);
			secondaryNode = rng.nextInt(nLeaves-1);
			if (secondaryNode >= primaryNode) secondaryNode++;
		}
		
		private double distance() {
			return distMatrix.get(primaryNode).get(secondaryNode);
		}
	}
	
	// If we have a required minimum number of successful reticulations, how many networks to make before we give up? 
	private static final int MAX_ATTEMPTS = 100; 
	private HybridNetworkParameters par;
	// Internal working variables:
	private double nextRateChangeTime;
	private ArrayList<ArrayList<Double>> distMatrix;
	private HybridNetwork network = null;
	private int nLeaves=0;
	private MersenneTwisterFast rng; 
	private ArrayList<NetworkNode> leafList = new ArrayList<NetworkNode>();
	private double currentTime = 0;
	private double speciationCurrentRate; 
	private double hybridizationBaseRate;
	private double introgressionBaseRate;
	private double currentReticulationThreshold;
	
	// These two for curiosity only.
	private int hybridizationAttempts = 0; 
	private int introgressionAttempts = 0;

/*	
	// Parameters (held in 'par'):
	StepwiseFunction<Double> speciationRateOverTime; // SPEC_RATE
	RateFunction speciationFunction; // SPEC_FUNC (linear = Yule process)
	StepwiseFunction<Double> hybridizationRateOverTime; // HYBR_RATE 
	RateFunction hybridizationFunction; // HYBR_FUNC (quadratic = proportional to number of pairs of taxa)
	StepwiseFunction<Double> introgressionRateOverTime; // INTR_RATE
	RateFunction introgressionFunction; // INTR_FUNC
	double hybridizationThreshold; // HYBR_THRESH, distance at which hybridization or introgression chance is zero. Linear function.
	// list of mix/probability pairs
	DiscreteDistribution<Double> hybridizationMixPDF; // H_MIX_PDF
	// list of mix/probability pairs
	DiscreteDistribution<Double> introgressionMixPDF; // I_MIX_PDF
	int numRandomTrees; // NUM_TREE
	int filoSitesPerTree; // FILO_PER_T: only affects output to Filo file.
	int dolloSitesPerTree; // DOLLO_PER_T: Number of Dollo characters to simulate = dolloSitesPerTree * numRandomTrees
	double dolloRate; // DOLLO_RATE: rate of loss of characters in Dollo simulation
	long seed; // RNT_SEED
	// Halting condition parameters. Simulation will halt when any one of these is reached.
	double haltTime; // HALT_TIME
	int maxHybridEvents; // MAX_HYBR
	int maxTaxa; // HALT_TAXA
*/	

	public HybridNetworkGenerator(HybridNetworkParameters params) {
		par = params;
		rng = new ExtRandom(params.seed);
	}
	
	public void setSeed(long seed) {
		rng.setSeed(seed);
	}
	
	private String nextLabel() {
		final byte A = (byte)'A';
		if (nLeaves < 26) {
			// single capital letter
			return java.lang.Character.toString((char)(A+nLeaves));
		} else {
			// capital letter plus a number
			char baseChar = (char)(A + nLeaves % 26);
			int count = nLeaves / 26;
			return String.format("%c%d", baseChar,count);
		}
	}
	
	/**
	 * Are we producing Filo sequences this run?
	 * @return
	 */
	public boolean haveFilo() {
		return par.filoSitesPerTree > 0;
	}
	
	/**
	 * Are we producing Dollo character data this run?
	 * @return
	 */
	public boolean haveDollo() {
		return par.dolloSitesPerTree > 0;
	}

	public void setHaltTime(double time) {
		par.haltTime = time;
	}
	public void setMaxHybridEvents(int max) {
		par.maxHybridEvents = max;
	}
	public void setMaxTaxa(int max) {
		par.maxTaxa = max;
	}
	
	private double getNextRateChangeTime(double t) {
		double next = Math.min(par.speciationRateOverTime.nextStep(t), par.hybridizationRateOverTime.nextStep(t));
		next        = Math.min(next, par.reticulatationThresholdOverTime.nextStep(t));
		return        Math.min(next, par.introgressionRateOverTime.nextStep(t));
	}
	private void init() {
		leafList.clear(); nLeaves=0; currentTime = 0;
		speciationCurrentRate = par.speciationRateOverTime.f(0);
		hybridizationBaseRate = par.hybridizationRateOverTime.f(0);
		introgressionBaseRate = par.introgressionRateOverTime.f(0);
		currentReticulationThreshold = par.reticulatationThresholdOverTime.f(0);
		nextRateChangeTime = getNextRateChangeTime(0);
		// Starting point: a root node with two leaf nodes, branch length zero
		SimpleNetworkNode root = new SimpleNetworkNode("root",0.0);
		SimpleNetworkNode leaf = new SimpleNetworkNode(nextLabel(),0.0);
		root.addChild(leaf);
		leafList.add(leaf);
		nLeaves++;
		leaf = new SimpleNetworkNode(nextLabel(),0.0);
		root.addChild(leaf);
		leafList.add(leaf);
		nLeaves++;	
		network = new HybridNetwork(root);
		distMatrix = new ArrayList<ArrayList<Double>>();
		distMatrix.add(new ArrayList<Double>());
		distMatrix.add(new ArrayList<Double>());
		distMatrix.get(0).add(0.0); distMatrix.get(0).add(0.0);
		distMatrix.get(1).add(0.0); distMatrix.get(1).add(0.0);
		introgressionAttempts = 0;
		hybridizationAttempts = 0;
	}
	
	
	private void implementEvent(Event event) {
		// Step one: Advance all leaf branch lengths by interval
		for (int i=0; i<nLeaves; i++) {
			Node leaf = network.getExternalNode(i);
			leaf.setBranchLength(leaf.getBranchLength()+event.interval);
			for (int j=i+1; j<nLeaves; j++) {
				double value = distMatrix.get(i).get(j);
				value += 2*event.interval;
				distMatrix.get(i).set(j, value);
				distMatrix.get(j).set(i, value);
			}
		}
		// Step two: Replace nodes as required by the event type. Will later add new node to 'newParent'.
		NetworkNode newParent = null;
		SimpleNetworkNode secondary;
		double wt = 0;
		int n = -1; // node number for the new or introgressed node
		switch (event.type) {
		case SPECIATION:
			newParent = new SimpleNetworkNode();
			insertNewNode(newParent, leafList.get(event.primaryNode));
			wt = 1;
			n = nLeaves;
			currentTime += event.interval;
			break;
		case HYBRIDIZATION: 
			HybridNetworkNode hybrid = new HybridNetworkNode();
			// Replace both leaves with an internal node which will parent the
			// replaced leaf and the hybridization node
			SimpleNetworkNode primary = new SimpleNetworkNode();
			insertNewNode(primary, leafList.get(event.primaryNode));
			secondary = new SimpleNetworkNode();
			insertNewNode(secondary, leafList.get(event.secondaryNode));
			hybrid.setCurrentParent(0);
			hybrid.setBranchLength(0.0);
			primary.addChild(hybrid);
			hybrid.setCurrentParent(1);
			hybrid.setBranchLength(0.0);
			secondary.addChild(hybrid);
			wt = par.hybridizationMixPDF.draw(rng);
			// 1-wt makes parent0 the dominant parent, as wt is expected to be <=0.5
			hybrid.setParent0Weight(1-wt);
			newParent = hybrid;
			n = nLeaves;
			currentTime += event.interval;
			break;
		case INTROGRESSION:
			HybridNetworkNode introgression = new HybridNetworkNode();
			insertNewNode(introgression, leafList.get(event.primaryNode));
			secondary = new SimpleNetworkNode();
			insertNewNode(secondary, leafList.get(event.secondaryNode));
			introgression.setCurrentParent(1);
			introgression.setBranchLength(0.0);
			secondary.addChild(introgression);
			wt = 1-par.introgressionMixPDF.draw(rng); // '1-' because parent0 is the introgressee. 
			introgression.setParent0Weight(wt);
			n = event.primaryNode;
			currentTime += event.interval;
			break;
		case HALT:
			currentTime += event.interval;
			network.setTreeDepth(currentTime);
			break;
		case EPOCH_CHANGE:
			currentTime = nextRateChangeTime;
			nextRateChangeTime = getNextRateChangeTime(currentTime);
			speciationCurrentRate = par.speciationRateOverTime.f(currentTime);
			introgressionBaseRate = par.introgressionRateOverTime.f(currentTime);
			hybridizationBaseRate = par.hybridizationRateOverTime.f(currentTime);
			currentReticulationThreshold = par.reticulatationThresholdOverTime.f(currentTime);
		}
		// Step 3: Add new leaf node (only for speciation and hybridizaiton events)
		if (newParent != null) {
			NetworkNode leaf = new SimpleNetworkNode(nextLabel(),0.0);
			newParent.addChild(leaf);
			leafList.add(leaf);
			nLeaves++;
			enlargeDistMatrix();
		}
		// Step 4: alter distance matrix
		if (event.type == EventType.SPECIATION 
				|| event.type == EventType.HYBRIDIZATION 
				|| event.type == EventType.INTROGRESSION) {
			int p = event.primaryNode;
			int s = event.secondaryNode;
			for (int i=0; i<nLeaves; i++) {
				if (i!=n) {
					double value = distMatrix.get(p).get(i)*wt + distMatrix.get(s).get(i)*(1-wt);
					distMatrix.get(i).set(n, value);
					distMatrix.get(n).set(i, value);
				}
			}
		}
		network.createNodeList();

		assert (network.getExternalNodeCount() == nLeaves);
	}
	
	/**
	 * Given a completed network, on the assumption that the network was generated by 
	 * this generator, why did it halt? 
	 */
	public String haltReason(HybridNetwork haltedNetwork) {
		StringBuffer reason = new StringBuffer();
		if (haltedNetwork.getExternalNodeCount()==par.maxTaxa) {
			reason.append("max taxa");
		}
		if (haltedNetwork.getHybridNodeCount()==par.maxHybridEvents) {
			if (reason.length()>0) reason.append(", ");
			reason.append("max hybridizations");
		}
		if ((par.haltTime-haltedNetwork.getTreeDepth())/par.haltTime < 1e-12) {
			if (reason.length()>0) reason.append(", ");
			reason.append("max time");
		}
		if (reason.length()==0) reason.append("halted for no reason");
		return reason.toString();
	}
	
	private void enlargeDistMatrix() {
		int n = distMatrix.size();
		for (int i=0; i<n; i++) distMatrix.get(i).add(0.0);
		ArrayList<Double> newRow = new ArrayList<Double>();
		newRow.ensureCapacity(n+1);
		for (int i=0; i<=n; i++) newRow.add(0.0);
		distMatrix.add(newRow);
	}
	
	/*
	 * Insert 'replacement' where 'leaf' is, such that 'replacement' has
	 * same branch length as 'leaf' did, then attach 'leaf' as a child
	 * to 'replacement' with zero branch length.
	 */
	private void insertNewNode(NetworkNode replacement, NetworkNode leaf) {
		replacement.setBranchLength(leaf.getBranchLength());
		//replacement.setNumber(-1); // keep non-negative numbers unique and only on leaves, for 'test()' routine.
		leaf.setBranchLength(0);
		Node parent = leaf.getParent();
		int childNum = 0;
		for (childNum=0; parent.getChild(childNum)!=leaf;childNum++){};
		parent.setChild(childNum, replacement);
		replacement.addChild(leaf);
	}
	
	/**
	 * Generate network with constraints on how many hybridization events it has.
	 * @param minHybrid minimum required number of hybridizations (no effect if <=0)
	 * @param maxHybrid maximum number of hybridizations - if necessary, will reduce to this number
	 *        (set to negative or very large number to have no effect.)
	 * @return Network
	 */
	public HybridNetwork generateNetwork(int minHybrid, int maxHybrid) {
		return generateNetwork(minHybrid,maxHybrid,null);
	}
	/*
	 * statsAccumulator[0] <- cumulative count of network generation attempts
	 * statsAccumulator[1] <- cumulative count of raw number of hybridizations
	 */
	public HybridNetwork generateNetwork(int minHybrid, int maxHybrid, int[] statsAccumulator) {
		HybridNetwork net;
		int count = 0;
		do {
			if (count++>MAX_ATTEMPTS) throw new RuntimeException("Couldn't reach min reticulate events in "+MAX_ATTEMPTS+" trials"); // guard against infinite loops
			net = generateNetwork();
		} while (net.getHybridNodeCount() < minHybrid);
		if (statsAccumulator!=null) {
			statsAccumulator[0]+=count;
			statsAccumulator[1]+=net.getHybridNodeCount();
		}
		if (maxHybrid >=0 && net.getHybridNodeCount() > maxHybrid) {
			// possibly should use separate RNG for this, so rng in later stages does not depend on 
			// number of nodes removed (for consistency between runs with different reduceHybridEventsTo values.)
			while (net.getHybridNodeCount() > maxHybrid) net.removeRandomHybridization(rng);
		}
		return net;
	}
	
	public HybridNetwork generateNetwork() {
		Event event = null;
		init();
		do {
			event = nextEvent();
			implementEvent(event);
		} while (event.type != EventType.HALT);
		// Convert format of distance matrix and store it in network
		network.setDistanceMatrix(convertDistMatrix());
		return network;
	}
	
	private DistanceMatrix convertDistMatrix() {
		
		SimpleIdGroup idGroup = new SimpleIdGroup(nLeaves);
		double[][] doubleDistMat = new double[nLeaves][nLeaves];
		for (int i=0; i<nLeaves; i++) {
			idGroup.setIdentifier(i, leafList.get(i).getIdentifier());
			for (int j=0; j<nLeaves; j++) {
				doubleDistMat[i][j] = distMatrix.get(i).get(j);
			}
		}
		return (new DistanceMatrix(doubleDistMat, idGroup));
	}
	
	private Event nextEvent() {
		boolean success = false;
		Event event=null;
		double accumulatedInterval = 0; // accumulate intervals previous attempts.
		while (!success) {
			event = nextEventAttempt();
			event.randomNodes();
			switch (event.type) {
			case HALT:
			case EPOCH_CHANGE:
				assert(false); // can't happen - test for this is later
			case SPECIATION:					
				success = true;
				break;
			case INTROGRESSION:
			case HYBRIDIZATION:
				success = hybridizationSuccess(event.distance());
				break;
			}
			accumulatedInterval += event.interval;
			// Now check if the proposed event would violate the max time halting condition or move past the
			// time at which rates change
			if (currentTime+accumulatedInterval > par.haltTime) {
				event.type = EventType.HALT;
				accumulatedInterval = par.haltTime - currentTime;
				success = true;
			} else if (currentTime+accumulatedInterval > nextRateChangeTime) {
				event.type = EventType.EPOCH_CHANGE;
				accumulatedInterval = nextRateChangeTime - currentTime;
				success = true;
			}
		}
		event.interval = accumulatedInterval;
		
		// Also check for halt due to number of leaves or max reticulation events
		if (network.getHybridNodeCount() == par.maxHybridEvents && 
				(event.type == EventType.HYBRIDIZATION || event.type == EventType.INTROGRESSION)) {
			event.type = EventType.HALT;
		}
		if (network.getExternalNodeCount() == par.maxTaxa && 
				(event.type == EventType.HYBRIDIZATION || event.type == EventType.SPECIATION)) {
			event.type = EventType.HALT;
			if (event.type == EventType.HYBRIDIZATION) hybridizationAttempts--; // Don't count this one.
		}
		return event;
	}
	
	private Event nextEventAttempt() {
		Event event = new Event();
		double sRate = speciationRate();
		double hRate = hybridizationRate();
		double iRate = introgressionRate();
		double rate = sRate+hRate+iRate;
		double interval = nextExponential(rate);
		event.interval = interval;
		double rand = rate*rng.nextDouble();
		if (rand < sRate) {
			event.type = EventType.SPECIATION;
		} else if (rand < sRate + iRate) {
			event.type = EventType.INTROGRESSION;
			introgressionAttempts++;
		} else {
			event.type = EventType.HYBRIDIZATION;
			hybridizationAttempts++;
		}
		return event;
	}
	
	public int getIntrogressionAttempts() {
		return introgressionAttempts;
	}
	public int getHybridizationAttempts() {
		return hybridizationAttempts;
	}
	
	private double hybridizationRate() {
		return calculateRate(hybridizationBaseRate, par.hybridizationLeafFunction);
	}
	private double introgressionRate() {
		return calculateRate(introgressionBaseRate, par.introgressionLeafFunction);
	}
	private double speciationRate() {
		return calculateRate(speciationCurrentRate, par.speciationLeafFunction);
	}
	private double calculateRate(double base, RateFunction shape) {
		double rate = 0;
		switch (shape) {
		case CONST:
			rate = base;
			break;
		case LINEAR:
			rate = base*nLeaves;
			break;
		case QUADRATIC:
			rate = base*nLeaves*(nLeaves-1);
			break;
		}
		return rate;
	}
	
	private boolean hybridizationSuccess(double dist) {
		// Linear success chance:
		return (rng.nextDouble() < par.reticFunction.probability(currentReticulationThreshold,dist));
	}
	// Should really be part of MersenneTwisterFast class
	private double nextExponential(double rate) {
		return -Math.log(rng.nextDouble())/rate;
	}
	
	// TODO: jUnit tests
	public static void test() {
		HybridNetworkGenerator netGen = new HybridNetworkGenerator(new HybridNetworkParameters());
		netGen.testGeneration();
	}
	public void testGeneration() {
		Event event;
		init();
		event = new Event(EventType.SPECIATION, 1.0, 0, 0); // A -> (A,C)
		implementEvent(event);
		event = new Event(EventType.SPECIATION, 1.0, 1, 1); // B -> (B,D)
		implementEvent(event);
		event = new Event(EventType.HYBRIDIZATION, 1.0, 1, 2); // B and C hybridize to make E
		implementEvent(event);
		event = new Event(EventType.INTROGRESSION, 1.0, 2, 4); // E introgresses into C
		implementEvent(event);
		event = new Event(EventType.INTROGRESSION, 1.0, 0, 3); // D introgresses into A
		implementEvent(event);
		event = new Event(EventType.HALT, 1.0, 0, 0);
		implementEvent(event);
		
		PrintWriter stdout = new PrintWriter(System.out);
		network.printSummary(stdout);
		stdout.close();	
	}
}
