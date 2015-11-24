package hybridseq;

import mdwUtils.StepwiseFunction;
import pal.math.MersenneTwisterFast;


public class DiscreteDistribution<T> {

	private StepwiseFunction<T> cumulativeFrequency;
	
	
	/**
	 * @param v: array of objects which can be result of random draw
	 * @param f: same length as v, contains relative frequencies/probabilities for drawing the objects
	 * @param r: A random number generator.
	 */
	public DiscreteDistribution(T[] values, double[] freqs) { 
		
		if (values.length != freqs.length) throw new IllegalArgumentException("Values and frequencies must have same length");
		double tot = 0;
		for (double f : freqs) tot += f;
		double[] cumFreq = new double[freqs.length-1];
		double cum = 0;
		for (int i=0; i<freqs.length-1; i++) {
			cum += freqs[i];
			cumFreq[i] = cum/tot;
		}
		cumulativeFrequency = new StepwiseFunction<T>(cumFreq, values);
	}
	
	/*
	 * TO DO: Figure out how to make constructor which takes a Collection. Following does not work
	public DiscreteDistribution(Collection<? extends T> collection, double[] freqs, MersenneTwisterFast r) {
		int n = collection.size();
		values = new T[n];     // This is where it fails
		collection.toArray(values);
		... etc ...
	}
	*/
	
	public T getValue(int i) {
		return cumulativeFrequency.getValue(i);
	}
	
	public int getNumberOfValues() {
		return cumulativeFrequency.numSteps()+1;
	}
	
	/**
	 * Draw a <T> at random from the distribution
	 * @return
	 */
	public T draw(MersenneTwisterFast rng) {
		double rand = rng.nextDouble();
		return cumulativeFrequency.f(rand);
	}
	
	/**
	 * Given a string of format "(<number>,<number>) " (repeated),
	 * parses into pairs, takes first number as value, second as weight.
	 * Returns a corresponding DiscreteDistribution.
	 * 
	 * NOTE: I have not taken the time to make this robust!
	 * 
	 * (Not sure how to generalize this to type <T>, given the need to
	 * parse. Possibly could do it so long as T is Serializable.)
	 * @param string
	 * @return
	 */
	public static DiscreteDistribution<Double> valueOf(String string) {
		final String pairDelim = "\\).*?\\(";
		String[] pairs = (")"+string+"(").split(pairDelim);
		int n=pairs.length-1;
		Double[] values = new Double[n];
		double[] freqs = new double[n];
		for (int i=0; i<n; i++) {
			final String numDelim = ",";
			String[] numbers = pairs[i+1].split(numDelim);
			values[i] = Double.valueOf(numbers[0]);
			freqs[i] = Double.valueOf(numbers[1]);
		}
		return new DiscreteDistribution<Double>(values,freqs);
	}
	
	public String toString() {
		StringBuffer buffer = new StringBuffer();
		int n = cumulativeFrequency.numSteps()+1;
		if (n==1) {
			buffer.append("(").append(cumulativeFrequency.getValue(0).toString()).append(",1)");
		} else {
			double[] freq = new double[n];
			freq[0] = cumulativeFrequency.getStep(0);
			freq[n-1] = 1 - cumulativeFrequency.getStep(n-2);
			for (int i=1; i<n-1; i++) {
				freq[i] = cumulativeFrequency.getStep(i)-cumulativeFrequency.getStep(i-1);
			}
			for (int i=0; i<n; i++) {
				buffer.append('(').append(cumulativeFrequency.getValue(i).toString()).append(',').append(freq[i]).append(')');
			}
		}
		return buffer.toString();
	}
	
	public boolean valuesInRange(T min, T max) {
		return cumulativeFrequency.valuesInRange(min, max);
	}
	public static void test() {
		//Double[] values = new Double[]{0.25,0.5,0.75};
		//double[] freqs = new double[]{1,2,1};
		MersenneTwisterFast rng = new MersenneTwisterFast();
		//DiscreteDistribution<Double> dd = new DiscreteDistribution<Double>(values,freqs,rng);
		DiscreteDistribution<Double> dd = valueOf("(0.25,1)(0.5,2),(0.75,1)");
		
		int[] count = new int[]{0,0,0};
		for (int i=0; i<1000; i++) {
			Double result = dd.draw(rng);
			if (result == 0.25) { count[0]++; }
			else if (result == 0.5) { count[1]++; }
			else { count[2]++; }
		}
		System.out.printf("Results of 1000 draws from distribution %s:\n",dd.toString());
		for (int i=0; i<3; i++) System.out.printf("Val %4f count=%d\n",dd.getValue(i),count[i]);
	}

}
