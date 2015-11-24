package hybridseq;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.bio.phylo.io.nexus.CharactersBlock;
import org.biojavax.bio.phylo.io.nexus.DataBlock;
import org.biojavax.bio.phylo.io.nexus.DistancesBlock;
import org.biojavax.bio.phylo.io.nexus.NexusBlockParser;
import org.biojavax.bio.phylo.io.nexus.NexusFile;
import org.biojavax.bio.phylo.io.nexus.NexusFileBuilder;
import org.biojavax.bio.phylo.io.nexus.NexusFileFormat;
import org.biojavax.bio.phylo.io.nexus.TaxaBlock;
import org.biojavax.bio.phylo.io.nexus.TreesBlock;

import biojavaExtensions.GenericBlockParser;
import biojavaExtensions.UseableUnknownBlockParser;

public class MainConfiguration {
	/**
	 * Holds all the configuration variables for Main.
	 * @author Michael Woodhams
	 *
	 */
	
	public HybridNetworkParameters params=null;
	public String inputFile = "input.nex";
	public String outputFile = "output.nex";
	public String filoOutputFile = "filo.out";
	public NexusFile nexusFile;
	
	public MainConfiguration(String[] args, PrintWriter out) {
		
		long cmdLineSeed = 0;
		String stringArg = null;
		String[][] parsed = null;
		try {
			parsed = mdwUtils.Strings.parseCmdLine(args);
		} catch (RuntimeException e) {
			usage();
			System.exit(1);
		}
		for (String[] arg : parsed) {
			char flag = arg[0].charAt(1);
			long count;
			@SuppressWarnings("unused")
			float real;
			// real/count for options which take float or int values
			// (currently float values not used by any option.)
			stringArg = arg[1];
			try {
				real = Float.valueOf(arg[1]);
			} catch (NumberFormatException e) {
				real = Float.NaN;
			}
			try {
				count = Long.valueOf(arg[1]);
			} catch (NumberFormatException e) {
				count = Long.MIN_VALUE;
			}
			
			// Possible flag values are 'real', 'count', 'stringArg'.
			
			switch (flag) {
			case 'i' :  inputFile = stringArg;
						break;
			case 'o' :  outputFile = stringArg;
						break;
			case 'f' :  filoOutputFile = stringArg;
						break;
			case 's' :  cmdLineSeed = count;
						break;
			default:	usage(); 	
						throw new IllegalArgumentException("Unrecognized flag: "+flag); 
			}
		}
		File file = new File(inputFile);
		NexusFileBuilder builder=new NexusFileBuilder();
		// The two custom Nexus blocks used by this program:
		builder.setBlockParser(HybridNetworkParameters.BLOCK_NAME, new GenericBlockParser(false,HybridNetworkParameters.HYBRIDSIM_VALID_KEYS));
		builder.setBlockParser("filo", new GenericBlockParser(true));
		/*
		 *  For Nexus block types NOT used by this program, use the simple 'unknown' parser.
		 *  This will parse pretty much anything, so we avoid issues with parse errors on
		 *  stuff that will never be used anyhow. 
		 */
		builder.setBlockParser(TaxaBlock.TAXA_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(TreesBlock.TREES_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(CharactersBlock.CHARACTERS_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(DataBlock.DATA_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(DistancesBlock.DISTANCES_BLOCK, new UseableUnknownBlockParser());
		builder.setBlockParser(NexusBlockParser.UNKNOWN_BLOCK,	new UseableUnknownBlockParser());
		
		try {
			NexusFileFormat.parseFile(builder, file);
		} catch (ParseException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		nexusFile = builder.getNexusFile();
		//params = new HybridNetworkParameters(nexusFile,out);
		params = new HybridNetworkParameters(nexusFile);
		if (cmdLineSeed!=0) params.seed = cmdLineSeed;
	}
	
	public void usage() {
		System.out.println("Command line options:");
		System.out.println("-i<file>: Input Nexus file");
		System.out.println("-o<file>: Output file");
		System.out.println("-f<file>: Filo output file");
		System.out.println("-s<seed>: Set random number generator seed (takes precedent over seed in input file)");
	}
}
