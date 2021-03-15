package ATACanalysis;

/*
 * 3-11-2019: Translating ATAC-chip analysis for Surajit Dhara from Matlab.
 * 
 * Copyright (c) 2019 by Bo Curry. All rights reserved.
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;


public class ATACanalysis {
	static final String version = "0.0.2";
	static final String[] nameFields = {"FeatureNum","ControlType","ProbeName","SystematicName"};
	static final String[] valueFields = {"gBGSubSignal","rBGSubSignal","LogRatio","gIsFeatNonUnifOL","rIsFeatNonUnifOL"};
	
	private String designdir;		// Full path to the array design directory
	private String exptdir;			// Full path to the working directory
	private String[] arrayNames;	// Simplified filenames of the array files
	private String[] sampleNames;	// From the HybSummary, if present, otherwise the arrayNames
	private int amadid;				// The design ID, extracted from the FeatureExtractor_Barcode header field.
	
	/* The actual data fields */
	private	String[] ctrltype, pnames, gnames;
	@SuppressWarnings("unused")
	private float[][] gsigs, rsigs; 	// [narrays][nfeatures]
	private float[][] logratios;	// [narrays][nfeatures]
	@SuppressWarnings("unused")
	private float[] bkg, bkgsd;
	
	// Combined data.
	private int narrays = 0;
	private int nfeatures, nprobes;
	private int[] probeno;
	
	private String[] genenames;
	private String[] probenames;
	private int[] controlType;
	private float[][] probeLogRatios;	// [iar][ip]
    private int[] ATACctrltype;			// [nprobes]-1= blue, 0=ctrl, 1= red
	
	// Probe mappings to ATAC regions
    private int[] probeRegion; 		// 1-based regions: [-1 0 N 1000+N] = [control, CGH, ATACctrl, ATACdiff]
	private String[] ATACctrlchr;
	private int[] ATACctrlstart;
	private int[] ATACctrlstop;
	private String[] ATACdiffchr;
	private int[] ATACdiffstart;
	private int[] ATACdiffstop;
	private String[] ATACdiffname;		// Name of the ATAC differential region
	private boolean[] ATACdiffisblue;	// is it blue?
    private int nredregions = 0, nblueregions = 0;
    
    // Log ratios by region
    private float[][] ATACdiffLR, ATACdiffSD, ATACctrlLR, ATACctrlSD;	// [narrays][nregions]
    
    
	/**
	 * The constructor just defines the working directory.
	 */
	public ATACanalysis(String exptdir)
	{
		this.exptdir= exptdir;
	}
	
	/* Reads in the raw data from the FE files. */
	public int readFEdata() throws IOException {
		if(exptdir == null) return(0);
		
		// Read the data
		FEdata mydata= new FEdata(exptdir, nameFields, valueFields);
		narrays= mydata.ReadFEdata();
		amadid= mydata.amadid;
		arrayNames= mydata.arrayNames;
		sampleNames= Arrays.copyOf(arrayNames, arrayNames.length);	// default
		System.out.format("Read data from %d FE files, design %d, in %s\n", narrays, amadid, exptdir);
				
		ctrltype= mydata.FEnameFieldValue(1);
		pnames= mydata.FEnameFieldValue("ProbeName");
		gnames= mydata.FEnameFieldValue("SystematicName");
		
		// Replace nonuniformity outliers with NaN
		mydata.filterValues("gBGSubSignal", "gIsFeatNonUnifOL");
		mydata.filterValues("rBGSubSignal", "rIsFeatNonUnifOL");
		// Log ratios get filtered if either channel is an outlier
		mydata.filterValues("LogRatio", "rIsFeatNonUnifOL");
		mydata.filterValues("LogRatio", "gIsFeatNonUnifOL");
		
		logratios= mydata.FEvalueFieldValue("LogRatio");
	
		bkg= mydata.backgroundLevels("gBGSubSignal", "rBGSubSignal");
		bkgsd= mydata.backgroundNoise("gBGSubSignal", "rBGSubSignal");
		//System.out.format("Backgrounds:");
		//for(int iar=0;iar<bkg.length;iar++) System.out.format("\t%.4f", bkg[iar]); 
		//System.out.format("\n");
		System.out.format("Background noise:\nGreen");
		for(int iar=0;iar<bkgsd.length/2;iar++) System.out.format(" %.4f", bkgsd[iar]); 
		System.out.format("\nRed  ");
		for(int iar=bkgsd.length/2;iar<bkgsd.length;iar++) System.out.format(" %.4f", bkgsd[iar]); 
		System.out.format("\n");
		
		return(narrays);
} // 	public int readFEdata() {

	/* Reads sample information from the "HybSum" file, if it exists in the 
	 * experimental directory. If this file exists, it's a tab delimited text file
	 * with columns headed "Array\tSample\t...\n". Here we care only about the association between
	 * array names and sample names. The array names must match those returned by FEdata.ReadFEdata(),
	 * namely "S10006_1_2", an S followed by the last 5 digits of the slide barcode followed by the
	 * subarray designation, if any. The Sample is a user defined name. Columns other than thee two
	 * are ignored.
	*/
	public void readHybSummary() throws IOException {
		if(exptdir == null) return;
		File file = new File(exptdir);
		String[] txtfiles = file.list();
        for(int i=0;i<txtfiles.length;i++) {
        	if(txtfiles[i].indexOf(".txt") < 0 || txtfiles[i].toLowerCase().indexOf("hybsum") < 0) continue;
        	// Open the file to see if it has the right format and amadid
        	// Open the txt file, and read lines 2 and 3 to get the amadid
    		BufferedReader in = new BufferedReader(new FileReader(exptdir+'\\'+txtfiles[i]));
    		String tline = in.readLine();	// header names
			String[] f = tline.split("\t");
			int arrayfield = -1, samplefield = -1;
        	for(int ifld=0;ifld<f.length;ifld++) {
        		if(f[ifld].equals("Array")) arrayfield= ifld;
        		if(f[ifld].equals("Sample")) samplefield= ifld;
        	}
        	if(arrayfield < 0 || samplefield < 0) {
        		System.out.format("%s is not a valid HybSummary file\n");	// txtfiles[i]
        		continue;	// but maybe another one is
        	}
        	// Read each line of the HybSummary file, an array name and a sample name.
       		while(true) {
       	   		tline = in.readLine();
    	        if(tline == null) break;
    	        else if(tline.isEmpty()) continue;
    	        f = tline.split("\t");
    	        String aname = f[arrayfield];
    	        // Find an exact match with one of the expt arrays. This is inefficient, but so what.
    	        for(int iar=0; iar<arrayNames.length;iar++) {
    	        	if(aname.equals(arrayNames[iar])) {
    	        		sampleNames[iar]= f[samplefield];
    	        		break;
    	        	}
    	        }
       		}
            in.close();
        }
	}	// readHybSummary
	
	
	/* Combines replicate features. */
	public void combineFeaturesToProbes() {
		
		if(narrays == 0) return;
		nfeatures = gnames.length;
		probeno= new int[nfeatures];
		
		/* 
		 * We need a sort that, like Matlab's, reports where the sorted values came from, so we can unsort.
		 */
		StringSorter sortnames = new StringSorter(gnames);
		int[] indices = sortnames.sort();	// We have to be able to unsort, too
		String[] sgnames= sortnames.sortedArray();
		
		ArrayList<String> genenamelist = new ArrayList<String>();
		genenamelist.add(sgnames[0]);
		probeno[indices[0]]= 0;
		for(int ifeat= 1, ip= 0;ifeat<nfeatures;ifeat++) {
			if(!sgnames[ifeat].equals(sgnames[ifeat-1])) {
				genenamelist.add(sgnames[ifeat]);
				ip++;
			}
			probeno[indices[ifeat]]= ip;
		}
		genenames= genenamelist.toArray(new String[0]);
		nprobes = genenames.length;

		//System.out.format("Found %d unique genes [%s, %s, ...] among %d features\n", nprobes, genenames[0], genenames[1], nfeatures);

		probenames = new String[nprobes];
		controlType = new int[nprobes];
		int[] nprobereps = new int[nprobes];	// This can stay local
		int maxreps = 0;
		for(int ifeat= 0;ifeat<nfeatures;ifeat++) {
			int ip= probeno[ifeat];
			probenames[ip]= pnames[ifeat];
			controlType[ip]= (int)Integer.valueOf(ctrltype[ifeat]);
			nprobereps[ip]++;
			if(controlType[ip] == 0 && nprobereps[ip] > maxreps) {
				maxreps= nprobereps[ip];
			}
		}
		System.out.format("Found %d unique genes with a maximum of %d reps\n", nprobes, maxreps);
		
		probeLogRatios = new float[narrays][nprobes];
		for(int iar=0;iar<narrays;iar++) {
			float[][] lrreps= new float[nprobes][maxreps];
			int[] nrepsfound= new int[nprobes];
			//if(iar==0) System.out.format("Uninitialized lrs= [%.2f,%.2f,...]\n", lrreps[0][0], lrreps[0][1]);
			for(int ifeat= 0;ifeat<nfeatures;ifeat++) {
				int ip= probeno[ifeat];
				if(controlType[ip] != 0) continue;	// Don't compute median  LR for control probes
				//if(nrepsfound[ip] >= maxreps)
				//	System.out.format("Oops, we found %d replicates of probe %s\n", nrepsfound[ip], genenames[ip]);
				float lr = logratios[iar][ifeat]/(float)Math.log10(2.);
				if(!Float.isNaN(lr)) lrreps[ip][nrepsfound[ip]++]= lr;
			}
			for(int ip=0;ip<nprobes;ip++) {
				// The probe LR is the median of the non-nan replicate lrs
				probeLogRatios[iar][ip]= median(lrreps[ip],nrepsfound[ip]);
			}
			if(iar<0) {
				// for debugging
				int ip = 6768;	// The first probe in ATACdiff[0]
				System.out.format("Lrs[0][%d] (%s)=", ip, genenames[ip]);
				for(int irep=0;irep<nrepsfound[ip];irep++) System.out.format(" %.4f", lrreps[ip][irep]);
				System.out.format("\nMedian= %.4f\n", probeLogRatios[iar][ip]);
			}
			
		}
		
		// *** Return or write out some measure of replicate reproducibility?
	}

	/* Read the mapping from genenames to significant regions from the design directory.
	 *  Probes are either:
		 * Control grid probes
		 * CGH backbone probes
		 * ATAC control region probes
		 * ATAC blue differential probes
		 * ATAC red differential probes
		 * These are distinguished by parsing their systematic names (genenames).
		 * Returns the total number of probes found in each category.
		*/
	public int[] readProbeRegions(String designpath) throws IOException {
		if(exptdir == null) return(null);
		designdir= designpath;
		final String ATACctrlfile= "ATACctrlRegions.txt";
		final String ATACdifffile= "ATACdiffRegions.txt";
		
		// Control file format: Chromosome	start	end	length (bp)
		BufferedReader in = new BufferedReader(new FileReader(designdir+'\\'+ATACctrlfile));
		int nlines = -1;	// it overcounts by 1
		// Read through once to get a line count
		String tline = in.readLine();	// header
		while(tline != null) {
			tline = in.readLine();
			nlines++;
		}
		ATACctrlchr = new String[nlines];
		ATACctrlstart = new int[nlines]; 
		ATACctrlstop = new int[nlines];
		in.close();
		//System.out.format("Found %d regions in %s\n", nlines, ATACctrlfile);
		in = new BufferedReader(new FileReader(designdir+'\\'+ATACctrlfile));
		tline = in.readLine();	// header
		for(int iline=0;iline<nlines;iline++) {
			tline = in.readLine();
			//System.out.format("Line[%d]: %s\n", iline, tline);
			String[] f = tline.split("\t");
			ATACctrlchr[iline]= f[0];
			ATACctrlstart[iline]= Integer.valueOf(f[1]);
			ATACctrlstop[iline]= Integer.valueOf(f[2]);
		}
		in.close();
		
		// Differential file format: seqnames	start	end	annot	transcript_id	symbol	bp	Isblue?
		in = new BufferedReader(new FileReader(designdir+'\\'+ATACdifffile));
		nlines = -1;	// it overcounts by 1
		// Read through once to get a line count
		tline = in.readLine();	// header
		while(tline != null) {
			tline = in.readLine();
			nlines++;
		}
		ATACdiffchr = new String[nlines];
		ATACdiffstart = new int[nlines]; 
		ATACdiffstop = new int[nlines];
		ATACdiffname = new String[nlines];
		ATACdiffisblue = new boolean[nlines];
		in.close();
		in = new BufferedReader(new FileReader(designdir+'\\'+ATACdifffile));
		tline = in.readLine();	// header
		for(int iline=0;iline<nlines;iline++) {
			tline = in.readLine();
			String[] f = tline.split("\t");
			ATACdiffchr[iline]= f[0];
			ATACdiffstart[iline]= Integer.valueOf(f[1]);
			ATACdiffstop[iline]= Integer.valueOf(f[2]);
			ATACdiffname[iline]= f[5];
			ATACdiffisblue[iline]= (f.length == 8 && f[7].equals("y"));
		}
		in.close();
		
	    // Map array probes to regions. Any non-control probe that doesn't map to an ATAC
	    // region is considered to be a CGH backbone probe.
		// These are encoded <0, 0, >0<1000 >1000 for Control, CGH, ATAC control, ATAC differential
		probeRegion = new int[genenames.length];	// 1-based
		int[] regioncount = {0, 0, 0, 0};	// [Control, CGH, ATAC ctrl, ATAC diff]
		for(int ip=0;ip<genenames.length;ip++) {
			if(controlType[ip] != 0) {
				probeRegion[ip]= -1;
				regioncount[0]++;
				continue;
			}
			int imcolon = genenames[ip].indexOf(':');
			int imdash = genenames[ip].indexOf('-');
			String chrstr= genenames[ip].substring(0, imcolon);
			int start = Integer.valueOf(genenames[ip].substring(imcolon+1,imdash));
			int end = Integer.valueOf(genenames[ip].substring(imdash+1));
	        // See if this probe is in an ATAC control region. If this is too slow, we
	        // can sort them to speed it up.
	        Integer[] onchr = strmatch(chrstr, ATACctrlchr, true);
	        for(int ireg=0;ireg<onchr.length;ireg++) {
		        // it counts if the probe overlaps the region at all
	        	if(ATACctrlstart[onchr[ireg]]<=end && ATACctrlstop[onchr[ireg]]>=start) {
	        			probeRegion[ip]= onchr[ireg] + 1;
	        			break;
	        	}
	        }
	        if(probeRegion[ip]!=0) {
				regioncount[2]++;
	        	continue;
	        }
	        // If this probe isn't in an ATAC control region, maybe it's in a differential region.
	        onchr = strmatch(chrstr, ATACdiffchr, true);
	        for(int ireg=0;ireg<onchr.length;ireg++) {
		        // it counts if the probe overlaps the region at all
	        	if(ATACdiffstart[onchr[ireg]]<=end && ATACdiffstop[onchr[ireg]]>=start) {
	        			probeRegion[ip]= 1001 + onchr[ireg];
	        			break;
	        	}
	        }
	        if(probeRegion[ip]!=0) regioncount[3]++;
	        else regioncount[1]++; // By default it's a CGH backbone probe
		}	// 		for(int ip=0;ip<genenames.length;ip++) {

		
	    // Annotate ATAC probes by ATAC control type: [-1, 0, 1]= [blue, ctrl, red]
		// There are multiple probes per region.
		ATACctrltype= new int[regioncount[2]+regioncount[3]];
		for(int ip=0, natac=0;ip<nprobes;ip++) {
			if(probeRegion[ip]<=0) continue;
			if(probeRegion[ip] < 1000) ATACctrltype[natac++]= 0;
			else if(ATACdiffisblue[probeRegion[ip]-1001]) ATACctrltype[natac++]= -1;
			else ATACctrltype[natac++]= 1;
		}

		return(regioncount);		
	} // public int[] readProbeRegions()
	
	/* Requires the prior identification of the CGH background probes (probeRegion==0). */
	public void centerLogratios() {
		if(probeRegion == null) return;
		float[] CGHlrs = new float[nprobes];
		float[] medLRs = new float[narrays];
		for(int iar=0;iar<narrays;iar++) {
			int nCGH = 0;
			for(int ip=0;ip<nprobes;ip++) {
				if(probeRegion[ip] == 0 && !Float.isNaN(probeLogRatios[iar][ip])) {
					CGHlrs[nCGH++]= probeLogRatios[iar][ip];
				}
			}
			medLRs[iar]= median(CGHlrs, nCGH);
			//System.out.format("Found %d CGH probes on array %d, averaging %.4f (median %.4f)\n", nCGH, iar, mean(CGHlrs), medLRs[iar]);
			for(int ip=0;ip<nprobes;ip++) probeLogRatios[iar][ip] -= medLRs[iar];
		}
			System.out.format("Median LRs of CGH probes:");
			for(int iar=0;iar<narrays;iar++) {
				System.out.format(" %.4f",  medLRs[iar]);
			}
			System.out.format("\n");
	}
	
	/* Must be called *after* the LRs have been centered
	 */
	public void combineProbesInRegions() {
		
		// Compute median LR and LR spread for ATAC regions. This is n^2.
	    ATACdiffLR = new float[narrays][ATACdiffchr.length];
	    ATACdiffSD = new float[narrays][ATACdiffchr.length];
	    ATACctrlLR = new float[narrays][ATACctrlchr.length];
	    ATACctrlSD = new float[narrays][ATACctrlchr.length];
	    
	    for(int ireg=0;ireg<ATACdiffchr.length;ireg++) {
	    	// Find all the probes in this region.
	    	int[] imin = new int[50];	// 50 should be enough
	    	int nlrs = 0;
	    	for(int ip=0;ip<nprobes;ip++)
	    		if(probeRegion[ip] == ireg+1001) imin[nlrs++]= ip;
	    	// Now average the LRs in this region for each array
	    	float[] lrs = new float[nlrs];
			for(int iar=0;iar<narrays;iar++) {
				for(int ip=0;ip<nlrs;ip++) lrs[ip]= probeLogRatios[iar][imin[ip]];
		    	ATACdiffLR[iar][ireg]= median(lrs);
		    	ATACdiffSD[iar][ireg]= iqr(lrs)/1.349F;	// /2/sqrt(2)/erfinv(.5)
			}
			if(ireg<0) {
				// For debugging
				System.out.format("Lrs[0][%d] (%s:%d-%d)=", ireg, ATACdiffchr[ireg], ATACdiffstart[ireg], ATACdiffstop[ireg]);
				for(int ip=0;ip<nlrs;ip++) System.out.format(" %.4f", probeLogRatios[0][imin[ip]]);
				System.out.format("\n");
			}
	    }

	    for(int ireg=0;ireg<ATACctrlchr.length;ireg++) {
	    	// Find all the probes in this region.
	    	int[] imin = new int[50];	// 50 should be enough
	    	int nlrs = 0;
	    	for(int ip=0;ip<nprobes;ip++)
	    		if(probeRegion[ip] == ireg+1) imin[nlrs++]= ip;
	    	// Now average the LRs in this region for each array
	    	float[] lrs = new float[nlrs];
			for(int iar=0;iar<narrays;iar++) {
				for(int ip=0;ip<nlrs;ip++) lrs[ip]= probeLogRatios[iar][imin[ip]];
				ATACctrlLR[iar][ireg]= median(lrs);
				ATACctrlSD[iar][ireg]= iqr(lrs)/1.349F;	// /2/sqrt(2)/erfinv(.5)
			}
	    }
	} // 	public void combineProbesInRegions() {


	public void writeRegions(String outfile) throws IOException {
		// *** PrintStream is *much* slower than BufferedWriter, but required for format()
		BufferedWriter outstream = new BufferedWriter(new FileWriter(exptdir+"\\"+outfile));
	    outstream.write("chr\tstart\tend\tsymbol\ttype");
	    for(int iar=0;iar<narrays;iar++) outstream.write("\t"+sampleNames[iar]+" L2R");
	    outstream.write("\n");
	    // Write out ATAC diff probes first, ATAC ctrl probes next, CGH probes last.
	    for(int ireg=0;ireg<ATACdiffchr.length;ireg++) {
	    	outstream.write(ATACdiffchr[ireg]+"\t"+ATACdiffstart[ireg]+"\t"+ATACdiffstop[ireg]+"\t");
	    	outstream.write(ATACdiffname[ireg]+"\t");
	    	outstream.write((ATACdiffisblue[ireg]) ? "blue" : "red");
	    	for(int iar=0;iar<narrays;iar++) outstream.write(String.format("\t%.4f", ATACdiffLR[iar][ireg]));
		    outstream.write("\n");
	    }
	    // Now the ATAC ctrl probes.
	    for(int ireg=0;ireg<ATACctrlchr.length;ireg++) {
	    	outstream.write(ATACctrlchr[ireg]+"\t"+ATACctrlstart[ireg]+"\t"+ATACctrlstop[ireg]+"\t\tctrl");
	    	for(int iar=0;iar<narrays;iar++) outstream.write(String.format("\t%.4f", ATACctrlLR[iar][ireg]));
		    outstream.write("\n");
	    }
	    // Finally, the CGH backbone probes
	    for(int ip=0;ip<nprobes;ip++) {
			if(probeRegion[ip] != 0) continue;
			String[] f = genenames[ip].split("[:-]");
	    	outstream.write(f[0]+"\t"+f[1]+"\t"+f[2]+"\t\tCGH");
	    	for(int iar=0;iar<narrays;iar++) outstream.write(String.format("\t%.4f", probeLogRatios[iar][ip]));
		    outstream.write("\n");
	    }
	    outstream.close();
	}

	/* Returns a boolean vector of the redness of each sample. */
	public boolean[] isSampleRed() {
		boolean[] isred = new boolean[sampleNames.length];
		if(nredregions == 0)
			for(int ireg=0;ireg<ATACdiffchr.length;ireg++) {
				if(ATACdiffisblue[ireg]) nblueregions++;
				else nredregions++;
			}
		float[] bluelrs = new float[nblueregions];
		float[] redlrs = new float[nredregions];
		System.out.format("Scoring %d blue and %d red regions\n", nblueregions, nredregions);
		for(int isamp=0;isamp<sampleNames.length;isamp++) {
			@SuppressWarnings("unused")
			int nblue = 0, nred = 0, nnan = 0;
			for(int ireg=0;ireg<ATACdiffchr.length;ireg++) {
				if(ATACdiffisblue[ireg] && !Float.isNaN(ATACdiffLR[isamp][ireg]))
					bluelrs[nblue++]= ATACdiffLR[isamp][ireg];
				else if(!Float.isNaN(ATACdiffLR[isamp][ireg]))
					redlrs[nred++]= ATACdiffLR[isamp][ireg];
				if(Float.isNaN(ATACdiffLR[isamp][ireg])) nnan++;
			}
			float medred = median(redlrs, nred);
			float medblue = median(bluelrs, nblue);
			//System.out.format("Sample %d averages %.4f in %d red and %.4f in %d blue regions (%d NaNs)\n", 
			//				isamp, medred, nred, medblue, nblue, nnan);
			isred[isamp]= (medred > medblue); 
		}
		return(isred);
	}	// public boolean[] isSampleRed() {
	
	/* Does a ttest on the red and blue regions, reports the mlogp of their being equal. */
	public float[] diffMlogp() {
		float[] mlogp = new float[sampleNames.length];
		if(nredregions == 0)
			for(int ireg=0;ireg<ATACdiffchr.length;ireg++) {
				if(ATACdiffisblue[ireg]) nblueregions++;
				else nredregions++;
			}
		float[] bluelrs = new float[nblueregions];
		float[] redlrs = new float[nredregions];
		for(int isamp=0;isamp<sampleNames.length;isamp++) {
			int nblue = 0, nred = 0;
			for(int ireg=0;ireg<ATACdiffchr.length;ireg++) {
				if(ATACdiffisblue[ireg] && !Float.isNaN(ATACdiffLR[isamp][ireg]))
					bluelrs[nblue++]= ATACdiffLR[isamp][ireg];
				else if(!Float.isNaN(ATACdiffLR[isamp][ireg]))
					redlrs[nred++]= ATACdiffLR[isamp][ireg];
			}
			double pval = ttest2(Arrays.copyOf(redlrs,nred), Arrays.copyOf(bluelrs,nblue));
			mlogp[isamp]= (float)-Math.log10(pval); 
		}
		return(mlogp);
	}

	//----------------------------- Utility methods ------------------------------------------

	/* One of the most useful Matlab functions, which Mathworks keeps threatening to take away. */
	static Integer[] strmatch(String target, String[] library, boolean isexact) {
		ArrayList<Integer> hitlist = new ArrayList<Integer>();
		for(int i=0;i<library.length;i++)
			if((!isexact && library[i].startsWith(target)) || (isexact && library[i].equals(target)))
				hitlist.add(i);
		Integer[] hits;
		hits= hitlist.toArray(new Integer[0]);
		return(hits);
	} // 	static int[] strmatch(String target, String[] library, boolean isexact) {

	static float median(float[] data) {
		return(median(data, data.length));
	}
	
	// Ignoring NaNs. Matlab returns NaN if any element is NaN, which seems wrong.
	// *** It gives the same results - why?
	static float median(float[] data, int length) {
		if(length == 0) return(Float.NaN);
		float[] sdata = Arrays.copyOf(data, length);
		Arrays.sort(sdata);	// NaNs sort high
		for(int i=sdata.length-1;i>=0&&Float.isNaN(sdata[i]);i--,length--);	// ignore NaNs at the end
		float medval;
		if(length <= 0) medval= Float.NaN;
		else if((length%2)==0) medval= 0.5F*(sdata[length/2-1]+sdata[length/2]);
		else			  medval= sdata[(length-1)/2];
		if(sdata.length == 0 && length < sdata.length) {
			// For debugging only
			int len= sdata.length;
			System.out.format("median() found %d NaNs of %d, reports %.4f not %.4f\n", len - length, len,
					medval, ((len%2)==0) ? 0.5F*(sdata[len/2-1]+sdata[len/2]) : sdata[(len-1)/2]);
		}
		return(medval);
	}

	static float iqr(float[] data) {
		return(iqr(data, data.length));
	}
	
	// *** Not quite right, but close.Ignore NaNs.
	static float iqr(float[] data, int length) {
		if(length == 0) return(Float.NaN);
		if(length == 1) return(0);
		float[] sdata = Arrays.copyOf(data, length);
		Arrays.sort(sdata);
		for(int i=sdata.length-1;i>=0&&Float.isNaN(sdata[i]);i--,length--);	// ignore NaNs, which sort at the end
		if(length <= 0) return(Float.NaN);
		if(length == 1) return(0);
		int q1= (int)((length+2)/4)-1;	// rounding
		int q3= (int)((3*length+2)/4)-1;	// rounding
		return(sdata[q3]-sdata[q1]);
	}

	static float mean(float[] data) {
		return(mean(data, data.length));
	}
	
	static float mean(float[] data, int length) {
		float sum = 0.F;
		for(int i=0; i<length;i++) sum += data[i];
		return(sum/length);
	}

	static float variance(float[] data) {
		return(variance(data, data.length));
	}
	
	static float variance(float[] data, int length) {
		float avg = 0.F, var = 0.F;
		for(int i=0; i<length;i++) avg += data[i];
		avg /= length;
		for(int i=0; i<length;i++)	var += (data[i]-avg)*(data[i]-avg);
		return(var/(length-1));
	}

	// From "Numerical Recipes in C"
	static double ttest2(float[] data1, float[] data2) {
		float avg1 = mean(data1);
		float avg2 = mean(data2);
		float var1 = variance(data1);
		float var2 = variance(data2);
		int df = data1.length + data2.length - 2;	// degrees of freedom
		float svar = ((data1.length-1)*var1+(data2.length-1)*var2)/df;
		double tstat = (double)(avg1-avg2)/Math.sqrt((double)svar*(1.0/data1.length+1.0/data2.length));
		double pval = betainc(df/(df+tstat*tstat),0.5*df,0.5);	// Incomplete beta function
		//System.out.format("betainc(%.4f,%.1f,%.1f) = %.4f\n", df/(df+tstat*tstat),0.5*df,0.5, pval);
		return(pval);
	}
	// Incomplete beta function, from "Numerical Recipes in C". betainc(x,z,w) in Matlab.
	static double betainc(double x, double a, double b) {
		double bt;
		if(x < 0.0 || x > 1.0) return(Double.NaN);	// Illegal parameters
		if(x == 0.0 || x == 1.0) return(x);
		//System.out.format("gammaln(%.1f) = %.4f\n", a+b, gammaln(a+b));
		// exp(gammaln(a+b)-gammaln(a)-gammln(b)) == 1/Beta(a,b)
		bt= Math.exp(gammaln(a+b)-gammaln(a)-gammaln(b)+a*Math.log(x)+b*Math.log(1.0-x));
		//System.out.format("beta(%.1f,%.1f) = %.4f\n", a, b, Math.exp(gammaln(a)+gammaln(b)-gammaln(a+b)));
		// so bt is x^a*(1-x)^b/Beta(a,b)
		if(x < (a+1.0)/(a+b+2.0)) return(bt*betacf(x,a,b)/a); // Use continued fraction directly
		else return(1.0 - bt*betacf(1.0-x,b,a)/b); 			// Use continued fraction after symmetry transformation
	}

	/* The log of the gamma function (gammaln in Matlab). Pg 168 in "Numerical Recipes in C".
	 * If N is a positive integer, gamma(n) == (n-1)! Since this can get big fast, we compute th e log. */
	static double gammaln(double xx) {
		double x, tmp, ser;
		final double[] cof = {76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5};
		
		x= xx - 1;
		tmp= x + 5.5;
		tmp -= (x+0.5)*Math.log(tmp);
		ser= 1.0;
		for(int j=0;j<cof.length;j++) {
			x += 1;
			ser += cof[j]/x;
		}
		return(-tmp+Math.log(2.50662827465*ser));
	}
	
	/* A quickly converging continued fraction computation of the iterated part of betainc,
	 *  from "Numerical Recipes in C".
	 *  Computes: 1/(1+d1/(1+d2/(1+d3/(1+... , where
	 *  d(2m+1)= (a+m)(a+b+m)x/(a+2m)/(a+2m+1)
	 *  d(2m)= m(b-m)x/(a+2m-1)/(a+2m)
	 *  */
	static double betacf(double x, double a, double b) {
		final int ITMAX = 100;			// max iterations
		final double EPSILON = 3e-7;	// stopping criterion
		double em, d;
		double bz, bm=1.0, bp, bpp;
		double az=1.0, am= 1.0, ap, app, aold;
		
		bz= 1 - (a+b)*x/(a+1);
		for(int m=1;m<=ITMAX;m++) {
			em= (double)m;
			d= em*(b-em)*x/(a+2*em-1)/(a+2*em);			// Even m (2,4,6,...)
			ap= az + d*am;
			bp= bz + d*bm;
			d= -(a+em)*(a+b+em)*x/(a+2*em)/(a+2*em+1);	// Odd m (3,5,7,...)
			app= ap + d*az;
			bpp= bp + d*bz;
			aold= az;
			am= ap/bpp;
			bm= bp/bpp;
			az= app/bpp;
			bz= 1.0;
			if(Math.abs(az-aold) < EPSILON*Math.abs(az)) {
				//System.out.format("betacf(%.4f,%.1f,%.1f) converged to %.4f in %d\n", x, a, b, az, m);
				return(az);
			}
		}
		return(Double.NaN);
	}
	
	//----------------------------- Main() and utility classes ------------------------------------------
	
	/* Called from the command line. 
	 * Arguments: [0] The pathname to the directory containing the probe design files
	 * 			  [1] The pathname to the directory containing the FE data files.
	 *  */
	public static void main(String[] args) throws IOException {
		if(args.length < 2) return;
		
		System.out.format("Running ATACanalysis ver. %s on data in %s\n", version, args[1]);
		System.out.format("   using ATAC regions in %s\n", args[0]);
		
		ATACanalysis thisanalysis= new ATACanalysis(args[1]);
		
		int nsamples = thisanalysis.readFEdata();
		if(nsamples==0)	System.out.format("Can't find any FE files in %s\n", args[1]);
		
		// Read the HybSummary.txt file, if there is one.
		thisanalysis.readHybSummary();
			
		thisanalysis.combineFeaturesToProbes();
		
		//for(int isamp=0;isamp<nsamples;isamp++)
		//System.out.format("%s\n", thisanalysis.sampleNames[isamp]);
		
		/* Now read the probe classifications. */
		int[] probesinRegions = thisanalysis.readProbeRegions(args[0]);
		System.out.format("Found [%d, %d, %d, %d] probes in [control, CGH, ATAC ctrl, ATAC diff] regions\n", 
				probesinRegions[0], probesinRegions[1], probesinRegions[2], probesinRegions[3]);
		
		// Center the CGH backbone probes
		thisanalysis.centerLogratios();
		
		/* Score the log ratios of the regions. */
		thisanalysis.combineProbesInRegions();
		
		// Save out this "just data" file. Order as ATAC diff, ATAC ctrl, CGH
		thisanalysis.writeRegions("ATACregionLRs.txt");

		// Now look for differential expression in the red and blue regions
		boolean[] isred = thisanalysis.isSampleRed();
		System.out.format("Red samples are");
		for(int isamp=0;isamp<nsamples;isamp++) if(isred[isamp]) 
			System.out.format(" %s", thisanalysis.sampleNames[isamp]);
		System.out.format("\n");
		
		// How significant is the separation?
		float[] mlogp = thisanalysis.diffMlogp();
		System.out.format("Differential significance (mlogp)");
		for(int isamp=0;isamp<nsamples;isamp++) System.out.format(" %.2f", mlogp[isamp]);
		System.out.format("\n");
		
		// Debugging gammaln, a la Matlab gammaln. It works!
		//for(int i=1;i<20;i+=2) System.out.format("gammaln(%d)= %.4f\n", i, gammaln((double)i));
		
		// Debugging betainc, a la Matlab betainc(.5,(0:10)',3). It works!
		//for(int i=0;i<=10;i++) System.out.format("betainc(0.5,%d,3)= %.12f\n", i, betainc(0.5,(double)i,3.0));
		
		// Write out the classification results
		PrintStream outstream = new PrintStream(thisanalysis.exptdir+"\\SampleClassification.txt");
	    outstream.format("Sample\tclass\tmlogp\n");
	    for(int isamp=0;isamp<nsamples;isamp++) {
	    	if(mlogp[isamp] < 1.0) outstream.format("%s\tno call\t%.2f\n", thisanalysis.sampleNames[isamp], mlogp[isamp]);
	    	else if(mlogp[isamp] < 2.0 && isred[isamp]) outstream.format("%s\tleans red\t%.2f\n", thisanalysis.sampleNames[isamp], mlogp[isamp]);
	    	else if(mlogp[isamp] < 2.0 && !isred[isamp]) outstream.format("%s\tleans blue\t%.2f\n", thisanalysis.sampleNames[isamp], mlogp[isamp]);
	    	else if(isred[isamp]) outstream.format("%s\tred\t%.2f\n", thisanalysis.sampleNames[isamp], mlogp[isamp]);
	    	else outstream.format("%s\tblue\t%.2f\n", thisanalysis.sampleNames[isamp], mlogp[isamp]);
	    }
	    outstream.close();
	}
	
	/* 
	 * We need a sort that, like Matlab's, reports where the sorted values came from, so we can unsort.
		I use this class so I can unsort a sorted array of strings. Could be generalized.
	 */
	class StringSorter {

		@SuppressWarnings("unused")
		private FirstIntCompare icompare = new FirstIntCompare();
		private FirstStringCompare scompare = new FirstStringCompare();
		private String[][] indexedArray;
				
		/* The constructor makes a copy of te list to be sorted, matched with aindex strings.
		 * *** Won't work for length > 9999999 *** */
		StringSorter(String[] strarray) {
			this.indexedArray= new String[strarray.length][2];
			for(int i=0;i<strarray.length;i++) {
				this.indexedArray[i][0]= strarray[i];
				this.indexedArray[i][1]= String.format("%07d", i);
			}
			//System.out.format("The first element of indexedArray is [%s, %s]\n", this.indexedArray[0][0], this.indexedArray[0][1]);
		}
		
		public int[] sort() {
			//System.out.format("Sorting indexedArray[%d][%d] begins [%s, %s]\n", 
			//		indexedArray.length, indexedArray[0].length, indexedArray[0][0], indexedArray[0][1]);
			Arrays.sort(indexedArray, scompare);
			// Return the original indices of the sorted array
			int[] indices = new int[indexedArray.length];
			for(int i=0;i<indexedArray.length;i++) indices[i]= Integer.valueOf(indexedArray[i][1]);
			return(indices);
		}

		public String[] sortedArray() {
			String [] sorted = new String[indexedArray.length];
			for(int i=0;i<indexedArray.length;i++) sorted[i]= indexedArray[i][0];
			return(sorted);
		}
		
		// Sort based on the first integr of the elements of the array
		class FirstIntCompare implements Comparator<int[]> {
			// "Note: this comparator imposes orderings that are inconsistent with equals.
			public int compare(int[] a, int[] b) {
				return a[0] - b[0];
			}
		}
		// Sort based on the first string of the array elements 
		class FirstStringCompare implements Comparator<String[]> {
			// "Note: this comparator imposes orderings that are inconsistent with equals.
			public int compare(String[] a, String[] b) {
				return (a[0].compareTo(b[0]));
			}
		}
	}
	
}
