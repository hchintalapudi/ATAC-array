package ATACanalysis;

/*
 * 3-20-2019: Reads and manipulates data from Agilent Feature Extraction output files.
 * Reads data from multiple output *.txt files of a single array design, in a single directory.
 * All *.txt files in the working directory are considered as candidates for FE output files,
 * but we only keep those that (a) are parseable as FE output files, and (b) share an AMADID with
 * the first such parseable file.
 * Translated from Matlab ReadFEtxtFile.m.
 * 
 * * Copyright (c) 2019 by Bo Curry. All rights reserved.
 *
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.io.IOException;
//import java.util.ArrayList;
import java.util.Arrays;

public class FEdata {
	static final String version = "0.0.1";
	static final String[] supportedTypes= {"integer","float","text","boolean"};
	static final String[] scanfTypeStr= {"%d","%f","%s","%d"};

	private String exptdir;			// Full path to the working directory
	private String[] FEfiles;		// The filenames of the *txt files
	private String[] nameFields;	// The FE name fields (common to all files of this design) to read
	private String[] valueFields;	// The FE value fields (unique for each file) to read
	
	private String[][] featureNames;	// [namefields,arrayfeatures] stored as integers or strings
	private float[][][] featureValues;	// [valuefields,arrays,arrayfeatures] all stored as floats
	private float[] bkg = {};			// [gbkg, rbkg] median signals of negative controls
	private float[] bkgsd = {};			// [gbkg, rbkg] spread of negative control signals
	
	public int amadid;				// The design ID, extracted from the FeatureExtractor_Barcode header field.
	public String[] arrayNames;			// Simplified filenames
	
	/*
	 * The constructor just defines the working directory and the desired fields.
	 */
	public FEdata(String exptdir, String[] nameFields, String[] valueFields)
	{
		this.exptdir= exptdir;
		//System.out.format("Looking for %d name fields and %d value fields in %s\n", nameFields.length, valueFields.length, exptdir);
		this.nameFields= new String[nameFields.length];
		this.valueFields= new String[valueFields.length];
		System.arraycopy(nameFields,0,this.nameFields,0,nameFields.length);
		System.arraycopy(valueFields,0,this.valueFields,0,valueFields.length);
		this.amadid= 0;	// It's virgin
	}
	
	/*
	 * Read in the data from whatever parseable files we can find. Returns the number of files read.
	 */
	public int ReadFEdata() throws IOException
	{
		File file = new File(exptdir);
		String[] fieldtypes;
		String[] fieldnames;
		String[] fieldvals;
		FEfiles = file.list();
		int nfiles= 0, nlines= 0;
		String tline;
		int[] namecols= new int[nameFields.length];
		int[] nametypes= new int[nameFields.length];
		int[] valuecols= new int[valueFields.length];
		
        // Look only at *.txt files
        for(int i=0;i<FEfiles.length;i++) {
        	if(FEfiles[i].indexOf(".txt") < 0) continue;
        	// Open the file to see if it has the right format and amadid
        	// Open the txt file, and read lines 2 and 3 to get the amadid
    		BufferedReader in = new BufferedReader(new FileReader(exptdir+'\\'+FEfiles[i]));
    		tline = in.readLine();	// header types
    		tline = in.readLine();	// header names
    		if(tline == null) continue;	// Some other weird file
			String[] f = tline.split("\t");
			int amadidfield= -1;
        	for(int ifld=0;ifld<f.length;ifld++)
        		if(f[ifld].equals("FeatureExtractor_Barcode")) {
        			amadidfield= ifld;
        			break;
        		}
        	if(amadidfield < 0) {
        		in.close();	// It's not a valid FE file
        		continue;
        	}
    		tline = in.readLine();	// header values
    		in.close();
			f = tline.split("\t");
			int amdid= (int)Integer.valueOf(f[amadidfield].substring(2,7));
        	if(amadid == 0) amadid= amdid;
        	if(amadid == amdid) {
        		FEfiles[nfiles++]= FEfiles[i];	// It's a keeper
        		//System.out.format("%s\\%s is a keeper\n", exptdir, FEfiles[i]);
        	}
        }
        
        // If we found at least one file, read lines 9 and 10 to match the field names
        if(nfiles == 0) return(0);
        
       FEfiles= Arrays.copyOf(FEfiles,nfiles);
       Arrays.sort(FEfiles);
        
        // Simplify the filenames for reference. This is the 5 digits of the filename
        // following the amadid, and the last 4 chars before the .txt. If I can't parse the 
        // filenames this way, name the samples sample1, etc.
        arrayNames= new String[nfiles];
        for(int ifile=0;ifile<nfiles;ifile++) {
        	String fstr= FEfiles[ifile];
        	int ind = fstr.indexOf(Integer.toString(amadid));
        	int indno = fstr.indexOf(".txt");
        	if(ind < 0 || indno < 0) arrayNames[ifile]= String.format("Sample_%d",ifile);
        	else arrayNames[ifile]= String.format("S%s%s",fstr.substring(ind+5,ind+10),fstr.substring(indno-4,indno));
        }
        
        // Read the name fields from file 0 only - they're all the same.
  		BufferedReader in = new BufferedReader(new FileReader(exptdir+'\\'+FEfiles[0]));
   		for(int iskip=0;iskip<8;iskip++) tline = in.readLine();
   		tline = in.readLine();	// field types
   		fieldtypes = tline.split("\t");
   		tline = in.readLine();	// field names
   		fieldnames = tline.split("\t");
   		// Count lines in the file to preallocate space
   		while(true) {
   	   		tline = in.readLine();
	        if(tline == null) break;
	        else if(tline.isEmpty()) continue;
	        else nlines++;
   		}
		in.close();
		
		// Decide which columns we're keeping, for namefields and valuefields, and if they are the right types
		for(int iname=0;iname<fieldnames.length;iname++) {
			//System.out.format("Trying to match field %s\n", fieldnames[iname]);
			for(int inamefld=0;inamefld<nameFields.length;inamefld++)
				if(nameFields[inamefld].equals(fieldnames[iname])) {
					namecols[inamefld]= iname;
					nametypes[inamefld]= -1;
					for(int itype=0;itype<supportedTypes.length;itype++)
						if(supportedTypes[itype].equals(fieldtypes[iname])) {
							nametypes[inamefld]= itype;
							//System.out.format("Name field %s is type %s\n", nameFields[inamefld], supportedTypes[itype]);
							break;
						}
				}
			for(int ivalfld=0;ivalfld<valueFields.length;ivalfld++)
				if(valueFields[ivalfld].equals(fieldnames[iname])) valuecols[ivalfld]= iname;
		}
		
		// Complain if I can't find a field or if it's the wrong type.
		for(int inamefld=0;inamefld<nameFields.length;inamefld++)
			if(namecols[inamefld]==0 || nametypes[inamefld]<0)
				System.out.format("Error: Can't match name field %s to a supported type\n", nameFields[inamefld]);
		for(int ivalfld=0;ivalfld<valueFields.length;ivalfld++)
			if(valuecols[ivalfld]==0)
				System.out.format("Error: Can't match value field %s\n", valueFields[ivalfld]);
				
		// Preallocate by type
		featureNames= new String[nameFields.length][nlines];
		featureValues	= new float[valueFields.length][nfiles][nlines];

		// Read the requested names from the first file
		// *** Reporting them all as Strings, because Object[] is too complicated.
  		in = new BufferedReader(new FileReader(exptdir+'\\'+FEfiles[0]));
   		for(int iskip=0;iskip<10;iskip++) tline = in.readLine();
   		for(int iline=0;iline<nlines;iline++) {
   	   		tline = in.readLine();
	        if(tline == null) break;
	   		fieldvals = tline.split("\t");
	   		for(int icol=0;icol<namecols.length;icol++)
	   			if(nametypes[icol]==2) featureNames[icol][iline]= fieldvals[namecols[icol]];
	   			else featureNames[icol][iline]= fieldvals[namecols[icol]]; // (int)Integer.valueOf(fieldvals[namecols[icol]]);
   		}
		in.close();		
		
		// Finally, read the values of the value fields for all the files. Treat them all as floats.
		for(int ifile=0;ifile<nfiles;ifile++) {
	  		in = new BufferedReader(new FileReader(exptdir+'\\'+FEfiles[ifile]));
	   		for(int iskip=0;iskip<10;iskip++) tline = in.readLine();
	   		for(int iline=0;iline<nlines;iline++) {
	   	   		tline = in.readLine();
		        if(tline == null) break;
		   		fieldvals = tline.split("\t");
		   		for(int icol=0;icol<valuecols.length;icol++)
		   			if(valuecols[icol] > 0)
		   				featureValues[icol][ifile][iline]= (float)Float.valueOf(fieldvals[valuecols[icol]]);
	   		}
			in.close();		
		}		
		
        return(nfiles);
	} // public int ReadFEdata()
	
	/*
	 * Give the caller selective access to the data from a single field.
	 */
	public String[] FEnameFieldValue(String fieldname)
	{
		for(int ifld=0;ifld<nameFields.length;ifld++)
			if(fieldname.equals(nameFields[ifld])) return(featureNames[ifld]);
		return(null);
	}	// 	public String[] FEnameFieldValue(String fieldname)
	public String[] FEnameFieldValue(int ifld) {
		 if(ifld < featureNames.length) return(featureNames[ifld]);
		 return(null);
	}

	public float[][] FEvalueFieldValue(String fieldname)
	{
		for(int ifld=0;ifld<valueFields.length;ifld++)
			if(fieldname.equals(valueFields[ifld])) return(featureValues[ifld]);
		return(null);
	}	// 	public float[][] FEvalueFieldValue(String fieldname)
	public float[][] FEvalueFieldValue(int ifld) {
		return(featureValues[ifld]);
	}
	
	/*
	 * Compute the background signal from the negative control features.
	 * *** Requires "ControlType" to be one of the specified nameFields,
	 * and at least one of gBGSubSignal or rBGSubSignal to be one of the value fields.
	 * Ignore any NaNs.
	 * *** Could read FE's estimates from the header fields gNegCtrlAveBGSubSig, gNegCtrlSDevBGSubSig,
			rNegCtrlAveBGSubSig, andrNegCtrlSDevBGSubSig
.	 */
	public float[] backgroundLevels(String gsigfield, String rsigfield)
	{
		if(bkg.length == 2*FEfiles.length) return(bkg);
		int[] sigfields= {-1, -1};
		for(int ifld=0;ifld<valueFields.length;ifld++) {
			if(valueFields[ifld].equals(gsigfield)) sigfields[0]= ifld;
			if(valueFields[ifld].equals(rsigfield)) sigfields[1]= ifld;
		}
		if(sigfields[0] < 0 && sigfields[1]<0) return(null);
		int narrays= featureValues[0].length;
		int nfeats = featureValues[0][0].length;
		for(int ifld=0;ifld<nameFields.length;ifld++)
			if(nameFields[ifld].equals("ControlType")) {
				bkg= new float[2*narrays];
				bkgsd= new float[2*narrays];
				int nneg = 0;
				int[] imneg = new int[nfeats];
				for(int ifeat=0;ifeat<nfeats;ifeat++) {
					if(Integer.valueOf(featureNames[ifld][ifeat]) < 0) imneg[nneg++]= ifeat;
				}
				float[] keepers = new float[nneg];
				System.out.format("We found %d negative controls, %d arrays\n", nneg, narrays);
				for(int iar=0;iar<narrays;iar++) {
					if(sigfields[0]>=0) {
						int nkept= 0;
						for(int ifeat=0;ifeat<nneg;ifeat++) {
							float value=  featureValues[sigfields[0]][iar][imneg[ifeat]];
							if(!Float.isNaN(value)) keepers[nkept++]= value;
						}
						float[] keptvals = Arrays.copyOf(keepers, nkept);
						//System.out.format("Array[%d] filtered out %d/%d green outliers\n", iar, nneg-nkept, nneg);
						//System.out.format("keptvals:");
						//for(int ip=0;ip<12;ip++) System.out.format("\t%.4f", keptvals[ip]); 
						//System.out.format("\n");
						Arrays.sort(keptvals);
						if((nkept%2)==0)	bkg[iar]= 0.5F*(keptvals[nkept/2-1]+keptvals[nkept/2]);
						else 			bkg[iar]= keptvals[(nkept-1)/2];
						float iqr = keptvals[3*nkept/4-1] - keptvals[nkept/4-1];	// use the floor
						bkgsd[iar]= iqr/1.349F;	//  2*sqrt(2)*erfinv(.5)
						//System.out.format("iar=%d, negs= [%.4f,%.4f,%.4f,%.4f]\n", iar, keepers[nneg/4], keepers[nneg/2], keepers[nneg/2+1], keepers[3*nneg/4]);
					}
					if(sigfields[1]>=0) {
						int nkept= 0;
						for(int ifeat=0;ifeat<nneg;ifeat++) {
							float value=  featureValues[sigfields[1]][iar][imneg[ifeat]];
							if(!Float.isNaN(value)) keepers[nkept++]= value;
						}
						float[] keptvals = Arrays.copyOf(keepers, nkept);
						//System.out.format("Array[%d] filtered out %d/%d red outliers\n", iar, nneg-nkept, nneg);
						Arrays.sort(keptvals);
						if((nkept%2)==0)	bkg[narrays+iar]= 0.5F*(keptvals[nkept/2-1]+keptvals[nkept/2]);
						else 			bkg[narrays+iar]= keptvals[(nkept-1)/2];
						float iqr = keptvals[3*nkept/4-1] - keptvals[nkept/4-1];	// use the floor
						bkgsd[narrays+iar]= iqr/1.349F;	//  2*sqrt(2)*erfinv(.5)
					}
				}
				return(bkg);
			}
		return(null);
	}	// 	public float[] backgroundLevels()
	
		public float[] backgroundNoise(String gsigfield, String rsigfield)
		{
			if(bkgsd.length == 2*FEfiles.length) return(bkgsd);
			//System.out.format("bgsd is %d long, there are %d FEfiles\n", bkgsd.length, FEfiles.length);
			@SuppressWarnings("unused")
			float[] mybkg = backgroundLevels(gsigfield, rsigfield);	// Computes bkgsd as a side effect
			return(bkgsd);
		} 	// 	public float[] backgroundNoise()
		
		/*
		 * Filter the data in the given value array, using the specified value field a a filter.
		 * Although the filter field is nominally a float, it should be either 0 or 1. 
		 * Any values that correspond to nonzero values of the filterField will be replaced with NaN.
		 */
		public void filterValues(String valueField, String filterField)
		{
			int ivalField = -1, ifiltField = -1;
			
			for(int ifld=0;ifld<valueFields.length;ifld++) {
				if(valueFields[ifld].equals(valueField)) ivalField= ifld;
				if(valueFields[ifld].equals(filterField)) ifiltField= ifld;
			}
			if(ivalField<0 || ifiltField<0) {
				System.out.format("Error: either %s or %s wasn't read\n", valueField, filterField);
				return;
			}
			
			float[][] values= featureValues[ivalField];
			float[][] filter= featureValues[ifiltField];
			int nfeats= values[0].length;
			int narrays= values.length;
			
			for(int iar=0;iar<narrays;iar++)
				for(int ifeat=0;ifeat<nfeats;ifeat++)
					if(filter[iar][ifeat] != 0) {
						values[iar][ifeat]= Float.NaN; 
						//System.out.format("Nanning out feature[%d][%d]\n", iar, ifeat);
					}
			
			featureValues[ivalField]= values;
		}	// public void filterValues(String valueField, String filterField)
		
}
