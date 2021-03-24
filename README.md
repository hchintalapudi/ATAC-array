# Analysis of Surajit Dhara's paper (Pancreatic cancer prognosis is predicted by a novel ATAC-array technology for assessing chromatin accessibility):
This reository has the analysis scripts for ATACseq and ATACarray analysis of the research project published and linked [here](https://www.biorxiv.org/content/10.1101/2021.01.21.427604v1) 

# ATAC-seq analysis:
ATAC-seq analysis pipeline and processing scripts authored by Sagar Chhangawala (sagar.cornell@gmail.com)




# ATAC-array analysis: 
ATAC-array analysis java code authored by  Bo Curry and implemented by Himanshu Chintalapudi for Surajit Dhara's research project.<br/>
The pipeline (src/ATACanalysis/ATACanalysis.java) is used to read the FE files, normalize and average the region data, and classify the samples.<br/>
It produces a file *"ATACregionLRs.txt"* which can be further processed through downstream analysis R scripts for visualization.<br/>

Input FE files in the folder 'Ag_FE_data' and more info in config/README.txt


**Requirements:**<br/>
JAVA - version 1.8

**Usage:**<br/>
java  ATACanalysis.ATACanalysis  ATACanalysis/config/  Ag_FE_data/* 
