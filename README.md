# ATAC-array analysis:
This pipeline (src/ATACanalysis/ATACanalysis.java) is used to read the FE files, normalize and average the region data, and classify the samples.<br/>
It produces a file *"ATACregionLRs.txt"* which can be further processed through downstream analysis R scripts for visualization.<br/>

Input FE files in the folder 'Ag_FE_data' and more info in config/README.txt


**Requirements:**<br/>
JAVA - version 1.8

**Usage:**<br/>
java  ATACanalysis.ATACanalysis  ATACanalysis/config/  Ag_FE_data/* 
