# ATAC-array analysis:
This pipeline (src/ATACanalysis/ATACanalysis.java) is used to read the FE files, normalize and average the region data, and classify the samples. 
It produces a file *"ATACregionLRs.txt"* which can be further processed through downstream analysis R scripts for visualization.

**Requirements:**<br/>
JAVA - version 1.8

**Usage:**<br/>
java ATACanalysis.ATACanalysis ATACanalysis/config/ Ag_FE_data/* 
