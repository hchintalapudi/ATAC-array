setwd("/media/sagar/DataDrive/Pancreatic_ATACseq/pancreaticcancer_atacseq_analysis/")
library(data.table)
library(GenomicFeatures)
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)
library(Vennerable)
library(pheatmap)
library(ggplot2)
library(plyr)
library(BiocParallel)
library(rtracklayer)
library(genefilter)
library(beanplot)
library(GenomicRanges)
library(RColorBrewer)

data.dir="Atlas_Final/"
outDir <- "Atlas_Final/"
plots.dir="Atlas_Final/plots/"
dir.create(plots.dir,recursive = T)
IDR_value=0.01

# create atlas.#################################################################################
# Inspired in EncodeEpigenomics/merged_peaks/merge_cell_type_peaks.R.

collapse.pairwise.celltype.peaks <- function (peaks1, peaks2, overlap.ratio.cutoff=0.75) {
  
  ## Function to find overlapping ratios
  find.overlap.ratio <- function (gr1, gr2) {
    
    ## Find Overlaps
    overlaps <- as.matrix (findOverlaps (gr1, gr2))
    
    ## Build ranges
    ranges <- cbind (start (gr1)[overlaps[,1]], end (gr1)[overlaps[,1]],
                     start (gr2)[overlaps[,2]], end (gr2)[overlaps[,2]])
    ranges <- t(apply (ranges, 1, sort))
    
    if(length(ranges==0))
      return (list (overlap.ratio = 0,
                    ranges = ranges,
                    overlaps = overlaps,
                    best.common = 0))
    
    ## Min widths
    widths <- pmin (width (gr1)[overlaps[,1]], width (gr2)[overlaps[,2]])
    
    ## Overlap ratios
    overlap.ratio <- (ranges[,3] - ranges[,2])/widths
    
    ## Find best in both  
    best.gr1 <- tapply (1:nrow (overlaps), overlaps[,1], function (x) { x[overlap.ratio[x]==max(overlap.ratio[x])][1] } )
    best.gr2 <- tapply (1:nrow (overlaps), overlaps[,2], function (x) { x[overlap.ratio[x]==max(overlap.ratio[x])][1] } )
    common <- intersect (best.gr1, best.gr2)
    
    return (list (overlap.ratio = overlap.ratio,
                  ranges = ranges,
                  overlaps = overlaps,
                  best.common = common))
  }
  
  ## Reset metadata
  values (peaks1) <- values (peaks2) <- NULL
  
  ## Overlapping peaks which exceed overlap.ratio.cutoff
  or <- find.overlap.ratio (peaks1, peaks2)
  common <- or$best.common[or$overlap.ratio[or$best.common] > overlap.ratio.cutoff]
  ## Create GRanges object with common regions
  union <- GRanges (seqnames (peaks1)[or$overlaps[common,1]],
                    IRanges (or$ranges[common,2], or$ranges[common,3]))
  
  ## Determine disjoint peaks
  peaks1 <- peaks1[countOverlaps (peaks1, union) == 0]
  peaks2 <- peaks2[countOverlaps (peaks2, union) == 0]
  
  ## Other overlapping peaks
  or <- find.overlap.ratio (peaks1, peaks2)
  common <- or$best.common
  ## Create GRanges with union of the regions
  union <- c(union,GRanges (seqnames (peaks1)[or$overlaps[common,1]],
                            IRanges (or$ranges[common,1], or$ranges[common,4])))
  
  ## Non overlapping peaks
  union <- c(union, peaks1[countOverlaps (peaks1, union) == 0])
  union <- c(union, peaks2[countOverlaps (peaks2, union) == 0])
  
  return (union)
}

# assumes function "collapse.pairwise.celltype.peaks".
mergeCTpeaks <- function(pkLst, GRfile, DTfile) {
  
  combined.peaks <- collapse.pairwise.celltype.peaks (pkLst[[1]], pkLst[[2]])
  for (pkCtr in 3:length(pkLst)) {
    combined.peaks <- collapse.pairwise.celltype.peaks (combined.peaks, pkLst[[pkCtr]])
  }
  
  ## annotate cell type where peak was called.
  
  for (cellType in names(pkLst)) {
    
    overlaps <- countOverlaps (combined.peaks, pkLst[[cellType]])
    mcols (combined.peaks)[,paste0(cellType, ".peak")] <- 0
    mcols (combined.peaks)[,paste0(cellType, ".peak")][overlaps > 0] <- 1
    
  }
  
  names(combined.peaks) <- 1:length(combined.peaks)
  combined.peaks$pattern <- apply(as.matrix(mcols(combined.peaks)), 1, function(row) {
    argsToPaste0 <- as.list(row)
    do.call(paste0, argsToPaste0)
  })
  
  combined.peaks.DT <- data.table(pk = names(combined.peaks),
                                  chr = as.character(seqnames(combined.peaks)),
                                  start = as.integer(start(combined.peaks)),
                                  end = as.integer(end(combined.peaks)),
                                  pattern = combined.peaks$pattern)
  
  saveRDS(combined.peaks, file = GRfile)
  saveRDS(combined.peaks.DT, file = DTfile)
  
}

#Cutoff of 25 M reads.
All_peaks=list.files(path=".",pattern="idrValues.txt$",recursive=T)

#Create Atlas with IDR peaks. Only Samples with > 19180 peaks (Determined by 1st quantile of # reproducible peaks /sample)
pkLst=list()
PTPks_names=c()
rep_peaks=data.frame(PT=character(),Peaks=integer())

for(i in All_peaks){
  Sample_name=unlist(strsplit(i,"/",fixed=T))[2]
  PT_num=unlist(strsplit(as.character(unlist(strsplit(i,"/",fixed=T))[2]),"_",fixed=T))[4]
  PT <- read.table(i)
  idr.results <- GRanges (PT[,'V1'], IRanges (PT[,'V2'], PT[,'V3']), IDR = 2^(PT[,'V5']/-125))
  PTPks <- idr.results[idr.results$IDR <= IDR_value]
  print(length(PTPks))
  rep_peaks=rbind(rep_peaks,data.frame(PT=PT_num,Peaks=length(PTPks)))
  if (length(PTPks) > 19180){
    PTPks_names=c(PTPks_names,PT_num)
    pkLst=c(pkLst,PTPks)
  } else {
    print(paste0(i," < 19,180"))
  }
}
names(pkLst) <- PTPks_names

mergeCTpeaks(pkLst = pkLst, 
             GRfile = paste0(data.dir, "combined.peaks.",IDR_value,".rds"), 
             DTfile = paste0(data.dir, "combined.peaks.DT.",IDR_value,".rds"))

combined.peaks <- readRDS(paste0(data.dir, "combined.peaks.",IDR_value,".rds"))
combined.peaks.DT <- readRDS(paste0(data.dir, "combined.peaks.DT.",IDR_value,".rds"))

plot.df <- data.frame(x=factor(rep_peaks$PT[(order(rep_peaks$Peaks))],levels=rep_peaks$PT[(order(rep_peaks$Peaks))]),y=rep_peaks$Peaks[(order(rep_peaks$Peaks))])
ggplot(plot.df,aes(x=x,y=y))+geom_point(size=2)+ylab("# Reproducible Peaks")+xlab("")+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust=.5))+geom_hline(yintercept=19180)
ggsave(paste0("Figures/Fig_S3d.pdf"),height=10,width=12)
write.table(plot.df,file="Source_tables/Fig_S3d.txt",sep="\t",row.names=F,col.names=T,quote=F)

####################################################
# annotate with genes and genomic features.

# Here's the protocol:
# If 2000 or less distance to TSS, promoter and assign corresponding gene.
# Else, if doesn't overlap with transcript, intergenic and assign closest gene.
# Else, if overlaps by more than 50% with exon, exon.
# Else, intron.

txdb <- makeTxDbFromGFF("~/Documents/transcriptome/hg19/Homo_sapiens.GRCh37.87.PConly.gtf",format="gtf",dataSource = "ensembl",organism = "Homo sapiens")

# Promoter peaks.
#################
# => WARNING: I might be missing some intron guys by not looking specifically
# at the 5'UTR.

transcripts = transcriptsBy(txdb)
transcriptsList = unlist(transcripts)
tssGR = resize(transcriptsList,1)

nearestTss <- nearest(x = combined.peaks,subject = tssGR)
distToNearestTss <- distance(x = combined.peaks, y = tssGR[nearestTss])
nearestTssTab <- data.table(peak_id = names(combined.peaks), transcript_id = tssGR[nearestTss]$tx_name, dist = distToNearestTss)

matchTranscript2PeakTab <- nearestTssTab[dist <= 2000, list(peak_id, transcript_id, annot = "promoter")]


# Intergenic peaks.
###################
overlapWithTranscr <- subsetByOverlaps(combined.peaks[!names(combined.peaks) %in% matchTranscript2PeakTab$peak_id], transcriptsList)
intergPks <- combined.peaks[!names(combined.peaks) %in% union(matchTranscript2PeakTab$peak_id, names(overlapWithTranscr))]

nearestTransc <- nearest(x = intergPks, subject = transcriptsList)
matchTranscript2PeakTab <- rbind(matchTranscript2PeakTab, data.table(peak_id = names(intergPks),
                                                                     transcript_id = transcriptsList[nearestTransc]$tx_name, annot = "intergenic"))



# Exon peaks.
#############
exons = exonsBy(txdb, by='tx', use.names=T)
exonsList = unlist(exons)
exonsList = exonsList[grep("_", as.character(seqnames(exonsList)), invert=T)]
names(exonsList) = paste(names(exonsList), exonsList$exon_id, sep="_")


ovlWithExon <- findOverlaps(overlapWithTranscr, exonsList)
ovlWithExonTab <- data.table(peak_id = names(overlapWithTranscr[queryHits(ovlWithExon)]),
                             exon_id = names(exonsList[subjectHits(ovlWithExon)]))


intersectRanges <- pintersect(overlapWithTranscr[ovlWithExonTab$peak_id], exonsList[ovlWithExonTab$exon_id])
ovlWithExonTab$intersectWidth <- width(intersectRanges)
ovlWithExonTab$peak_width <- width(overlapWithTranscr[ovlWithExonTab$peak_id])
ovlWithExonTab$intersectRatio <- ovlWithExonTab$intersectWidth/ovlWithExonTab$peak_width
exonPksTab <- ovlWithExonTab[, list(exon_id = exon_id[which.max(intersectRatio)], intersectRatio = max(intersectRatio)), by = peak_id][intersectRatio >= .5]
exonPksTab$transcript_id <- unlist(lapply(names(exonsList[exonPksTab$exon_id]), function(x){paste(strsplit(x, "_")[[1]][1:2], collapse="_")}))
exonPks <- overlapWithTranscr[exonPksTab$peak_id]
exonPks$transcript_id <- exonPksTab$transcript_id

matchTranscript2PeakTab <- rbind(matchTranscript2PeakTab,
                                 data.table(peak_id = names(exonPks), transcript_id = exonPks$transcript_id, annot = "exon"))


# Intron peaks.
###############
intronPks <- combined.peaks[setdiff(names(combined.peaks),matchTranscript2PeakTab$peak_id)]

nearestTransc <- nearest(x = intronPks, subject = transcriptsList)
matchTranscript2PeakTab <- rbind(matchTranscript2PeakTab,
                                 data.table(peak_id = names(intronPks), transcript_id = transcriptsList[nearestTransc]$tx_name, annot = "intron"))



geneLookup = read.table("hg19_Ensembl_conversion.txt", sep="\t", stringsAsFactors=F,header=T)
geneLookup = geneLookup[,c(2,3)]
geneLookup = unique(geneLookup[geneLookup[,1] != "n/a",])

matchTranscript2PeakTab$symbol = unlist(lapply(matchTranscript2PeakTab$transcript_id, function(x){
  paste(geneLookup[geneLookup[,1] == x,2], collapse=",")
}))

temp <- as.data.frame(combined.peaks)
temp$seqnames <- paste0("chr",temp$seqnames)
temp$seqnames <- sub("chrMT","chrM",temp$seqnames)
combined.peaks <- GRanges (temp$seqnames, IRanges (temp$start, temp$end), temp$strand,temp[,grep(".peak",colnames(temp))], pattern=temp$pattern)
names(combined.peaks) <- 1:length(combined.peaks)

setkey(matchTranscript2PeakTab, peak_id)
combined.peaks$annot <- matchTranscript2PeakTab[names(combined.peaks)]$annot
combined.peaks$transcript_id <- matchTranscript2PeakTab[names(combined.peaks)]$transcript_id
combined.peaks$symbol <- matchTranscript2PeakTab[names(combined.peaks)]$symbol
combined.peaks.DT$annot <- matchTranscript2PeakTab[combined.peaks.DT$pk]$annot
combined.peaks.DT$transcript_id <- matchTranscript2PeakTab[combined.peaks.DT$pk]$transcript_id
combined.peaks.DT$symbol <- matchTranscript2PeakTab[combined.peaks.DT$pk]$symbol

saveRDS(combined.peaks, file = paste0(data.dir, "combined.peaks.",IDR_value,".rds"))
saveRDS(combined.peaks.DT, file =paste0(data.dir, "combined.peaks.DT.",IDR_value,".rds"))

###############################################################
# alright, let's do counts.
names(combined.peaks) <- paste0("Peak",names(combined.peaks))
combined.peaks.DT[,pk:=paste0("Peak",pk)]

# List of bam file to count
bam.files = list.files(".", pattern="*.rmdup.bam$", full.names = T,recursive=T)
bam.files <- bam.files[grep(regex,bam.files,invert = T)]

do.counts = function(combined.peaks, file){
  peak.counts = tryCatch({
    bamFiles <- BamFileList(file, yieldSize=2000000)
    names(bamFiles) <- unlist(strsplit(file,split="/"))[2]
    counts = summarizeOverlaps(combined.peaks, bamFiles, "IntersectionNotEmpty", ignore.strand=TRUE)
  }, warning = function(war){
    print(paste("Error in summarizeOverlaps:", war))
    print("Trying to read BAM directly")
    
    param <- ScanBamParam (which = combined.peaks, flag = scanBamFlag (isProperPair = TRUE), tag=c('IH', 'NH'))
    tags <- readBamGappedAlignments (bam.files, param=param)
    counts = summarizeOverlaps(combined.peaks, tags, "IntersectionNotEmpty", ignore.strand=TRUE, obeyQname=TRUE)
    return(counts)
  }, error = function(err){
    print(paste("Error reading BAM:", err))
  })
  save(counts, file=paste0(outDir,unlist(strsplit(file,split="/"))[2],".counts_",IDR_value,".rds"))
}

register(MulticoreParam(workers=10))
bplapply(bam.files, do.counts, BPPARAM=MulticoreParam(), combined.peaks=combined.peaks)

for(file in bam.files){
  load(paste0(outDir,unlist(strsplit(file,split="/"))[2],".counts_",IDR_value,".rds"))
  command = paste0("combined.peaks$",unlist(strsplit(file,split="/"))[2], " = unlist(assays(counts)$counts)")
  eval(parse(text=command))
}

combined.peaks.DT = readRDS(paste0(outDir, "combined.peaks.DT.",IDR_value,".rds"))
combined.peaks = readRDS(paste0(outDir, "combined.peaks.",IDR_value,".rds"))

setkey(combined.peaks.DT, pk)
combined.peaks.DT <- cbind(combined.peaks.DT[names(combined.peaks)], as.data.frame(combined.peaks)[,grep("^Sample_ATAC",colnames( as.data.frame(combined.peaks)))])
colnames(combined.peaks.DT) = gsub("Sample_ATAC_", "", colnames(combined.peaks.DT))

saveRDS(combined.peaks, file = paste0(data.dir, "combined.peaks.",IDR_value,".rds"))
saveRDS(combined.peaks.DT, file =paste0(data.dir, "combined.peaks.DT.",IDR_value,".rds"))

annot <- as.data.frame(combined.peaks.DT[, list(cnt = .N), by = annot][order(cnt, decreasing = TRUE)])
annot$percent <- round(annot$cnt/sum(annot$cnt)*100)
annot$label <- paste0(annot$annot,"\n",annot$cnt," peaks","\n",annot$percent,"%")

pdf("Figures/Fig_S3b.pdf")
pie(annot$cnt,annot$label,col=brewer.pal(length(annot$annot),"Set1"))
dev.off()
write.table(annot[,c("annot","cnt","percent")],file="Source_tables/Fig_S3b.txt",sep="\t",row.names=F,col.names=T,quote=F)

#beanplot
combined.peaks.df <- as.data.frame(combined.peaks)
combined.peaks.df <- combined.peaks.df[which(combined.peaks.df$annot=="exon"),]
beanplot.df <- data.frame(Sum=rowSums(combined.peaks.df[,grep(".peak",colnames(combined.peaks.df))]),annot="Exon",row.names = NULL)

combined.peaks.df <- as.data.frame(combined.peaks)
combined.peaks.df <- combined.peaks.df[which(combined.peaks.df$annot=="intergenic"),]
beanplot.df <- rbind(beanplot.df,data.frame(Sum=rowSums(combined.peaks.df[,grep(".peak",colnames(combined.peaks.df))]),annot="Intergenic",row.names = NULL))

combined.peaks.df <- as.data.frame(combined.peaks)
combined.peaks.df <- combined.peaks.df[which(combined.peaks.df$annot=="intron"),]
beanplot.df <- rbind(beanplot.df,data.frame(Sum=rowSums(combined.peaks.df[,grep(".peak",colnames(combined.peaks.df))]),annot="Intron",row.names = NULL))

combined.peaks.df <- as.data.frame(combined.peaks)
combined.peaks.df <- combined.peaks.df[which(combined.peaks.df$annot=="promoter"),]
beanplot.df <- rbind(beanplot.df,data.frame(Sum=rowSums(combined.peaks.df[,grep(".peak",colnames(combined.peaks.df))]),annot="Promoter",row.names = NULL))

combined.peaks.df <- as.data.frame(combined.peaks)
beanplot.df <- rbind(beanplot.df,data.frame(Sum=rowSums(combined.peaks.df[,grep(".peak",colnames(combined.peaks.df))]),annot="All",row.names = NULL))

pdf(paste0("Figures/Fig_S3c.pdf"),height=8,width=10)
beanplot(Sum ~ annot, data=beanplot.df,bw="nrd",what=c(0,1,0,0),col=list(c("chartreuse","white","white","white"),c("brown1","white","white","white"),c("deepskyblue","white","white","white"),c("deeppink","white","white","white"),c("chocolate1","white","white","white")))
dev.off()
write.table(beanplot.df,file="Source_tables/Fig_S3c.txt",sep="\t",row.names=F,col.names=T,quote=F)



