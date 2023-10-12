##instalations 
BiocManager::install("plyranges", force=TRUE)
options(timeout = 800) 
register(SerialParam())
gc()
## libraries 
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rbowtie2)
library(ggplot2)
library(DiffBind)
library(GenomicAlignments)
library(dplyr)
library(BiocParallel)
library(DiffBind)
library(profileplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(profileplyr)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(csaw)
library(PCAtools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


##### set dir
setwd("D:/MiguelW12/Documents/gsd1a_atac_analysis/peak_data")
output_dir <- "res"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
### load samples to dba
##GSD1a
GSD1a_vs_HC<-dba.peakset(NULL,
                         peaks="762.sorted_no_mt_noDup/762.sorted_no_mt_noDup_peaks_peaks.xls",
                         peak.caller="macs", sampID="GSD1a1",condition = "GSD1a", replicate=1,
                         bamReads = "bam_ready_qc/762.sorted_no_mt_noDup.bam")
GSD1a_vs_HC<-dba.peakset(GSD1a_vs_HC,
                         peaks="6894.sorted_no_mt_noDup/6894.sorted_no_mt_noDup_peaks_peaks.xls",
                         peak.caller="macs", sampID="GSD1a2",condition = "GSD1a", replicate=2,
                         bamReads = "bam_ready_qc/6894.sorted_no_mt_noDup.bam")
##HC                     
GSD1a_vs_HC<-dba.peakset(GSD1a_vs_HC,
                         peaks="498.sorted_no_mt_noDup/498.sorted_no_mt_noDup_peaks_peaks.xls",
                         peak.caller="macs", sampID="HC1",condition = "HC", replicate=1,
                         bamReads = "bam_ready_qc/498.sorted_no_mt_noDup.bam")

GSD1a_vs_HC<-dba.peakset(GSD1a_vs_HC,
                         peaks="17507.sorted_no_mt_noDup/17507.sorted_no_mt_noDup_peaks_peaks.xls",
                         peak.caller="macs", sampID="HC2",condition = "HC", replicate=2,
                         bamReads = "bam_ready_qc/17507.sorted_no_mt_noDup.bam")
## verify sample load
plot(GSD1a_vs_HC)
##set thresh
GSD1a_vs_HC$config$th = 0.1
## counts 
GSD1a_vs_HC_counts <- dba.count(GSD1a_vs_HC, bParallel = TRUE, score=DBA_SCORE_NORMALIZED,
                                bUseSummarizeOverlaps=TRUE, summits=250)
GSD1a_vs_HC_counts$config$th = 0.1
######QC
peakdata <- dba.show(GSD1a_vs_HC_counts)$Intervals
peakdata
info <- dba.show(GSD1a_vs_HC_counts)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes
##
dba.plotVenn(GSD1a_vs_HC_counts,GSD1a_vs_HC_counts$masks$GSD1a, main = "Open chromatic region overlaps in GSD1a")
dba.plotVenn(GSD1a_vs_HC_counts,GSD1a_vs_HC_counts$masks$HC, main = "Open chromatic region overlaps in HC")
##occupancy
HCvsGSD1a_occupancy<-dba.peakset(GSD1a_vs_HC, consensus = DBA_CONDITION,minOverlap = 0.5)
HCvsGSD1a_occupancy
##pca
sset <- dba(HCvsGSD1a_occupancy,bSummarizedExperiment=TRUE)
pca_dataFrame <- as.data.frame(assay(sset))
rownames(pca_dataFrame) <- rowRanges(sset)
p <- pca(pca_dataFrame, metadata = colData(sset), removeVar = 0.5)
biplot(p, showLoadings = FALSE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

## normaliation
GSD1a_vs_HC_counts_n <- dba.normalize(GSD1a_vs_HC_counts, library = DBA_LIBSIZE_PEAKREADS,
                                      method = DBA_DESEQ2)

normlibs <- cbind(FullLibSize=GSD1a_vs_HC_counts_n$norm$DESeq2$lib.sizes, NormFacs=GSD1a_vs_HC_counts_n$norm$DESeq2$norm.facs,
                  NormLibSize=round(GSD1a_vs_HC_counts_n$norm$DESeq2$lib.sizes/GSD1a_vs_HC_counts_n$norm$DESeq2$norm.facs))
rownames(normlibs) <- info$ID
normlibs
##differntial
GSD1a_vs_HC_contrast <- dba.contrast(GSD1a_vs_HC_counts_n,categories=DBA_CONDITION,minMembers = 2)
GSD1a_vs_HC_contrast$config$th = 0.1
dif_GSD1a_vs_HC_counts <- dba.analyze(GSD1a_vs_HC_contrast)
dif_GSD1a_vs_HC_counts
dba.show(dif_GSD1a_vs_HC_counts, bContrasts=TRUE)
res_deseq_hc_vs_gsd <- dba.report(dif_GSD1a_vs_HC_counts, contrast = 1, method=DBA_DESEQ2, th=1)
write.csv(as.data.frame(res_deseq_hc_vs_gsd), file.path(output_dir, "res_deseq_hcvsgsd_peaks_DF_noanno.csv"), row.names = TRUE)
### anno 
GenomeInfoDb::seqlevels(res_deseq_hc_vs_gsd)
GenomeInfoDb::seqlevels(txdb)
##convert chr names 
newnames <- paste0(c("chr1","chr2","chr3","chr4", "chr5","chr6","chr7", "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT"))
names(newnames) <- paste0(c("NC_000001.11","NC_000002.12","NC_000003.12","NC_000004.12", "NC_000005.10","NC_000006.12","NC_000007.14", "NC_000008.11","NC_000009.12","NC_000010.11","NC_000011.10","NC_000012.12","NC_000013.11","NC_000014.9","NC_000015.10","NC_000016.10","NC_000017.11","NC_000018.10","NC_000019.10","NC_000020.11","NC_000021.9","NC_000022.11","NC_000023.11","NC_000024.10","NC_012920.1"))
##
res_deseq_hc_vs_gsd_level <- renameSeqlevels(res_deseq_hc_vs_gsd,newnames)
GenomeInfoDb::seqlevels(res_deseq_hc_vs_gsd_level)
GenomeInfoDb::seqlevels(txdb)
length(seqlevels(res_deseq_hc_vs_gsd))
length(seqlevels(txdb))
res_deseq_hc_vs_gsd_level2 <- keepSeqlevels(res_deseq_hc_vs_gsd_level, standardChromosomes(res_deseq_hc_vs_gsd_level)[1:22], pruning.mode = "coarse")
seqlevels(res_deseq_hc_vs_gsd_level2)
# Annotate peaks
anno_df <- annotatePeak(res_deseq_hc_vs_gsd_level2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db')
anno_df
p <- as.GRanges(anno_df)
p
write.csv(as.data.frame(anno_df), file.path(output_dir, "res_deseq_hcvsgsd_peaks_DF_withanno.csv"), row.names = TRUE)
write.csv(as.data.frame(p), file.path(output_dir, "res_deseq_hcvsgsd_peaks_DF_withanno-granges.csv"), row.names = TRUE)

#vis
plotAnnoPie(anno_df)
plotAnnoBar(anno_df)
upsetplot(anno_df)

## vis res
dba.plotMA(dif_GSD1a_vs_HC_counts, method = DBA_DESEQ2, bNormalized = TRUE, sub = "HCvsGSD1a")
dba.plotVolcano(dif_GSD1a_vs_HC_counts, method = DBA_DESEQ2)
##
sum(res_deseq_hc_vs_gsd$Fold<0)
sum(res_deseq_hc_vs_gsd$Fold>0)

pvals <- dba.plotBox(dif_GSD1a_vs_HC_counts, contrast = 1)
pvals

### GSEA
gsea_df <- as.data.frame(anno_df)
original_gene_list <- gsea_df$Fold
names(original_gene_list) <- gsea_df$ENSEMBL
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gsea_res <- gseGO(geneList=gene_list, 
                  ont ="ALL", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, 
                  maxGSSize = 1000, 
                  pvalueCutoff = 0.1, 
                  verbose = TRUE, 
                  OrgDb = 'org.Hs.eg.db',
                  pAdjustMethod = "BH"
)

#save file
write.csv(as.data.frame(gsea_res), file.path(output_dir, "gseago-gsd1a-atac.csv"), row.names = TRUE)
## subtract res and vis
id_list <- c("GO:0099503",
             "GO:0005768",
             "GO:1901135",
             "GO:0005764",
             "GO:0030163",
             "GO:0031667",
             "GO:0006914",
             "GO:0003682",
             "GO:0000302",
             "GO:0007005", 
             "GO:0005740", 
             "GO:0005759", 
             "GO:0006325") 

gsea_res@result = gsea_res@result[gsea_res@result$ID %in% id_list, ]
## vis
edo <- setReadable(gsea_res, 'org.Hs.eg.db', 'ENSEMBL')
cnetplot(edo, node_label="all", categorySize="pvalue", cex_label_gene = 0.3, cex_label_category = 1,
         cex_category = 1, cex_gene = 0.3, colorEdge = TRUE, showCategory=10) 

##gviz ( for sample plotting generate granges per grouping+ chr plotting available)
dTrack <- DataTrack(p, name = "uniform")
colnames(mcols(p))
chromosome(dTrack) <- "chr9"
dTrack <- dTrack[1:4, ]
plotTracks(dTrack)
atrack <- AnnotationTrack(anno_range, name = "CpG")
##
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr9")
plotTracks(ideoTrack, from = 127397138, to = 127407855)
##
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr2")
plotTracks(ideoTrack, from = 127397138, to = 127407855)

##not used 
plotTracks(dtrack, type = "boxplot", showSampleNames = TRUE, 
           cex.sampleNames = 0.6)

plotTracks(dTrack, type = "histogram", groups = rownames(values(dTrack)),
           showSampleNames = TRUE, cex.sampleNames = 0.6, separator = 1)





#################################################################################







