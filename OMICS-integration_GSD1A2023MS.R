### analysis code in this script written by Hadasa Kaufman
##libraries
library (rio)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggVennDiagram)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(tidytree)
## data import
metil_data<-import("G:/My Drive/PhD/Uri_Data/DMP_data_GSD1A.csv") 
atac_data<-import("G:/My Drive/PhD/Uri_Data/peak_anno_peaks_DF.csv")

#count representations per gene api
gdata<-atac_data %>%
  count(SYMBOL) 
nrow(gdata[gdata$n==1,])
nrow(gdata[gdata$n>1,])

hist(gdata$n[gdata$n<=30],breaks=60, ) 


#######order introns annotations
atac_intron<-atac_data %>%
  filter(str_detect(annotation, "intron"))

atac_first_intron<-atac_intron %>%
  filter(str_detect(annotation, "1 of"))

atac_con_intron<-atac_intron %>%
  filter(!str_detect(annotation, "1 of"))

atac_first_intron[,"annotation"]<-"intron 1"

atac_con_intron[,"annotation"]<-"intron >=2"

atac_not_intron<-atac_data %>%
  filter(!str_detect(annotation, "intron"))

####order exon annotations

atac_not_intron_or_exon<-atac_not_intron %>%
  filter(!str_detect(annotation, "Exon"))

atac_exon<-atac_not_intron %>%
  filter(str_detect(annotation, "Exon"))
atac_exon[,"annotation"]<-'exon'

#merge all the sub datasets
atac_data_new_annotations<-rbind(atac_not_intron_or_exon,atac_exon,atac_con_intron,atac_first_intron)

#create annotation dataset
atac_annotation_chart<-atac_data_new_annotations%>%
  count(annotation)


names(atac_annotation_chart)<-c("annotation","count")


#create annotatins distrebution plot
atac_annotation_chart[,"value"]<-round((atac_annotation_chart[,"count"]*100)/sum(atac_annotation_chart[,"count"]),1)

df2 <- atac_annotation_chart %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

anno_plot<-ggplot(df2, aes(x = "" , y = value, fill = fct_inorder(annotation))) +
  geom_bar(stat="identity", width=1) +
  geom_col(color = "black",width = 1) +
  coord_polar("y", start=0)+
  #geom_text(aes(label = count),position = position_stack(vjust = 0.5), color = "white", size=3) +
geom_label_repel(data = df2,
                 aes(y = pos, label = paste0(value, "%")),
                 size = 2.5, nudge_x = 1, show.legend = FALSE)+
  theme_void() # remove background, grid, numeric labels

# Print plots to a pdf file
pdf("G:/My Drive/PhD/Uri_Data/plots/atac_annotations_distrebution.pdf")
print(anno_plot)     # Plot 1 --> in the first page of PDF
dev.off() 

#atac_relevant_annotations_data<-atac_data_new_annotations %>%
#  filter(annotation %in% c('Promoter (1-2kb)','Promoter (2-3kb)','Promoter (<=1kb)','intron 1'))

atac_relevant_annotations_data<-atac_data_new_annotations %>%
  filter(annotation!="Distal Intergenic" )

atac_relevant_annotations_data_by_gene<-atac_relevant_annotations_data %>%
  count(SYMBOL) 
nrow(gdata[atac_relevant_annotations_data_by_gene$n==1,])
nrow(gdata[atac_relevant_annotations_data_by_gene$n>1,])

hist(atac_relevant_annotations_data_by_gene$n[atac_relevant_annotations_data_by_gene$n<=30],breaks=60, ) 

hist(atac_relevant_annotations_data$FDR[atac_relevant_annotations_data$FDR<1],breaks = 50) 
sum(atac_relevant_annotations_data$FDR==0)

atac_relevant_annotations_data['FDR_']<-ifelse(atac_relevant_annotations_data$FDR<0.1,"<0.1",">=0.1")
ggplot(atac_relevant_annotations_data,aes(Fold,-log10(FDR)),color='black') +
         geom_point()+
          geom_point(atac_relevant_annotations_data[atac_relevant_annotations_data$FDR<0.1&atac_relevant_annotations_data$Fold>0,],mapping=aes(Fold,-log10(FDR)),color="red")+
          geom_point(atac_relevant_annotations_data[atac_relevant_annotations_data$FDR<0.1&atac_relevant_annotations_data$Fold<0,],mapping=aes(Fold,-log10(FDR)),color="blue")

atac_relevant_annotations_data<-atac_relevant_annotations_data[!is.na(atac_relevant_annotations_data$SYMBOL),]

atac_sagn<-atac_relevant_annotations_data[atac_relevant_annotations_data$FDR<0.1,]

######################3

write.csv(atac_sagn,"G:/My Drive/PhD/Uri_Data/sagn_atac.csv")

#create data_frame of the most extream log2FC in the atac data
pgdata_atac<-atac_relevant_annotations_data %>%
  group_by(SYMBOL) %>%
  mutate(atac_max=ifelse(min(Fold)*(-1)<max(Fold),(1),(-1))*max(abs(Fold)),min_p_atac=min(FDR))
nrow(pgdata_atac)

pgdata_atac_extream<-pgdata_atac[!duplicated(pgdata_atac$SYMBOL),]



###########################integration of met+atac#######################
##########################################################################################

##order met data
library("reshape")
metil_data<-import("G:/My Drive/PhD/Uri_Data/DMP_data_GSD1A.csv")

names(metil_data)

#separation of the genes list
metil_data[c('First_gene', 'gene2','gene3','gene4','gene5','gene6','gene7','gene8','gene9','gene10')] <- str_split_fixed(metil_data$UCSC_RefGene_Name, ';', 10)
names(metil_data)

#reorder the dataframe
metil_data_melt<-melt(metil_data[c('First_gene', 'gene2','gene3','gene4','gene5','gene6','gene7','gene8','gene9','gene10','logFC','adj.P.Val','Name')],id=c('logFC','adj.P.Val','Name'))
metil_data_melt<-metil_data_melt[c('logFC','adj.P.Val','value','Name')]
names(metil_data_melt)=c('logFC','adj.P.Val','SYMBOL','Name')

metil_data_melt<-metil_data_melt[!duplicated(metil_data_melt$SYMBOL)|!duplicated(metil_data_melt$adj.P.Val)|!duplicated(metil_data_melt$logFC),]

#create significant names vector
metil_data_names_sagn<-metil_data_melt[metil_data_melt$adj.P.Val<0.1,]$Name
length(metil_data_names_sagn)

metil_data_melt_na<-metil_data_melt[metil_data_melt$SYMBOL!="",]

#count rows per gene (SYMBOL)
gdata<-metil_data_melt_na %>%
  count(SYMBOL) 
nrow(gdata[gdata$n==1,])
nrow(gdata[gdata$n>1,])

#hists
hist(gdata$n[gdata$n<=30],breaks=60, ) 
hist(metil_data_melt_na[metil_data_melt_na$adj.P.Val<0.2,"adj.P.Val"],breaks = 50) 

metil_data_melt_na['mat_FDR_cat']<-ifelse(metil_data_melt_na$adj.P.Val<0.1,"<0.1",">=0.1")
ggplot(metil_data_melt_na,aes(logFC,-log10(adj.P.Val)),color='black') +
  geom_point()+
  geom_point(metil_data_melt_na[metil_data_melt_na$adj.P.Val<0.1&metil_data_melt_na$logFC>0,],mapping=aes(logFC,-log10(adj.P.Val)),color="red")+
  geom_point(metil_data_melt_na[metil_data_melt_na$adj.P.Val<0.1&metil_data_melt_na$logFC<0,],mapping=aes(logFC,-log10(adj.P.Val)),color="blue")

met_sagn<-metil_data_melt_na[metil_data_melt_na$adj.P.Val<0.1,]


#merge atac met
merge_atac_met<-merge(metil_data_melt_na,atac_relevant_annotations_data,by="SYMBOL",how="inner")

ggplot(merge_atac_met[merge_atac_met$FDR<1,],aes(x=Fold,y=logFC))+
  geom_point()+xlab("atac seq log2FC")+ylab("methylation log2FC")


#create data_frame of the most extream log2FC in the metylation data
pgdata_met<-metil_data_melt_na %>%
  group_by(SYMBOL) %>%
  mutate(met_max=ifelse(min(logFC)*(-1)<max(logFC),(1),(-1))*max(abs(logFC)),min_p_met=min(adj.P.Val))
nrow(pgdata_met)

pgdata_met['mat_FDR_cat']<-ifelse(pgdata_met$min_p_met<0.1,"<0.1",">=0.1")
ggplot(pgdata_met,aes(met_max,-log10(min_p_met)),color='black') +
  geom_point()+
  geom_point(pgdata_met[pgdata_met$min_p_met<0.1&pgdata_met$met_max>0,],mapping=aes(met_max,-log10(min_p_met)),color="red")+
  geom_point(pgdata_met[pgdata_met$min_p_met<0.1&pgdata_met$met_max<0,],mapping=aes(met_max,-log10(min_p_met)),color="blue")

pgdata_met_extream<-pgdata_met[!duplicated(pgdata_met$SYMBOL),]


merge_extream<-merge(pgdata_met_extream,pgdata_atac_extream,by="SYMBOL",how="inner")

ggplot(merge_extream[merge_extream$FDR<2,],aes(x=atac_max,y=met_max))+
  geom_point()+xlab("atac seq log2FC")+ylab("methylation log2FC")

# List of items
x<-list("ATAC seq"=atac_sagn$SYMBOL,"methyltion data"=met_sagn$SYMBOL)
# 2D Venn diagram
ggVennDiagram(x)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

methyl_sagn_from_merge_extreame_dataset<-merge_extream[merge_extream$min_p_met<0.1,]
atac_sagn_from_merge_extream_dataset<-merge_extream[merge_extream$min_p_atac<0.1,]

overlap_genes<-merge(methyl_sagn_from_merge_extreame_dataset,atac_sagn_from_merge_extream_dataset,by="SYMBOL",how="inner")
write.csv(overlap_genes,"G:/My Drive/PhD/Uri_Data/atac_methyl_sagn_genes.csv")

x<-list("methyltion data"=methyl_sagn_from_merge_extreame_dataset$SYMBOL,"ATAC-seq data"=atac_sagn_from_merge_extream_dataset$SYMBOL)
# 2D Venn diagram
venn_diagram<-ggVennDiagram(x,label = "percent")+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

# Print plots to a pdf file
pdf("G:/My Drive/PhD/Uri_Data/plots/extream_atac_mathyl_sagn_overllap_ven_diagram.pdf")
print(venn_diagram)     # Plot 1 --> in the first page of PDF
dev.off() 



#Permutation test for getting 47 sagnificant overlap genes
#we have 490 sagnificant genes in atac, if we take random group mrom the extream mathylation data, is 74 overlapp genes is high or low overllap. 

permutation.test <- function(met_data, met_sagn,ATAC_genes, n){
  distribution=c()
  result=0
  
  for(i in 1:n){
    met_sample_index <- sample(1:nrow(met_data),nrow(met_sagn))
    #print(RNA_sample_index[1:5])
    #print(length(RNA_sample_index))
    met_sapmle_genes<-met_data[met_sample_index,]
    # print(nrow(RNA_sapmle_genes))
    #print(RNA_sapmle_genes[1:2,])
    
    
    distribution[i]<-nrow(merge(ATAC_genes,met_sapmle_genes,by="SYMBOL",how="inner"))
    
  }
  result=sum(distribution >= nrow(overlap_genes))/n
  return(list(result, distribution))
}



test1 <- permutation.test(merge_extream,methyl_sagn_from_merge_extreame_dataset, atac_sagn_from_merge_extream_dataset, 1000)
hist(test1[[2]], breaks=30, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=nrow(overlap_genes), lwd=3, col="red")
test1[[1]]

pdf("G:/My Drive/PhD/Uri_Data/plots/permutation_test_for_overlap_sagnificant_genes.pdf")
hist(test1[[2]], breaks=30, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=nrow(overlap_genes), lwd=3, col="red")
dev.off() 
###############################################################################################
################################################################################################
###############################################################################################
########################################heatmap clasification################################
################################################################################################



methyl_samples<-read.csv("G:/My Drive/PhD/Uri_Data/met_integration.csv")
atac_samples<-read.csv("G:/My Drive/PhD/Uri_Data/peaks_noanno.csv")
names(methyl_samples)
methyl_samples<-methyl_samples[,c("HC.2","HC.3","GSD1A.6","GSD1A.7","X")]
names(methyl_samples)<-c("HC1","HC2","GSD1a1","GSD1a2","Name")

#take only the sagnificant rows from atac
#methyl_samples_0.25<-methyl_samples[methyl_samples$Name %in% metil_data_names_sagn,]
#methyl_samples_0.25["data"]<-"blue"

###################normelize the methyl_samples
methyl_samples['var']<-0

for (i in 1:ncol(methyl_samples[c("HC1","HC2","GSD1a1","GSD1a2")]))
{
  print(sum(methyl_samples[,i]))
  methyl_samples[,i]<-(methyl_samples[,i]/sum(methyl_samples[,i]))*((nrow(methyl_samples)+nrow(atac_samples))/nrow(atac_samples))
  print(sum(methyl_samples[,i]))
}

for (i in 1:nrow(methyl_samples))
{
  #print(as.numeric(methyl_samples[i,c("HC1","HC2","GSD1a1","GSD1a2")]))
  methyl_samples[i,'var']<-var(as.numeric(methyl_samples[i,c("HC1","HC2","GSD1a1","GSD1a2")]))
}

hist(log10(methyl_samples$var),breaks = 100)
sum(methyl_samples$var>0.000000000027)

high_var_methyl<-methyl_samples[methyl_samples$var>0.000000000027,]

nrow(methyl_samples)
methyl_random_indexs<-sample(1:nrow(methyl_samples),10000, replace = FALSE)
methyl_samples_rand<-methyl_samples[methyl_random_indexs,]

names(atac_samples$GSD1a2)
atac_samples["X"]<-paste0(substr(atac_samples$CHR, start = 8, stop = 10),atac_samples$START)
atac_samples<-atac_samples[c("HC1","HC2","GSD1a1","GSD1a2")]
##atac_samples["data"]<-"red"

atac_samples['var']<-0

###################normelize the atac_samples
for (i in 1:ncol(atac_samples[c("HC1","HC2","GSD1a1","GSD1a2")]))
{
  print(sum(atac_samples[,i]))
  atac_samples[,i]<-(atac_samples[,i]/sum(atac_samples[,i]))*((nrow(methyl_samples)+nrow(atac_samples))/nrow(methyl_samples))
  print(sum(atac_samples[,i]))
}
nrow(atac_samples)

for (i in 1:nrow(atac_samples))
{
  print(as.numeric(methyl_samples[i,c("HC1","HC2","GSD1a1","GSD1a2")]))
  print(var(as.numeric(methyl_samples[i,c("HC1","HC2","GSD1a1","GSD1a2")])))
  atac_samples[i,'var']<-var(as.numeric(methyl_samples[i,c("HC1","HC2","GSD1a1","GSD1a2")]))
  print(atac_samples[i,'var'])
}
hist(log10(atac_samples$var),breaks = 100)
sum(atac_samples$var>0.0000000000085)

high_var_atac<-atac_samples[atac_samples$var>0.0000000000085,]

atac_random_indexs<-sample(1:nrow(atac_samples),10000, replace = FALSE)
atac_rand_samples_data<-atac_samples[atac_random_indexs,]

##########

samples_data<-rbind(high_var_methyl[c("HC1","HC2","GSD1a1","GSD1a2")],high_var_atac[c("HC1","HC2","GSD1a1","GSD1a2")])
nrow(samples_data)


data<-as.matrix(samples_data[c("HC1","HC2","GSD1a1","GSD1a2")])



pdf("G:/My Drive/PhD/Uri_Data/plots/heatmap_classification.pdf")
heatmap(data)
dev.off() 


### GSEA
integ_data = read.csv('D:/MiguelW12/Documents/integ_do.csv', header=TRUE)
names(integ_data)
gene_list_entrez = integ_data$SYMBOL
original_gene_list_ENTREZ <- integ_data$Fold
eg = bitr(gene_list_entrez, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(original_gene_list_ENTREZ) <- eg$ENTREZID
original_gene_list_ENTREZ = sort(original_gene_list_ENTREZ, decreasing = TRUE)

###run
y <- gseDO(original_gene_list_ENTREZ,
            minGSSize     = 1,
            pvalueCutoff  = 1,
            pAdjustMethod = "BH",
            verbose       = FALSE)

write.csv(as.data.frame(y), file = "GSEA-DO-integ.csv")
#subset
do_id_list <- c("DOID:4971", "DOID:3369", "DOID:201", "DOID:26", "DOID:0060100")
y@result = y@result[y@result$ID %in% do_id_list, ]
do_mean <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
#vis
dotplot(do_mean, showCategory=5, font.size = 8,label_format = 60, color = "p.adjust", title = "dotplot_integration", split=".sign") + facet_grid(.~.sign)

######################################################