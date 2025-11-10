#load libraries
libs <- c("DESeq2","BiocParallel","ggplot2","tximport","sva","ggvenn","ggVennDiagram",
          "ggpubr","reshape2","car","ggtext","lme4","WGCNA","clusterProfiler","enrichplot",
          "SummarizedExperiment","biomaRt","lmerTest","GenomicRanges","org.Mm.eg.db","rtracklayer")

library(dplyr)
library(ggplot2)
library(ggvenn)
library(rtracklayer)
library(GenomicRanges)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggpubr)
library("reshape2")
library(lme4)


lapply(libs, library, character.only=T)
#setting parallel cores

register(MulticoreParam(10))

#setting working directory
setwd("~/Documents/phd/comp/proj2/23042025/")

#####
#initial analysis with all but natural borrelia infections and tlr2 genotype 4
#read metadata file
metaData <- read.table("metadata.tsv", header=T)

#removing genotype c3/c3 and natural infections of borrelia
metaData <- subset(metaData, (Genotype!=4 & Nat.Inf.!=1))

metaData <- metaData[order(metaData$NGI.ID),]

#converting to factor
metaData <- as.data.frame(lapply(metaData, as.factor))

#creating extra column in metaData for combos of condition+ gentoype + sex
metaData$group <- as.factor(paste(metaData$Condition, metaData$Genotype, sep="."))
metaData$genoSex <- as.factor(paste(metaData$Genotype, metaData$Sex, sep="."))
metaData$groupSex <- as.factor(paste(metaData$group, metaData$Sex, sep="."))
metaData$conditionSex <- as.factor(paste(metaData$Condition, metaData$Sex, sep="."))



#coding genotype as a linear variable
metaData$Genotype.linear <- c(-1, 0, 1)[match(metaData$Genotype, c(1, 2, 3))]


#checking strcuture metaData
str(metaData)


# ####Reading in count estimates from rsem and batch correcting########
# #creating file names
# files <- file.path("~/Documents/phd/comp/proj2/rsem_res", paste(metaData$NGI.ID,"genes.results",sep="."))
# 
# names(files) <- metaData$NGI.ID
# 
# #importing counts
# txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
#                 txIdCol = "transcript_id", abundanceCol = "TPM",
#                 countsCol = "expected_count", lengthCol = "effective_length")
# 
# #setting 0 read count to 1
# txi$length[txi$length == 0] <- 1
# 
# 
# #importing protein gene IDs
# prot_id <- read.table("../protein_ids.txt", header = F)
# 
# #retaining only protein coding genes
# txi <- lapply(c(1:3), function(x){
#   txi[[x]][rownames(txi[[x]]) %in% prot_id$V1,]})
# 
# names(txi) <- c("abundance","counts","length")
# txi$countsFromAbundance <- "no"
# 
# #import toga-ncbi matched annotations to improve annotations
# annot <- read.table("../bv_chromium.ncbi_toga.cds.intersected.summary.filt3.txt", header = F)
# 
# loc.annot <- annot[grep("LOC", annot$V1), ]
# 
# matched_indices <- match(rownames(txi$abundance), loc.annot$V1)
# 
# # Replace row names where matches are found
# rownames(txi$abundance)[!is.na(matched_indices)] <- loc.annot$V11[matched_indices[!is.na(matched_indices)]]
# rownames(txi$counts)[!is.na(matched_indices)] <- loc.annot$V11[matched_indices[!is.na(matched_indices)]]
# rownames(txi$length)[!is.na(matched_indices)] <- loc.annot$V11[matched_indices[!is.na(matched_indices)]]
# 
# txi_cp <- txi
# 
# #removing low gene counts
# keep <- rowMedians(txi_cp$abundance) >= 0.1
# 
# txi_cp <- lapply(c(1:3), function(x){
#   txi_cp[[x]] <- txi_cp[[x]][keep,]})
# 
# names(txi_cp) <- c("abundance","counts","length")
# txi_cp$countsFromAbundance <- "no"
# 
# covar_mat <- cbind(cbind(metaData$Condition, metaData$Sex), metaData$Genotype)
# 
# #batch correction with combatseq
# txi_cp$counts <- ComBat_seq(as.matrix(txi_cp$counts), batch=metaData$DOE, covar_mod = covar_mat)
# 
# #exporting batch corrected estimates
# write.table(cbind(gene = rownames(txi_cp$counts), txi_cp$counts),"counts.rsem.tsv",sep="\t", quote = F, row.names = F, col.names = T)
# write.table(cbind(gene = rownames(txi_cp$abundance), txi_cp$abundance),"abundance.rsem.tsv",sep="\t", quote = F, row.names = F, col.names = T)
# write.table(cbind(gene = rownames(txi_cp$length), txi_cp$length),"length.rsem.tsv",sep="\t", quote = F, row.names = F, col.names = T)

#making deseqobject

set.seed(42)

#importing rsem counts
files <- file.path("~/Documents/phd/comp/proj2/23042025/", paste(c("abundance","counts","length"),"rsem","tsv",sep="."))

txi_cp <- lapply(files, function(x){
  df <- read.table(x, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  mat <- as.matrix(df[ , -1])
  rownames(mat) <- df$gene
  return(mat)
  })
names(txi_cp) <- c("abundance","counts","length")
txi_cp$countsFromAbundance <- "no"

#creating deseq object
dds <- DESeqDataSetFromTximport(txi_cp, colData = metaData, design = ~1)

#setting control as the reference level for further analysis
dds$Condition <- relevel(dds$Condition, ref = c("Unstimulated"))

#transforming the normalised read counts using regularised log
vst.dds <- vst(dds, blind=TRUE)



#PCA plot with Infection as the condition of interest; plots top 500 genees by variance
pca <- plotPCA(vst.dds, intgroup=c("Condition","Genotype","Sex"), returnData=T) 
pca$Condition <- factor(pca$Condition, levels=c("Borrelia","GAS","Unstimulated"))
percentVar <- round(100 * attr(pca, "percentVar"))

pdf("pca_batchCorrected_allSamples.pdf", width=10, height=7)
ggplot(pca, aes(-PC1, PC2, color=Condition, shape=Genotype, alpha=Sex)) +
geom_point(inherit.aes = T,aes(shape=Genotype), size=4) + xlab(paste("PC1: ",percentVar[1],"% variance")) +
ylab(paste("PC2: ",percentVar[2],"% variance")) + scale_alpha_manual(values = c(0.5, 1)) +
coord_fixed() + scale_shape_manual(labels=c("c1/c1","c1/c2","c2/c2"), values=c(16,17,15)) + scale_color_manual(labels = c("*B. afzelii*", "*S. pyogenes*","Unstimulated"), values=c("#78A71F","#71a3d7","#FFB419")) + theme_linedraw() +
theme(axis.text=element_text(size=18), axis.title=element_text(size=20),  legend.text = element_markdown(size=17), legend.background = element_blank(), legend.position = "right",legend.title = element_text(size=18),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 1, color="black")) +
guides(color=guide_legend(override.aes = list(shape = 18, size=5, linetype=NA)), alpha=guide_legend(override.aes = list(shape = 18, size=5, linetype=NA)))
dev.off()

##########

##DGE BY CONDITION##

#remaking deseq object
dds.con <- DESeqDataSetFromTximport(txi_cp, colData = metaData, design =~ Condition)

dds.con$Condition <- relevel(dds.con$Condition, ref = c("Unstimulated"))

#calculating de
dds.con.de <- DESeq(dds.con, parallel = T)

#fetching contrast level
resultsNames(dds.con.de)

#degs b/w borrelia and unstimulated
res.bor <- results(dds.con.de, alpha=0.05, name="Condition_Borrelia_vs_Unstimulated")
res.bor <- na.omit(res.bor)
res.bor.de <- res.bor[with(res.bor, padj< 0.05 & abs(log2FoldChange) >=1),]
summary(res.bor.de)



#degs b/w gas and unstimulated
res.gas <- results(dds.con.de, alpha=0.05, name="Condition_GAS_vs_Unstimulated")
res.gas <- na.omit(res.gas)
res.gas.de <- res.gas[with(res.gas, padj< 0.05 & abs(log2FoldChange) >=1),]
summary(res.gas.de)


#venn digaram
bor.gas.up <- list("B. afzelii"=rownames(res.bor.de[with(res.bor.de, log2FoldChange>=1),]), "S. pyogenes"=rownames(res.gas.de[with(res.gas.de, log2FoldChange>=1),]))
venn.bor.gas.up <- ggvenn(bor.gas.up,fill_color = c("#78A71F","#71a3d7"), text_size = 6, set_name_size = 0, fill_alpha=0.7)
venn.bor.gas.up <- venn.bor.gas.up + theme( title = element_text(size=18), plot.title = element_text(hjust = 0.5, vjust=-2.9),plot.tag.position = c(0.05, 0.88), plot.tag = element_markdown(size=16))+ labs(title="Up-regulated genes")


bor.gas.down<- list("B. afzelii"=rownames(res.bor.de[with(res.bor.de, log2FoldChange<=-1),]), "S. pyogenes"=rownames(res.gas.de[with(res.gas.de,log2FoldChange<=-1),]))
venn.bor.gas.down <- ggvenn(bor.gas.down, fill_color = c("#78A71F","#71a3d7"), text_size=6, set_name_size = 0, fill_alpha=0.7)
venn.bor.gas.down <- venn.bor.gas.down + annotate("text", x = c(-0.7,0.7), y = c(-1.1,-1.1), label = c(expression(italic("B. afzelii")), expression(italic("S. pyogenes"))), size=6)+  theme(title = element_text(size=18), plot.title = element_text(hjust = 0.5, vjust=-2.9),plot.tag.position = c(0.05, 0.88))+ labs(title="Down-regulated genes") 


##GO TERM analysis
#list of degs
bor.gas.degs <- list("B. afzelii"=rownames(res.bor.de), "S. pyogenes"=rownames(res.gas.de))

#comparing common go terms
bor.gas.degs.go <- compareCluster(bor.gas.degs, universe=rownames(dds.con) ,fun = enrichGO, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff = 0.05, pAdjustMethod="fdr", keyType = "SYMBOL")

#finding pairwise terms
bor.gas.degs.go.plot <- dotplot(bor.gas.degs.go, showCategory=10, font.size=15)

#saving plot
pdf("degs_condition_venn_GO.pdf", height = 8, width=13)
ggarrange(ggarrange(venn.bor.gas.up, venn.bor.gas.down, nrow =2,common.legend = T),  bor.gas.degs.go.plot, ncol=2, widths=c(0.8, 1), labels = "AUTO", font.label = list(size=20))
dev.off()



#linear modelling for condition borrelia
#subsetiing degs across borrelia infections
vst.bor.cond.degs <- vst.dds[rownames(vst.dds) %in% bor.gas.degs$`B. afzelii` ,colnames(vst.dds) %in% metaData$NGI.ID[metaData$Bor.paired==1]]

#calcualting prinicipal compoennts for linear modeling
rv.pca.bor.cond <- prcomp(t(assay(vst.bor.cond.degs)))

#creating df for the loadings
pca.bor.cond.loadings <- as.data.frame(cbind( "NGI.ID"= rownames(rv.pca.bor.cond$x),rv.pca.bor.cond$x))

#merging laoding with metaData
pca.bor.cond.loadings <- merge(pca.bor.cond.loadings, metaData, by="NGI.ID")

#converting loading to numeric for linear modeling later
pca.bor.cond.loadings$PC1 <- as.numeric(pca.bor.cond.loadings$PC1)

#pc1 boxplot
box.pc1.bor.cond.degs.all <-  ggplot(pca.bor.cond.loadings) +geom_boxplot(aes(x=group, y=-PC1, group=group, fill=Condition))+
  scale_fill_manual(labels=c("*B. afzelii*", "Unstimulated"),values=c("#78A71F","#FFB419")) + geom_point(aes(x=group, y=-PC1, group=groupSex, shape=Sex), position = position_jitterdodge(), size=1.5, alpha=1/2,colour = "black", fill = "white", stroke = 1.25) +
  scale_shape_manual(labels=c("F","M"), values=c(21,23)) + theme_bw() + scale_x_discrete(labels=c("c1/c1", "c1/c2","c2/c2","c1/c1", "c1/c2","c2/c2")) +
  theme (axis.text = element_text(size=18), axis.title = element_text(size=20), legend.text = element_markdown(size=17), legend.title = element_text(size=18),legend.background = element_blank(), legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 1, color="black")) +
  labs(x="Genotype", y="PC1") + guides(shape=guide_legend(override.aes = list(size=3)))


pdf("degs_bor_cond_boxAll.pdf", height = 6, width=8)
box.pc1.bor.cond.degs.all
dev.off()


#checking for association between genotype, condition and pc1 and progressively eliminating non-significant factors 

pc1.bor.cond.lmer <- lmer(PC1 ~ Condition  + Genotype.linear + Sex + Genotype.linear:Condition+ Condition:Sex+ Genotype.linear:Sex + Condition:Sex:Genotype.linear+ (1|Vole), data=pca.bor.cond.loadings)
car::Anova(pc1.bor.cond.lmer, type="3", test.statistic = "F")


pc1.bor.cond.lmer <- lmer(PC1 ~ Condition  + Genotype.linear + Sex + Genotype.linear:Condition+ Condition:Sex+ Genotype.linear:Sex + (1|Vole), data=pca.bor.cond.loadings)
car::Anova(pc1.bor.cond.lmer, type="3", test.statistic = "F")

pc1.bor.cond.lmer <- lmer(PC1 ~ Condition  + Genotype.linear + Sex + Genotype.linear:Condition+ Condition:Sex+ (1|Vole), data=pca.bor.cond.loadings)
car::Anova(pc1.bor.cond.lmer, type="3", test.statistic = "F")
ranova(pc1.bor.cond.lmer)



#linear modelling for condition gas
#subsetiing degs common to GAS infections
vst.gas.cond.degs <- vst.dds[rownames(vst.dds) %in% bor.gas.degs$`S. pyogenes` ,colnames(vst.dds) %in% metaData$NGI.ID[metaData$GAS.paired==1]]


#calculating principal components
rv.pca.gas.cond <- prcomp(t(assay(vst.gas.cond.degs)))

#fetching loading
pca.gas.cond.loadings <- as.data.frame(cbind( "NGI.ID"= rownames(rv.pca.gas.cond$x),rv.pca.gas.cond$x))

#merging with metadata
pca.gas.cond.loadings <- merge(pca.gas.cond.loadings, metaData, by="NGI.ID")

#converting to numeric
pca.gas.cond.loadings$PC1 <- as.numeric(pca.gas.cond.loadings$PC1)

#pc1 boxplot
box.pc1.gas.cond.degs.all <-  ggplot(pca.gas.cond.loadings) +geom_boxplot(aes(x=group, y=-PC1, group=group, fill=Condition))+
  scale_fill_manual(labels=c("*S.pyogenes*", "Unstimulated"),values=c("#71a3d7","#FFB419")) + geom_point(aes(x=group, y=-PC1, group=groupSex, shape=Sex), position = position_jitterdodge(), size=1.5, alpha=1/2,colour = "black", fill = "white", stroke = 1.25) +
  scale_shape_manual(labels=c("F","M"), values=c(21,23)) + theme_bw() + scale_x_discrete(labels=c("c1/c1", "c1/c2","c2/c2","c1/c1", "c1/c2","c2/c2")) +
  theme (axis.text = element_text(size=18), axis.title = element_text(size=20), plot.tag = element_markdown(size=18), legend.text = element_markdown(size=18), legend.title = element_text(size=18),legend.background = element_blank(), legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 1, color="black")) +
  labs(x="Genotype", y="PC1") + guides(shape=guide_legend(override.aes = list(size=3)))

pdf("degs_gas.cond_boxAll.pdf", height = 6, width=8)
box.pc1.gas.cond.degs.all
dev.off()

#checking for association between genotype, condition and pc1 and progressively removing non-significant terms 
pc1.gas.cond.lmer <- lmer(PC1 ~ Condition  + Genotype.linear + Sex + Genotype.linear:Condition+ Condition:Sex+ Genotype.linear:Sex + Condition:Sex:Genotype.linear+ (1|Vole), data=pca.gas.cond.loadings)
car::Anova(pc1.gas.cond.lmer, type="3", test.statistic = "F")

pc1.gas.cond.lmer <- lmer(PC1 ~ Condition  + Genotype.linear + Sex + Genotype.linear:Condition+ Condition:Sex+ Genotype.linear:Sex + (1|Vole), data=pca.gas.cond.loadings)
car::Anova(pc1.gas.cond.lmer, type="3", test.statistic = "F")

pc1.gas.cond.lmer <- lmer(PC1 ~ Condition  + Genotype.linear + Sex + Genotype.linear:Condition+ Condition:Sex+ (1|Vole), data=pca.gas.cond.loadings)
car::Anova(pc1.gas.cond.lmer, type="3", test.statistic = "F")

pc1.gas.cond.lmer <- lmer(PC1 ~ Condition  + Genotype.linear + Sex + Genotype.linear:Condition+ (1|Vole), data=pca.gas.cond.loadings)
car::Anova(pc1.gas.cond.lmer, type="3", test.statistic = "F")

pc1.gas.cond.lmer <- lmer(PC1 ~ Condition  + Genotype.linear + Sex +(1|Vole), data=pca.gas.cond.loadings)
car::Anova(pc1.gas.cond.lmer, type="3", test.statistic = "F")
ranova(pc1.gas.cond.lmer)

###########################################
##DGE by genotpye in each condition########
condition <- c("Borrelia", "GAS", "Unstimulated")

dds.geno <- DESeqDataSetFromTximport(txi_cp, colData = metaData, design = ~ Genotype)

geno.de.list <- list()
res.geno.de.list <- list()


for(i in 1:length(condition))
{
  #object for unstimulated, borrelia, gas separatel
  de.geno <- dds.geno[,colnames(dds.geno) %in% metaData$NGI.ID[metaData$Condition==condition[i]]]
  
  #relevel by genotype 1
  de.geno$Genotype <- relevel(de.geno$Genotype, ref = c("1"))
  
  #calculating deseq
  de.geno.de <- DESeq(de.geno, parallel = T)
  resultsNames(de.geno.de)
  
  #geno2 vs 1 results
  res.2VS1 <- results(de.geno.de,alpha=0.05, name="Genotype_2_vs_1")
  res.2VS1 <- res.2VS1[!(is.na(res.2VS1$padj)),]
  res.2VS1.de <- res.2VS1[with(res.2VS1, padj< 0.05 & abs(log2FoldChange) >=1),] 
  
  #geno3 vs 1 results
  res.3VS1 <- results(de.geno.de,alpha=0.05, name="Genotype_3_vs_1")
  res.3VS1 <- res.3VS1[!(is.na(res.3VS1$padj)),]
  res.3VS1.de <- res.3VS1[with(res.3VS1, padj< 0.05 & abs(log2FoldChange) >=1),] 
  
  #geno3 vs 2 results
  res.3VS2 <- results(de.geno.de,alpha=0.05, contrast = c("Genotype", "3","2"))
  res.3VS2 <- res.3VS2[!(is.na(res.3VS2$padj)),]
  res.3VS2.de <- res.3VS1[with(res.3VS2, padj< 0.05 & abs(log2FoldChange) >=1),] 
  
  #writing de onject
  geno.de.list[[condition[i]]] <- de.geno.de
  
  #writing results
  res.geno.de.list[[paste(condition[i],"2VS1", sep=".")]] <- res.2VS1.de
  res.geno.de.list[[paste(condition[i],"3VS1", sep=".")]] <- res.3VS1.de
  res.geno.de.list[[paste(condition[i],"3VS2", sep=".")]] <- res.3VS2.de
}


##number of genes for each condition
lapply(res.geno.de.list , nrow)

##GO TERM analysis
#list of degs
geno.degs <- list("B. afzelii c1/c2 vs c1/c1"=rownames(res.geno.de.list$Borrelia.2VS1), "B. afzelii c2/c2 vs c1/c1"=rownames(res.geno.de.list$Borrelia.3VS1),"B. afzelii c2/c2 vs c1/c2"=rownames(res.geno.de.list$Borrelia.3VS2),"S. pyogenes c2/c2 vs c1/c1"=rownames(res.geno.de.list$GAS.3VS1), "S. pyogenes c2/c2 vs c1/c2"=rownames(res.geno.de.list$GAS.3VS2),"S. pyogenes c1/c2 vs c1/c1"=rownames(res.geno.de.list$GAS.2VS1), "Unstimulated c2/c2 vs c1/c1"=rownames(res.geno.de.list$Unstimulated.3VS1), "Unstimulated c2/c2 vs c1/c2"=rownames(res.geno.de.list$Unstimulated.3VS2), "Unstimulated c1/c2 vs c1/c1"=rownames(res.geno.de.list$Unstimulated.2VS1))

#comparing common go terms
geno.degs.go <- compareCluster(geno.degs, universe=rownames(dds) ,fun = enrichGO, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff = 0.05, pAdjustMethod="fdr", keyType = "SYMBOL")

#finding pairwise terms
geno.degs.go <- pairwise_termsim(geno.degs.go)


pdf("genotype_GO.pdf", height = 12, width=15)
dotplot(geno.degs.go, showCategory=5,font.size=15)
dev.off()

########################################
#########   WGCNA   ####################
########################################

allowWGCNAThreads()
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#wgcna requires genes as columns and rows as samples;transforming
datExpr0 <- as.data.frame(t(assay(vst.dds)))

#importing trait data
traits <-metaData

#checking for outlier genes; checks if any samples or genes are outliers in expression
gsg <- goodSamplesGenes(datExpr0, verbose = 3)

#inspecting gsg; goodGenes is logical for genes and goodSamples logical for samples
str(gsg)

#are all genes and samples ok? if FALSE, the need to be removed
gsg$allOK


##step by step module identification
# Choose a set of soft-thresholding powers
powers = c(c(1:15), seq(from = 16, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType = "signed", corFnc = WGCNA::cor,verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(2,1))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Co-expression similarity and adjacency

softPower <- 10

#splitting into blocks

bwCor <- blockwiseModules(datExpr0, power = softPower, networkType = "signed",
                          deepSplit = 2,
                          minModuleSize = 30,
                          maxBlockSize = 20000,
                          mergeCutHeight = 0.25,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          numericLabels = T,
                          verbose = 3,
                          corType = "bicor",
                          maxPOutliers = 0.05)


#saving Eigene Genes as separate object; eigen genes are representative values for each module (==PC1)
bwEigenGenes <- bwCor$MEs

#changing column names
colnames(bwEigenGenes) <- sub("ME","Module.", colnames(bwEigenGenes))

#no. of samples
nSamples <- nrow(datExpr0)

#no. of genes
nGenes <- ncol(datExpr0)

#correlating eigen genes with expression data to obtain genes with high module membership
moduleGeneMembership <- WGCNA::cor(datExpr0, bwEigenGenes, use='p')

moduleGeneMembershipPvalue <- corPvalueStudent(moduleGeneMembership, nSamples) 
moduleGeneMembershipPvalue.adj <- apply(moduleGeneMembershipPvalue, 2, function(x){p.adjust(x, method="fdr")})



#setting levels for the traits
traits$Condition <- relevel(traits$Condition, ref=c("Unstimulated"))

#converting categorical variable into binaries; binarise compares factors against the refernce level, so creating separate column for the reference
#condition
cond.bin <- binarizeCategoricalColumns(traits$Condition)
cond.bin$data.Unstimulated.vs.all <- abs(rowSums(cond.bin[,1:2])-1)


trait.bin <- cond.bin
colnames(trait.bin) <- c("B. afzelii", "S. pyogenes", "Unstimulated")

#setign rownames to user id
rownames(trait.bin) <- traits$NGI.ID



#calcualting correlation between eigen genes and traits
moduleTraitCor <- WGCNA::cor(bwEigenGenes, trait.bin, use='p')

#estimating correlations between modules and trai and pvalues for correlation
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

moduleTraitPvalue.adj <- apply(moduleTraitPvalue, 2, function(x){p.adjust(x, method="fdr")})


#creating heatmap data
heatmap.data <- cbind(merge(bwEigenGenes,trait.bin,by='row.names'), Tlr2=datExpr0[,c("Tlr2")])

heatmap.data <- heatmap.data[,c("Row.names",paste("Module", 0:(ncol(bwEigenGenes)-1), sep="."), colnames(trait.bin),"Tlr2")]  

moduleTraitCor.sorted <- moduleTraitCor[paste("Module", 0:(ncol(bwEigenGenes)-1), sep="."),]

#adding association with Tlr2 expression
moduleTraitCor.sorted <- cbind(moduleTraitCor.sorted, Tlr2=moduleGeneMembership[rownames(moduleGeneMembership)=="Tlr2",paste("Module", 0:(ncol(bwEigenGenes)-1),sep=".")])

moduleTraitCorPval.sorted <- moduleTraitPvalue.adj[paste("Module", 0:(ncol(bwEigenGenes)-1),sep="."), ]

moduleTraitCorPval.sorted <- cbind(moduleTraitCorPval.sorted, Tlr2=moduleGeneMembershipPvalue.adj[rownames(moduleGeneMembershipPvalue.adj)=="Tlr2",paste("Module", 0:(ncol(bwEigenGenes)-1),sep=".")])


# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor.sorted, 2), "\n(",
                   signif(moduleTraitCorPval.sorted, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor.sorted)


#mapping gens with module contributions
moduleGeneMap <- data.frame(
  gene_id = names(bwCor$colors),
  colors = labels2colors(bwCor$colors),
  modules= bwCor$colors
)


#creating labels for plotting
module.table <- as.data.frame(table(paste("Module",moduleGeneMap$modules, sep=".")))

module.table <- module.table[match(grep("Module",colnames(heatmap.data), value = T), module.table$Var1),]

colnames(heatmap.data[,grep("Module",colnames(heatmap.data), value = T)]) <- paste(grep("Module",colnames(heatmap.data), value = T), module.table$Freq,sep="\n")

colnames(heatmap.data)[2:(ncol(bwEigenGenes)+1)] <- paste(grep("Module",colnames(heatmap.data), value = T), "\n(", module.table$Freq, ")",sep="")



### conditon + tlr2 heatmap
pdf("wgcna_moduleTrait_heatmap.pdf", height = 6, width=26)

par(mar = c(10, 10, 3, 3))
labeledHeatmap(Matrix = t(moduleTraitCor.sorted),
              yLabels = names(heatmap.data)[(ncol(bwEigenGenes)+2):ncol(heatmap.data)],
             xLabels = names(heatmap.data)[2:(ncol(bwEigenGenes)+1)],
            xSymbols = names(heatmap.data)[2:(ncol(bwEigenGenes)+1)],
           colorLabels = F,
          colors = blueWhiteRed(10),
         textMatrix = t(textMatrix),
        setStdMargins = FALSE,
       cex.text = 1.25,
      cex.lab = 1.5,
     zlim = c(-1,1),
    cex.main=2,
   main = paste("Module-trait relationship"))

dev.off()


##plot pc1 for all modules for all conditions#####

boxPlot.allModule <- lapply(0:(nrow(module.table)-1), function(x){
  vst.mod <- vst.dds[rownames(vst.dds) %in% moduleGeneMap$gene_id[moduleGeneMap$modules==x],]
  #calculating prinicipal component
  rv.mod <- prcomp(t(assay(vst.mod)))
  
  #making dataframe with all elements
  mod.pca <- as.data.frame(cbind( "NGI.ID"= rownames(rv.mod$x),rv.mod$x))
  
  #merging with metadata
  mod.pca <- merge(mod.pca, metaData, by="NGI.ID")
  
  #setting pc1 to numeric
  mod.pca$PC1 <- as.numeric(mod.pca$PC1)
  
  module <- paste("Module",x,sep=".")
  
  #boxplot 
  ggplot(mod.pca, aes(x=group, y=-PC1, fill=Condition, group=group)) + geom_boxplot(show.legend=F) + 
    scale_x_discrete(labels=c("c1/c1", "c1/c2","c2/c2","c1/c1", "c1/c2","c2/c2","c1/c1", "c1/c2","c2/c2")) + 
    scale_fill_manual(labels=c("*B. afzelii*","*S. pyogenes*", "Unstimulated"),values=c("#78A71F","#71a3d7","#FFB419")) + theme_bw() + 
    theme (axis.text = element_markdown(size=12), axis.title = element_markdown(size=22), legend.text = element_markdown(size=16), legend.title = element_text(size=17),legend.background = element_rect(linewidth = 0.5, color = "black"), legend.position = "right") + 
    labs(x="",y="", title = module)
  
})

names(boxPlot.allModule) <- paste("Module", 0:(nrow(module.table)-1),sep=".")

goTerm.Module <- lapply(0:(nrow(module.table)-1), function(x){
  vst.mod <- vst.dds[rownames(vst.dds) %in% moduleGeneMap$gene_id[moduleGeneMap$modules==x],]
  mod.go <- enrichGO(rownames(vst.mod), universe=names(datExpr0) ,OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff = 0.05, pAdjustMethod="fdr", keyType = "SYMBOL")
  name <- paste("Module",x,sep=".")
  if(nrow(mod.go) > 0){
    print(x)
    #mod.go.plot <- pairwise_termsim(goTerm.Module[[name]])
    dotplot(mod.go, showCategory=10, title=name)
    
  }
})

goTerm.Module.plot <- Filter(Negate(is.null), goTerm.Module)



pdf("allModule_pc1_boxplot.pdf", height = 25, width=20)
do.call(ggpubr::ggarrange, c(boxPlot.allModule, ncol = 4, nrow = 6, common.legend = TRUE))
dev.off()

pdf("allModule_go.pdf", height = 26, width=24)
do.call(ggpubr::ggarrange, c(goTerm.Module.plot, nrow = ceiling(length(goTerm.Module.plot)/ 4), ncol=4))
dev.off()

moduleGeneID <-  lapply(c("2","9","15"), function(x){moduleGeneMap$gene_id[moduleGeneMap$modules==x]})
names(moduleGeneID) <- c("Module.2","Module.9","Module.15")

##genes significantly associated with borrelia condition
geneTraitSignificance.bor <- WGCNA::cor(datExpr0, trait.bin$`B. afzelii`, use = "p")
GS.pvalue.bor <- corPvalueStudent(geneTraitSignificance.bor, nSamples)
GS.pval.adj.bor <- as.matrix(p.adjust(GS.pvalue.bor, method="fdr"))
rownames(GS.pval.adj.bor) <- rownames(GS.pvalue.bor)

#finding key genes driving pattern in module 2; gene association with module and condition > 0.7  
geneDrivers.module2.bor <- data.frame(Module.mem = moduleGeneMembership[, "Module.2"], Module.mem.padj= moduleGeneMembershipPvalue.adj[,"Module.2"],Trait.mem=geneTraitSignificance.bor, Trait.mem.padj=GS.pval.adj.bor )

key.geneDrivers.module2.bor <- geneDrivers.module2.bor[with(geneDrivers.module2.bor, rownames(geneDrivers.module2.bor) %in% moduleGeneID$Module.2 & Module.mem > 0.7 & Trait.mem > 0.55),]

#finding key genes driving pattern in module 9; gene association with module and condition > 0.7  
geneDrivers.module9.bor <- data.frame(Module.mem = moduleGeneMembership[, "Module.9"], Module.mem.padj= moduleGeneMembershipPvalue.adj[,"Module.9"],Trait.mem=geneTraitSignificance.bor, Trait.mem.padj=GS.pval.adj.bor )

key.geneDrivers.module9.bor<- geneDrivers.module9.bor[with(geneDrivers.module9.bor, rownames(geneDrivers.module9.bor) %in% moduleGeneID$Module.9 & Module.mem > 0.7 & Trait.mem > 0.7),]

#GO term for module 9 key genes
key.geneDrivers.module9.bor.go <- enrichGO(rownames(key.geneDrivers.module9.bor), universe=names(datExpr0) ,OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff = 0.05, pAdjustMethod="fdr", keyType = "SYMBOL")

key.geneDrivers.module9.bor.go.plot <- dotplot(key.geneDrivers.module9.bor.go, showCategory=10, title="Module 9 key genes", font.size=18)

pdf("keyModule9_keyGenes_borrelia.pdf", height = 8, width=8)
key.geneDrivers.module9.bor.go.plot
dev.off()

##repating for GAS; gene-trait association
geneTraitSignificance.gas <- WGCNA::cor(datExpr0, trait.bin$`S. pyogenes`, use = "p")
GS.pvalue.gas <- corPvalueStudent(geneTraitSignificance.gas, nSamples)
GS.pval.adj.gas <- as.matrix(p.adjust(GS.pvalue.gas, method="fdr"))
rownames(GS.pval.adj.gas) <- rownames(GS.pvalue.gas)

#finding key genes driving pattern in module 2; gene association with module and condition > 0.7  
geneDrivers.module2.gas <- data.frame(Module.mem = moduleGeneMembership[, "Module.2"], Module.mem.padj= moduleGeneMembershipPvalue.adj[,"Module.2"],Trait.mem=geneTraitSignificance.gas, Trait.mem.padj=GS.pval.adj.gas )
key.geneDrivers.module2.gas <- geneDrivers.module2.gas[with(geneDrivers.module2.gas, rownames(geneDrivers.module2.gas) %in% moduleGeneID$Module.2 & Module.mem > 0.7 & Trait.mem > 0.55),]

#finding key genes driving pattern in module 15; gene association with module and condition > 0.7  
geneDrivers.module15.gas <- data.frame(Module.mem = moduleGeneMembership[, "Module.15"], Module.mem.padj= moduleGeneMembershipPvalue.adj[,"Module.15"],Trait.mem=geneTraitSignificance.gas, Trait.mem.padj=GS.pval.adj.gas )
key.geneDrivers.module15.gas<- geneDrivers.module15.gas[with(geneDrivers.module15.gas, rownames(geneDrivers.module15.gas) %in% moduleGeneID$Module.15 & Module.mem > 0.7 & Trait.mem > 0.7),]


#Linear modelling of module 2 genes for borrelia, genotype as linear effect

vst.mod2.bor <- vst.dds[rownames(vst.dds) %in% moduleGeneID$Module.2, colnames(vst.dds) %in% metaData$NGI.ID[metaData$Bor.paired==1]]

#calculating prinicipal component
rv.mod2.pca.bor <- prcomp(t(assay(vst.mod2.bor)))

#making dataframe with all elements
mod2.pca.loading.bor <- as.data.frame(cbind( "NGI.ID"= rownames(rv.mod2.pca.bor$x),rv.mod2.pca.bor$x))


#merging with metadata
mod2.pca.loading.bor <- merge(mod2.pca.loading.bor, metaData, by="NGI.ID")

#setting pc1 to numeric
mod2.pca.loading.bor$PC1 <- as.numeric(mod2.pca.loading.bor$PC1)

#making condition as factor
mod2.pca.loading.bor$Condition <- as.factor(as.character(mod2.pca.loading.bor$Condition))

#linear modelling for module 2 affecting borrelia; successively removing non-sign terms
lmer.mod2.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ Sex:Condition:Genotype.linear+ (1|Vole), data=mod2.pca.loading.bor)
car::Anova(lmer.mod2.linear.pc1.bor, type="3", test.statistic="F")

lmer.mod2.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+(1|Vole), data=mod2.pca.loading.bor)
car::Anova(lmer.mod2.linear.pc1.bor, type="3", test.statistic="F")

lmer.mod2.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Condition+(1|Vole), data=mod2.pca.loading.bor)
car::Anova(lmer.mod2.linear.pc1.bor, type="3", test.statistic="F")

lmer.mod2.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Sex:Condition+(1|Vole), data=mod2.pca.loading.bor)
car::Anova(lmer.mod2.linear.pc1.bor, type="3", test.statistic="F")

lmer.mod2.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + (1|Vole), data=mod2.pca.loading.bor)
car::Anova(lmer.mod2.linear.pc1.bor, type="3", test.statistic="F")

lmer.mod2.linear.pc1.bor <- lmer(PC1 ~ Condition + Sex + (1|Vole), data=mod2.pca.loading.bor)
car::Anova(lmer.mod2.linear.pc1.bor, type="3", test.statistic="F")


#Linear model for module 9 borrelia, linear genotype
vst.mod9.bor <- vst.dds[rownames(vst.dds) %in% moduleGeneID$Module.9, colnames(vst.dds) %in% metaData$NGI.ID[metaData$Bor.paired==1]]

#calculating prinicipal component
rv.mod9.pca.bor <- prcomp(t(assay(vst.mod9.bor)))

#making dataframe with all elements
mod9.pca.loading.bor <- as.data.frame(cbind( "NGI.ID"= rownames(rv.mod9.pca.bor$x),rv.mod9.pca.bor$x))


#merging with metadata
mod9.pca.loading.bor <- merge(mod9.pca.loading.bor, metaData, by="NGI.ID")

#setting pc1 to numeric
mod9.pca.loading.bor$PC1 <- as.numeric(mod9.pca.loading.bor$PC1)


box.pc1.bor.mod9.all <- ggplot(mod9.pca.loading.bor) +geom_boxplot(aes(x=group, y=-PC1, group=group, fill=Condition))+
  scale_fill_manual(labels=c("*B. afzelii*", "Unstimulated"),values=c("#78A71F","#FFB419")) + geom_point(aes(x=group, y=-PC1, group=groupSex, shape=Sex), position = position_jitterdodge(), size=1.5, alpha=1/2,colour = "black", fill = "white", stroke = 1.25) +
  scale_shape_manual(labels=c("F","M"), values=c(21,23)) + theme_bw() + scale_x_discrete(labels=c("c1/c1", "c1/c2","c2/c2","c1/c1", "c1/c2","c2/c2")) +
  theme (axis.text = element_text(size=18), axis.title = element_text(size=20), plot.tag = element_markdown(size=16), legend.text = element_markdown(size=18), legend.title = element_text(size=18),legend.background = element_blank(), legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(linewidth = 1, color="black")) +
  labs(x="Genotype", y="PC1") + guides(shape=guide_legend(override.aes = list(size=3)))



pdf("mod9_boxplot_bor.pdf", height = 6, width=8)
box.pc1.bor.mod9.all
dev.off()


#making condition as factor
mod9.pca.loading.bor$Condition <- as.factor(as.character(mod9.pca.loading.bor$Condition))

lmer.mod9.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ Sex:Condition:Genotype.linear+ (1|Vole), data=mod9.pca.loading.bor)
car::Anova(lmer.mod9.linear.pc1.bor, type="3", test.statistic="F")

lmer.mod9.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ (1|Vole), data=mod9.pca.loading.bor)
car::Anova(lmer.mod9.linear.pc1.bor, type="3", test.statistic="F")

lmer.mod9.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Condition+ (1|Vole), data=mod9.pca.loading.bor)
car::Anova(lmer.mod9.linear.pc1.bor, type="3", test.statistic="F")
ranova(lmer.mod9.linear.pc1.bor)

#Linear model of key genes of module 9 in borrelia, linear genotype
vst.key.mod9.bor <- vst.dds[rownames(vst.dds) %in% rownames(key.geneDrivers.module9.bor), colnames(vst.dds) %in% metaData$NGI.ID[metaData$Bor.paired==1]]

#calculating prinicipal component
rv.key.mod9.pca.bor <- prcomp(t(assay(vst.key.mod9.bor)))

#making dataframe with all elements
key.mod9.pca.loading.bor <- as.data.frame(cbind( "NGI.ID"= rownames(rv.key.mod9.pca.bor$x),rv.key.mod9.pca.bor$x))


#merging with metadata
key.mod9.pca.loading.bor <- merge(key.mod9.pca.loading.bor, metaData, by="NGI.ID")

#setting pc1 to numeric
key.mod9.pca.loading.bor$PC1 <- as.numeric(key.mod9.pca.loading.bor$PC1)


#making condition as factor
key.mod9.pca.loading.bor$Condition <- as.factor(as.character(key.mod9.pca.loading.bor$Condition))

#modelling
lmer.key.mod9.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ Sex:Condition:Genotype.linear+ (1|Vole), data=key.mod9.pca.loading.bor)
car::Anova(lmer.key.mod9.linear.pc1.bor, type="3", test.statistic="F")

lmer.key.mod9.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ (1|Vole), data=key.mod9.pca.loading.bor)
car::Anova(lmer.key.mod9.linear.pc1.bor, type="3", test.statistic="F")

lmer.key.mod9.linear.pc1.bor <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Condition+ (1|Vole), data=key.mod9.pca.loading.bor)
car::Anova(lmer.key.mod9.linear.pc1.bor, type="3", test.statistic="F")


###Modelling module 2 and 9 for GAS
#Module 2 genes for GAS, genotype as linear factor
vst.mod2.gas <- vst.dds[rownames(vst.dds) %in% moduleGeneID$Module.2, colnames(vst.dds) %in% metaData$NGI.ID[metaData$GAS.paired==1]]

#calculating prinicipal component
rv.mod2.pca.gas <- prcomp(t(assay(vst.mod2.gas)))

#making dataframe with all elements
mod2.pca.loading.gas <- as.data.frame(cbind( "NGI.ID"= rownames(rv.mod2.pca.gas$x),rv.mod2.pca.gas$x))

#merging with metadata
mod2.pca.loading.gas <- merge(mod2.pca.loading.gas, metaData, by="NGI.ID")

#setting pc1 to numeric
mod2.pca.loading.gas$PC1 <- as.numeric(mod2.pca.loading.gas$PC1)

#making conditoin as factor
mod2.pca.loading.gas$Condition <- as.factor(as.character(mod2.pca.loading.gas$Condition))

lmer.mod2.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ Sex:Condition:Genotype.linear+  (1|Vole), data=mod2.pca.loading.gas)
car::Anova(lmer.mod2.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod2.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ (1|Vole), data=mod2.pca.loading.gas)
car::Anova(lmer.mod2.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod2.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Condition+ (1|Vole), data=mod2.pca.loading.gas)
car::Anova(lmer.mod2.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod2.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + (1|Vole), data=mod2.pca.loading.gas)
car::Anova(lmer.mod2.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod2.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + (1|Vole), data=mod2.pca.loading.gas)
car::Anova(lmer.mod2.linear.pc1.gas, type="3", test.statistic = "F")

#Module 9 genes for GAS, genotype as linear factor
###linear modeling of module 9 genes for gas
vst.mod9.gas <- vst.dds[rownames(vst.dds) %in% moduleGeneID$Module.9, colnames(vst.dds) %in% metaData$NGI.ID[metaData$GAS.paired==1]]

#calculating prinicipal component
rv.mod9.pca.gas <- prcomp(t(assay(vst.mod9.gas)))

#making dataframe with all elements
mod9.pca.loading.gas <- as.data.frame(cbind( "NGI.ID"= rownames(rv.mod9.pca.gas$x),rv.mod9.pca.gas$x))

#merging with metadata
mod9.pca.loading.gas <- merge(mod9.pca.loading.gas, metaData, by="NGI.ID")

#setting pc1 to numeric
mod9.pca.loading.gas$PC1 <- as.numeric(mod9.pca.loading.gas$PC1)

#plotting
# box.pc1.gas.mod9.all <- ggplot(mod9.pca.loading.gas) +geom_boxplot(aes(x=group, y=-PC1, group=group, fill=Condition))+
#   scale_fill_manual(labels=c("*S. pyogenes*", "Unstimulated"),values=c("#71a3d7","#FFB419")) + geom_point(aes(x=group, y=-PC1, group=groupSex, shape=Sex), position = position_jitterdodge(), size=1.5, alpha=1/2,colour = "black", fill = "white", stroke = 1.25) +
#   scale_shape_manual(labels=c("F","M"), values=c(21,23)) + theme_bw() + scale_x_discrete(labels=c("c1/c1", "c1/c2","c2/c2","c1/c1", "c1/c2","c2/c2")) +
#   theme (axis.text = element_text(size=18), axis.title = element_text(size=20), plot.tag = element_markdown(size=16), legend.text = element_markdown(size=18), legend.title = element_text(size=18),legend.background = element_blank(), legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   labs(x="Genotype", y="PC1") + guides(shape=guide_legend(override.aes = list(size=3)))
# 


#making condition as factor
mod9.pca.loading.gas$Condition <- as.factor(as.character(mod9.pca.loading.gas$Condition))

lmer.mod9.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ Sex:Condition:Genotype.linear+ (1|Vole), data=mod9.pca.loading.gas)
car::Anova(lmer.mod9.linear.pc1.gas, type="3", test.statistic="F")

lmer.mod9.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ (1|Vole), data=mod9.pca.loading.gas)
car::Anova(lmer.mod9.linear.pc1.gas, type="3", test.statistic="F")

lmer.mod9.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+  (1|Vole), data=mod9.pca.loading.gas)
car::Anova(lmer.mod9.linear.pc1.gas, type="3", test.statistic="F")

lmer.mod9.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Sex:Genotype.linear+  (1|Vole), data=mod9.pca.loading.gas)
car::Anova(lmer.mod9.linear.pc1.gas, type="3", test.statistic="F")


#GAS specific module 15, as linear module
vst.mod15.gas <- vst.dds[rownames(vst.dds) %in% moduleGeneID$Module.15, colnames(vst.dds) %in% metaData$NGI.ID[metaData$GAS.paired==1]]

#calculating prinicipal component
rv.mod15.pca.gas <- prcomp(t(assay(vst.mod15.gas)))

#making dataframe with all elements
mod15.pca.loading.gas <- as.data.frame(cbind( "NGI.ID"= rownames(rv.mod15.pca.gas$x),rv.mod15.pca.gas$x))


#merging with metadata
mod15.pca.loading.gas <- merge(mod15.pca.loading.gas, metaData, by="NGI.ID")

#setting pc1 to numeric
mod15.pca.loading.gas$PC1 <- as.numeric(mod15.pca.loading.gas$PC1)

#making conditoin as factor
mod15.pca.loading.gas$Condition <- as.factor(as.character(mod15.pca.loading.gas$Condition))

lmer.mod15.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ Sex:Condition:Genotype.linear+  (1|Vole), data=mod15.pca.loading.gas)
car::Anova(lmer.mod15.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod15.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ Sex:Condition+ (1|Vole), data=mod15.pca.loading.gas)
car::Anova(lmer.mod15.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod15.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Condition:Genotype.linear + Sex:Genotype.linear+ (1|Vole), data=mod15.pca.loading.gas)
car::Anova(lmer.mod15.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod15.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + Sex:Genotype.linear+ (1|Vole), data=mod15.pca.loading.gas)
car::Anova(lmer.mod15.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod15.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ Sex + (1|Vole), data=mod15.pca.loading.gas)
car::Anova(lmer.mod15.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod15.linear.pc1.gas <- lmer(PC1 ~ Condition + Genotype.linear+ (1|Vole), data=mod15.pca.loading.gas)
car::Anova(lmer.mod15.linear.pc1.gas, type="3", test.statistic = "F")

lmer.mod15.linear.pc1.gas <- lmer(PC1 ~ Condition + (1|Vole), data=mod15.pca.loading.gas)
car::Anova(lmer.mod15.linear.pc1.gas, type="3", test.statistic = "F")

#########Modelling effect of tlr2 expression for each condition

#extracting tlr2 expression
tlr2.vst <-  data.frame(vst=assay(vst.dds)[rownames(assay(vst.dds))=="Tlr2"], NGI.ID=colnames(assay(vst.dds)))
tlr2.vst <- merge(tlr2.vst, metaData, by="NGI.ID")

#violin plot for tlr2 expression
tlr2.violin <- ggplot(tlr2.vst, aes(x= Genotype, y=vst, fill=Condition)) + geom_violin(inherit.aes = T, trim = F, show.legend = F, alpha=0.5) + 
  theme_bw() + theme(axis.text = element_markdown(size=18), axis.title=element_markdown(size=20), strip.text.x = element_markdown(size=18), plot.tag = element_markdown(size=17, face="bold"), plot.title = element_markdown(size=20, face = "bold.italic",hjust=0.5),  legend.text = element_markdown(size=16), legend.position = "right", legend.title = element_text(size=17)) + 
  geom_boxplot(show.legend = F,width=0.1) + facet_grid(~Condition, labeller = labeller(Condition=c("Borrelia"="*B. afzelii*","GAS"="*S. pyogenes*" ,"Unstimulated"="Unstimulated"))) + scale_fill_manual(values=c("#78A71F","#71a3d7","#FFB419")) +
  scale_x_discrete(labels=c("c1/c1","c1/c2","c2/c2")) + labs(x="Genotype", y="Gene Expression (vst)") +
  geom_point(aes(x=Genotype, y=vst, group=genoSex, shape=Sex), position = position_jitterdodge(), size=1.5, alpha=1/2,colour = "black", fill = "white", stroke = 1.25) +
  scale_shape_manual(labels=c("F","M"), values=c(21,23)) 

pdf("tlr2_violin.pdf", height = 6, width=8.5)
tlr2.violin
dev.off()


#modelling effect of tlr2 expression for borrelia
tlr2.bor.lm <- lm(vst ~ Genotype.linear + Sex + Genotype.linear:Sex, data=tlr2.vst[tlr2.vst$Condition=="Borrelia",])
car::Anova(tlr2.bor.lm, type="3", test.statistic = "F")

tlr2.bor.lm <- lm(vst ~ Genotype.linear + Sex + Genotype.linear:Sex, data=tlr2.vst[tlr2.vst$Condition=="Borrelia",])
car::Anova(tlr2.bor.lm, type="3", test.statistic = "F")

#modelling effect of tlr2 expression for gas

tlr2.gas.lm <- lm(vst ~ Genotype.linear + Sex + Genotype.linear:Sex, data=tlr2.vst[tlr2.vst$Condition=="GAS",])
car::Anova(tlr2.gas.lm, type="3", test.statistic = "F")

#modelling effect of tlr2 expression for unstimulated
tlr2.uns.lm <- lm(vst ~ Genotype.linear + Sex + Genotype.linear:Sex, data=tlr2.vst[tlr2.vst$Condition=="Unstimulated",])
car::Anova(tlr2.uns.lm, type="3", test.statistic = "F")

tlr2.uns.lm <- lm(vst ~ Genotype.linear + Sex , data=tlr2.vst[tlr2.vst$Condition=="Unstimulated",])
car::Anova(tlr2.uns.lm, type="3", test.statistic = "F")


########################################
######Allele-specific expression########
########################################


#saving file names to import
# files <- paste(getwd(), list.files(path = "../ase/alleleCounts_snp_aseReadCounter", pattern = ".csv", full.names = T), sep="/")
# 
# #importing all files into a list
# files.list <- lapply(files, function(x){read.csv(x, header=T)})
# 
# #renaming list names
# names(files.list) <- gsub(".remapped.keep.rmdup.addRG.sorted.bam.aseReadCount.csv","", gsub("/Users/mridula/Documents/phd/comp/proj2/23042025/../ase/alleleCounts_snp_aseReadCounter/[0-9]+.","", files,perl = T), perl = T)
# 
# #filtering file list to have only uninfected samples
# 
# files.list <- files.list[names(files.list) %in% metaData$NGI.ID]
# 
# #creating and populating a new column called Loci
# files.list <- lapply(files.list, function(x){cbind(x, "Loci"=paste(unlist(x[1]), unlist(x[2]),sep=":"))})
# 
# 
# #creating a master df with all snps across all samples
# loci.dat <- data.frame("Loci"=unique(unlist(sapply(files.list, function(x){paste(unlist(x[1]), unlist(x[2]),sep=":")}))))
# 
# #merging all snps with each df in list
# allCounts <- lapply(files.list, function(x){merge(loci.dat, x[c("refCount","altCount","totalCount","Loci")], by="Loci", all.x=T)})
# 
# #pulling altcounts, ref and total in matrix form
# altCounts <- as.data.frame(bind_cols(lapply(allCounts, "[", , "altCount")))
# rownames(altCounts) <- allCounts$P20454_140$Loci
# 
# refCounts <- as.data.frame(bind_cols(lapply(allCounts, "[", , "refCount")))
# rownames(refCounts) <- allCounts$P20454_140$Loci
# 
# totCounts <- as.data.frame(bind_cols(lapply(allCounts, "[", , "totalCount")))
# rownames(totCounts) <- allCounts$P20454_140$Loci
# 
# #setting loci with < 7 reads expressed as NA ie loci not expressed
# totCounts[totCounts < 7] <- NA
# 
# ####
# #changing reads with counts =0 or < 3 as homozygotes, so don't consider
# altCounts[altCounts ==0] <- NA
# refCounts[refCounts ==0] <- NA
# 
# 
# ####
# #setting NA in ref and alt for thos totCounts with NA
# altCounts[is.na(totCounts)] <- NA
# refCounts[is.na(totCounts)] <- NA
# 
# #changing reads with counts =0 or < 3 as homozygotes, so don't consider
# altCounts[altCounts < 3] <- NA
# refCounts[refCounts < 3] <- NA
# 
# 
# ##update matrices to have NA wherevere necessary; can pull out homozygotes from totCounts
# altCounts[is.na(refCounts)] <- NA
# refCounts[is.na(altCounts)] <- NA
# 
# write.table(cbind(snp=rownames(altCounts), altCounts), "altCounts_snp.tsv", col.names=T, row.names = T, quote = F, sep="\t")
# write.table(cbind(snp=rownames(refCounts), refCounts), "refCounts_snp.tsv", col.names=T, row.names = T, quote = F, sep="\t")

#importing counts for alternate counts
altCounts <- read.table("altCounts_snp.tsv", header=T, row.names = 1)
altCounts <- altCounts[,-1]


#importing counts for reference counts
refCounts <- read.table("refCounts_snp.tsv", header=T, row.names = 1)
refCounts <- refCounts[,-1]

## individuals for which all 3 conditons are present 
pairedVole <- metaData$Vole[metaData$Bor.paired ==1 & metaData$GAS.paired==1 & metaData$Condition=="Unstimulated"]

#splitting by condition
condition <- levels(metaData$Condition)

createDFs <- function(cond){
  listDF <- list()
  for (c in 1:length(cond))
  {
    #subsetting metaData to reorder columns to match individuals in each condition
    met <- metaData[metaData$Condition== cond[c] & metaData$Vole %in% pairedVole, c("NGI.ID","Vole")]
    
    #create reference count df for each condition
    refcounts <- refCounts[, colnames(refCounts) %in% met$NGI.ID]
    
    #reorder by metadata columns
    refcounts <-refcounts[match(met$NGI.ID, colnames(refcounts))]
    
    #renaming colnames
    colnames(refcounts) <- met$Vole
    
    #create alt count df for each condition
    altcounts <- altCounts[, colnames(altCounts) %in% met$NGI.ID]
    
    #reorder column names
    altcounts <- altcounts[match(met$NGI.ID, colnames(altcounts))]
    
    #renaming colnames
    colnames(altcounts) <- met$Vole
    
    #create "missingness" column with sum of NAs; how many loci has missing ("NA") data in totality
    refcounts$missingness <-  rowSums(is.na(refcounts))
    altcounts$missingness <-  rowSums(is.na(altcounts))
    
    #number of samples for each condition
    n.cond = length(met$Vole)
    
    
    #create df for filtered counts; filtering for too too few heterozygotes by removing MAF < 0.15 (so 2pq=0.255)
    refcountsFilt <- refcounts[refcounts$missingness < (n.cond - 0.255*n.cond),]
    altcountsFilt <- altcounts[altcounts$missingness < (n.cond - 0.255*n.cond),]
    
    #creating list names for all 4 DFs
    reflistName <- paste("refCounts", cond[c], sep=".")
    altlistName <- paste("altCounts", cond[c], sep=".")
    
    refFilt.listName <- paste("refFiltCounts", cond[c], sep=".")
    altFilt.listName <- paste("altFiltCounts", cond[c], sep=".")
    
    #assigning dfs to list
    listDF[[reflistName]] <- refcounts
    listDF[[altlistName]] <- altcounts
    listDF[[refFilt.listName]] <- refcountsFilt
    listDF[[altFilt.listName]] <- altcountsFilt
  }
  return(listDF)
}

counts.list <- createDFs(condition)



#check should be TRUE
all(rownames(counts.list$altFiltCounts.Borrelia) %in% rownames(counts.list$refFiltCounts.Borrelia))
all(rownames(counts.list$altFiltCounts.GAS) %in% rownames(counts.list$refFiltCounts.GAS))
all(rownames(counts.list$altFiltCounts.Unstimulated) %in% rownames(counts.list$refFiltCounts.Unstimulated))


#finding common snps to all conditions
commonSNPS <- intersect(intersect(rownames(counts.list$refFiltCounts.Unstimulated), rownames(counts.list$refFiltCounts.GAS)), rownames(counts.list$refFiltCounts.Borrelia))

#filtering for onlycommon SNPs
for (i in condition){
  
  #storing df name
  refdfName <- paste("refFiltCounts",i, sep=".")
  altdfName <- paste("altFiltCounts",i, sep=".")
  
  #fetching only common snps for each df and resaving filtered df
  counts.list[[refdfName]] <- counts.list[[refdfName]][rownames(counts.list[[refdfName]]) %in% commonSNPS,]
  counts.list[[altdfName]] <- counts.list[[altdfName]][rownames(counts.list[[altdfName]]) %in% commonSNPS,]
  
  #sorting by snp
  counts.list[[refdfName]] <- counts.list[[refdfName]][sort(rownames(counts.list[[refdfName]])),]
  counts.list[[altdfName]] <- counts.list[[altdfName]][sort(rownames(counts.list[[altdfName]])),]
  
}

#importing gtf; genes.gtf is subset of genomic.gtf and contains only lines containing "gene"; even pseduogenes are included
bv.gtf <- rtracklayer::import("~/Documents/phd/comp/genome/genes.gtf")

#converting snps grnages object
commonSNPS.gr <- GRanges(seqnames = commonSNPS)

#finding overlaps between common snps and genes
overlaps.snp.gtf <- findOverlaps(commonSNPS.gr, bv.gtf)

#subsetting to contain genes overlapping common snps
commonSNPs.gtf <- unlist(split(bv.gtf[subjectHits(overlaps.snp.gtf)], 1:length(overlaps.snp.gtf)))

#subsetting to contain list of snps present in geneic regions
matchedSNPs <- unlist(split(commonSNPS.gr[queryHits(overlaps.snp.gtf)], 1:length(overlaps.snp.gtf)))


#getting number of snps for each gene
snp.gene.table <- table(commonSNPs.gtf$gene_id)


#matching snps with genes 
snpGene.keyValue <- as.data.frame(cbind("snp"=paste(seqnames(matchedSNPs), ranges(matchedSNPs), sep=":"), "gene"=commonSNPs.gtf$gene_id))



##Added up all the read counts for the snps that separate the two haplotypes 
#and then performed t-test to check for allelic fold change ( afc = log2 (alternate read count/reference read count)) deviation from 0
#Reports 3 lists (one for each condition), with mean AFC and corresponding pvalue

#snps that separate the two hapltoypes; both syn and non.syn snps
key.snp <- c("NW_025965322.1:3120521","NW_025965322.1:3120616","NW_025965322.1:3120958", "NW_025965322.1:3121533","NW_025965322.1:3121588","NW_025965322.1:3121641","NW_025965322.1:3121783","NW_025965322.1:3121807","NW_025965322.1:3121889","NW_025965322.1:3121918","NW_025965322.1:3121976","NW_025965322.1:3122152","NW_025965322.1:3122347")


#restricting to only heterozygous individuals from genotype
het.Voles <- as.character(metaData$Vole[metaData$Vole %in% pairedVole & metaData$Genotype =="2" & metaData$Condition =="Unstimulated"])

tlr2.hap.dat.list <- list()

#subsetting ref and alt read counts for tlr2 snps
for(x in 1:length(condition)){
  #storing df name
  refdfName <- paste("refFiltCounts",condition[x], sep=".")
  altdfName <- paste("altFiltCounts",condition[x], sep=".")
  
  #sum reads for all snps within sample
  ref.tlr2 <- colSums(counts.list[[refdfName]][rownames(counts.list[[refdfName]]) %in% key.snp, colnames(counts.list[[refdfName]])!="missingness"], na.rm = T)
  alt.tlr2 <- colSums(counts.list[[altdfName]][rownames(counts.list[[altdfName]]) %in% key.snp, colnames(counts.list[[altdfName]])!="missingness"], na.rm = T)
  
  #restricting to only heterozygous individuals from genotype
  ref.tlr2 <- ref.tlr2[names(ref.tlr2) %in% het.Voles]
  alt.tlr2 <- alt.tlr2[names(alt.tlr2) %in% het.Voles]
  
  tlr2.hap.dat.list[[refdfName]] <- ref.tlr2
  tlr2.hap.dat.list[[altdfName]] <- alt.tlr2
  
}


tlr2.hap.dat <- as.data.frame(do.call(rbind, tlr2.hap.dat.list))

tlr2.hap.dat$hap <- factor(matrix(unlist(strsplit(rownames(tlr2.hap.dat), split =".", fixed = T)), ncol=2, byrow = T)[,1], levels=c("refFiltCounts", "altFiltCounts")) 

tlr2.hap.dat$Condition <- as.factor(matrix(unlist(strsplit(rownames(tlr2.hap.dat), split =".", fixed = T)), ncol=2, byrow = T)[,2]) 

tlr2.hap.data.melt <- melt(tlr2.hap.dat)
names(tlr2.hap.data.melt)[3:4] <- c("Vole", "readCount")

condition.lab <- c("Borrelia"="*B. afzelii*", "GAS"="*S. pyogenes*","Unstimulated"="Unstimulated")

#lineplot for ref and alt counts by each snp
tlr2.hap.plot <- ggplot(data=tlr2.hap.data.melt, aes(x=hap,y=readCount, color=Condition, group=paste(Condition,Vole))) + geom_point(size=2.5, show.legend = F, alpha=0.65) + geom_line(show.legend = F, linewidth=1, alpha=0.65) +
  labs(x="Haplotype", y="Read count in heterozygotes")+ scale_x_discrete(labels=c("Cluster 1", "Cluster 2"))+
  theme_bw() + theme(axis.title= element_markdown(size=20), axis.text = element_markdown(size=18), strip.text.x = element_markdown(size=18), plot.title = element_markdown(size=20, face = "bold.italic",hjust=0.5))  +
  facet_grid(~Condition, labeller = labeller(Condition=condition.lab)) + scale_color_manual(values=c("#78A71F","#71a3d7","#FFB419"))



##t-test for ase in tlr2 
tlr2.afcFun <- function(cond){
  listDF <- list()
  for (c in 1:length(cond))
  {
    #storing df name
    refdfName <- paste("refFiltCounts",cond[c], sep=".")
    altdfName <- paste("altFiltCounts",cond[c], sep=".")
    
    #splitting counts list df for ref and alt
    ref.counts <- tlr2.hap.dat[refdfName,]
    alt.counts <- tlr2.hap.dat[altdfName,]

    #creating lfc df
    afc <- log2(alt.counts[,1:25]/ref.counts[,1:25])
    
    #calculating mean lfc across all samples
    afc$meanAFC <- rowMeans(afc, na.rm = T)
    
    #t-test to check if lfc across samples is significantly different than 0 (mu=0) and calculating adjusted p-values
    afc$pvalue <- apply(afc[,names(afc) != "meanAFC"], 1, function(x) {t.test(x, mu=0)$p.value})
    afc$padj <- p.adjust(afc$pvalue, method = "fdr")
    
    #creating lfc list name 
    afcListName <- paste("afc", cond[c], sep=".")
    
    #assigning dfs to list
    listDF[[afcListName]] <- afc
  }
  return(listDF)
}

tlr2.afc.list <- tlr2.afcFun(condition)

tlr2.afc.melt <-  bind_rows(tlr2.afc.list)[1:25]

tlr2.afc.melt$Condition <- rownames(tlr2.afc.melt)

tlr2.afc.melt <- melt(tlr2.afc.melt)

tlr2.afc.melt$Condition <- as.factor(tlr2.afc.melt$Condition)

names(tlr2.afc.melt)[2:3] <- c("Vole", "aFC")

tlr2.afc.melt$Condition <- gsub("altFiltCounts.",tlr2.afc.melt$Condition, replacement = "")

tlr2.afc.melt <- merge(tlr2.afc.melt, metaData[metaData$Condition=="Unstimulated",c("Vole","Sex")],by="Vole")

tlr2.hap.boxplot <- ggplot(data=tlr2.afc.melt) + geom_boxplot(aes(x=Condition,y=aFC, fill=Condition), show.legend = F) + 
  labs(x="Condition", y="Log 2 allelic fold change")+ scale_fill_manual(values=c("#78A71F","#71a3d7","#FFB419"))+ scale_x_discrete(labels=c("*B. afzelii*", "*S. pyogenes*","Unstimulated"))+
  geom_point(aes(x=Condition, y=aFC, group=paste(Condition, Sex), shape=Sex), position = position_jitterdodge(), size=2.25, alpha=1/2,colour = "black", fill = "white", stroke = 1.25) +
  scale_shape_manual(labels=c("F","M"), values=c(21,23)) + theme_bw() +
  theme (axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_markdown(size=18), axis.title = element_markdown(size=20), legend.text = element_markdown(size=18), legend.title = element_text(size=18),legend.background = element_blank(), legend.position = "top") +
  guides(shape=guide_legend(override.aes = list(size=3))) 

pdf("ase_tlr2.pdf",height = 10, width=8)
ggarrange(tlr2.hap.plot, tlr2.hap.boxplot, nrow=2, align = "h", heights = c(1,0.75))
dev.off()





