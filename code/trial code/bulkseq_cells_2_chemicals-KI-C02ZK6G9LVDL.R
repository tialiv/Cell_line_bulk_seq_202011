#####################################################2 chemicals########################################################

---
title: "Cell lines bulk RNAseq data analysis - DES and KTZ"
author: "Tianyi"
date: "'r format(Sys,time(), '%d %B%, %Y')'"
output:html_document

---


setwd("/Users/tili/Desktop/Results/PCR/bulk seq_primary_KGN_COV/DESeq2/")

info <- read.csv("forskolin info.csv", header = T, sep = ";")
cells <- read.delim("subreadCounts.txt")
cells <- read.delim("subreadCounts_ens_unstr_new.txt")
cells <- read.delim("forskolin.txt")
cells[is.na(cells)] <- 0 

library(dplyr)
info <- info %>% filter(time != "7d") %>% filter(cell_line == "COV434")
gene <- cells[,1] 
rownames(cells) <- gene 
cells <- cells[,-1] 
# reorder the count matrx using match() function
genomic_idx <- match(info$sample, colnames(cells))
genomic_idx
cells_ordered  <- cells[ ,genomic_idx]

all(info$name == colnames(cells)) 

# adult <- as.data.frame(apply(adult, 2, as.numeric))
# rownames(adult) <- gene

info <- info %>% filter(chemical != "DES")
cells <- re_order(info, cells)

############DES dataset construction#################
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=cells, colData = info, design = ~ cell_line + group) 
hist(log(counts(dds))) #if it is a skewed dataset, needs transformation
colSums(counts(dds)) # how many reads were sequenced for each sample (library size)
d <- info
d$read <- colSums(counts(dds)) # how many reads were sequenced for each sample (library size)
library(ggplot2)
ggplot(d, aes(x=sample, y=read, fill = cell_line)) + geom_bar(stat = "identity") + theme +
  ggtitle("Samples Sequenced Read") + xlab("")
boxplot(log2(counts(dds) + 1), notch = TRUE) 

#########remove low counts##########
keep <- rowSums(counts(dds)) >= 1 
dds <- dds[keep,]
boxplot(log2(counts(dds) + 1), notch = TRUE) 

############Checking the data cluster and plot PCA########
library(ggplot2)
r <- vst(dds, blind = F)
theme <- theme(text=element_text(size=9),
               axis.text.x=element_text(color="black", angle = 90),
               axis.text.y=element_text(color="black"),
               axis.ticks=element_blank(),
               strip.background = element_rect(colour=NA, fill=NA),
               plot.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               axis.line.x = element_line(size=0.4),
               axis.line.y = element_line(size=0.4),
               panel.background = element_blank(),
               legend.position = "bottom")
pca <- plotPCA(r, intgroup = c("cell_line","group"), returnData = T) 
percentVar <- round(100 * attr(pca, "percentVar")) 
ggplot(pca, aes(PC1, PC2, color = group, shape = cell_line)) + ggtitle("PCA") + theme +
  geom_point(size = 3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

# To explore more PCs since plotPCA only plot PC1 and PC2
r_mat <- assay(r) #extract the vst matrix from the object
PCA <- prcomp(t(r_mat))
df <- cbind(info, PCA$x)
ggplot(df) + geom_point(aes(x=PC2, y = PC3, color = chemical, shape = cell_line)) + theme

ggsave(filename="cell bulk seq PCA 201208.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)

##############QC: hierarchical clustering########
r_cor <- cor(r_mat)
library(pheatmap)
df <- as.data.frame(colData(dds)[,c("cell_line","forskolin","media","time")])
pheatmap(r_cor, annotation = df, scale = "column", show_rownames = F,
         clustering_distance_cols = "manhattan")

#############umap###########
library(umap)
library(dplyr)
set.seed(1)
normalized_counts <- assay(r) %>% t()
umap_results <- umap::umap(normalized_counts)
out <- data.frame(umap_results$layout)
umap_plot_df <- data.frame(umap_results$layout) %>%
  tibble::rownames_to_column("sample") %>% 
  dplyr::inner_join(info, by = "sample")

ggplot(umap_plot_df, aes(x = X1, y = X2, color = cell_line, shape = forskolin)) + geom_point(size = 2.5) + theme +
  ggtitle("UMAP")

ggsave(filename="cell bulk seq UMAP 201208.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)

########remove batch effects##########
write.csv(info, file="info.csv",row.names=T)
batch2 <- as.data.frame(pca)
info2 <- read.csv("info.csv", header = T, sep = ";")
removed <- assay(r)
removed <- limma::removeBatchEffect(removed, info2$batch)
assay(r) <- removed
pca_removed <- plotPCA(r, intgroup = "group", returnData = T) 
ggplot(pca_removed, aes(PC1, PC2, color = group)) + ggtitle("PCA") + theme +
  geom_point(size = 3)

ggsave(filename="adult COCs PCA_removed batch effect 201126.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)

########setting reference group##########
dds$group <- relevel(dds$group, ref = "DMSO")
dds$forskolin <- relevel(dds$forskolin, ref = "no")

############ voom transformation ################
cell_line <- factor(info$cell_line)
group <- factor(info$group)
design <- model.matrix(~0 + cell_line + group)
voom_norm <- limma::voom(cells, design, lib.size = colSums(cells))

voom_cells <- voom_norm$E

###############differential expression analysis######
dds <-DESeq(dds)
plotDispEsts(dds)
res <- results(dds, contrast = c("group", "DES_10-10M", "DMSO"), independentFiltering = T, pAdjustMethod = "fdr")
res <- results(dds, contrast = c("forskolin", "yes", "no"), independentFiltering = T, pAdjustMethod = "fdr")
res <- results(dds, contrast = c("media", "clear", "red"), independentFiltering = T, pAdjustMethod = "fdr")

comparison <- list()
for (i in chemical) {
    files <- paste0("/Users/tili/Desktop/Cell_line_bulk_seq_202011/data/processed/info_",                        
                    i,"_",j,".csv")
    info <- read.csv(file = files, sep = ",")
    info$group <- factor(info$group)
    groups <- factor(levels(info$group)[levels(info$group) != "DMSO"]) 
    for (k in groups) {
      print(paste0("Differential expressed genes for"," ", k, " ", "vs DMSO", " ", i))
      res <- get_res(k)
      comparison[[paste0("res_", k, " ", i)]] <- check_sig(res, 0.1) %>% as.data.frame() %>% 
        mutate(geneid = rownames(.))
  }
} 

contrast <- c("group", "DES_10-10M", "DMSO") # control group should be placed in the last 

# for the functional analysis such as GO analysis, shrinkage of the fold changes is recommended
res_shrunk <- results(dds, contrast = contrast, alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = 2, type = "apeglm", lfcThreshold = 1)
plotMA(res)
plotMA(res_shrunk)

res.sig <- res[which(res$padj<0.05),]
res.sig <- res.sig[order(res.sig$padj),]
summary(res.sig)
write.table(res.sig,file="significant_diff_expr.txt", sep="\t",quote=F,row.names=T)

res.sig.shrunk <- data.frame(res_shrunk)
res.sig.shrunk <- res_shrunk[which(res_shrunk$svalue<0.05),]
res.sig.shrunk <- res.sig.shrunk[order(res.sig.shrunk$svalue),]
summary(res.sig.shrunk)

###########getting common gene symbols and names if you only have the gene id#####
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
ensembl_genes <- rownames(res.sig)
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), mart = mart)
res.sig <- as.data.frame(res.sig)
res.sig$geneID <- rownames(res.sig)
final <- merge(res.sig, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)

ensembl_genes <- rownames(res.sig.shrunk)
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart = mart)
res.sig.shrunk <- as.data.frame(res.sig.shrunk)
res.sig.shrunk$geneID <- rownames(res.sig.shrunk)
final <- merge(res.sig.shrunk, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)

# listFilters(mart)[1:10,]
# listFilterValues(mart, filter = "chromosome_name")[1:20]
# res.sig$GeneSymbol <- output$external_gene_name


############Functional analysis########

library(org.Hs.eg.db)

# Choose to look at ontology CC(cellular component), MF(molecular function) or BP(biological process)

library(DOSE)
library(pathview)
library(clusterProfiler)

all_genes <- data.frame(res_shrunk)
ensembl_genes <- rownames(all_genes)
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart = mart)
all_genes <- as.data.frame(all_genes)
all_genes$geneID <- rownames(all_genes)
all <- merge(all_genes, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)
summary(all)

ego <- enrichGO(gene = final$geneID, universe = all$geneID, keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, readable = TRUE)

cluster_summary <- data.frame(ego)
# To save the object as .rda file, load it later with load() function
# save(ego, file="results/ego.rda")
# ego <- load(file="results/ego.rda")

dotplot(ego, showCategory=50) # dotplot showes the number of genes associated with the first 50 terms and p values

emapplot(ego, showCategory = 50) # enrichment GO plot to show relationship between most significantly enriched GO terms

foldchange <- final$log2FoldChange
names(foldchange) <- final$external_gene_name

#gene set analysis#
gse <- gseGO(geneList = foldchange, ont = "ALL", keyType = "ENSEMBL",
             nPerm= 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05,
             verbose = T, pAdjustMethod = "fdr")

# category netplot showes relationship between genes associated with most significant GO terms and their fold changes
cnetplot(ego, categorySize="pvalue", showCategory = 5, foldChange=foldchange, vertex.label.font=6)

# Manully choose the interested GO term from the original dataset by subsetting, only choose no.1,2,5 and no.6
ego2 <- ego
ego2@result <- ego@result[c(1,2,5,6),]
cnetplot(ego2, categorySize="pvalue", showCategory = 5, foldChange=foldchange, vertex.label.font=6)

########Gene set analysis#######3
res.sig.shrunk <- data.frame(res_shrunk)
res.sig.shrunk$geneID <- rownames(res.sig.shrunk)
foldchange <- res.sig.shrunk$log2FoldChange
names(foldchange) <- res.sig.shrunk$geneID
foldchange <- na.omit(foldchange)
foldchange <- sort(foldchange, decreasing = T)

#gene set analysis#
gse <- gseGO(geneList = foldchange, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.1,  minGSSize = 100,
             maxGSSize = 500, verbose = T, pAdjustMethod = "none", OrgDb = org.Hs.eg.db)
gseaGO_results <- gse@result

gseaplot(gse, geneSetID = 'GO:0007423')
dotplot(gse, showCategory = 50, split = ".sign") + facet_grid(.~.sign) 

############# Msigdb gene set analysis ###########3

library(msigdbr)
msigdbr_species()
msigdbr_collections()

msigHs_c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  select(gs_name, entrez_gene) %>% 
  as.data.frame()

msigHs_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, entrez_gene) %>% 
  as.data.frame()

library(biomaRt)
mart <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id","description"), mart = mart)
res.sig <- as.data.frame(res.sig)
res.sig$geneID <- rownames(res.sig)
final <- merge(res.sig, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)

enrich <- enricher(gene = gene$entrezgene_id, TERM2GENE = msigHs_h,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 500)
enrich <- data.frame(enrich)
write.csv(enrich, file = "gse with hallmark_0_05_BH_with 4107 identified gene.csv")

genelist <- gene %>% 
  dplyr::select(entrezgene_id, log2FoldChange) 

List <- genelist$log2FoldChange
names(List) <- as.character(genelist$entrezgene_id)
List <- sort(List, decreasing = T)

edo <- GSEA(geneList = List, minGSSize = 10, maxGSSize = 500,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            TERM2GENE = msigHs_h)
top10 <- data.frame(edo)

gseaplot(edo, geneSetID = 1, by = "runningScore", title = edo$Description[1])
gseaplot(edo, geneSetID = 1, by = "preranked", title = edo$Description[1])
enrichplot::gseaplot2(edo, geneSetID = 1:4, title = edo$Description[1:4])


###############visualization##############
library(pheatmap)
NR <- read.csv("nuclear receptor.csv", sep = ";", header = T)

select <- NR$ensembl_gene_id
library(tidyverse)
normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  left_join(output, by=c("gene" = "ensembl_gene_id"))

id <- match(NR$geneID, normalized_counts$gene)
id
normalized_counts <- normalized_counts[id,] 
d1 <- na.omit(normalized_counts)
d3 <- data.frame(d1[,2:43]) 
rownames(d3) <- d1$external_gene_name
df <- as.data.frame(colData(dds)[,c("forskolin","cell_line","media")])
pheatmap(d3, cluster_rows = F, annotation = df,
         scale = "none", fontsize_row = 10, clustering_distance_cols = "manhattan",
         height = 20, main = "Nuclear receptors expression")

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("cell_line","group")])
pheatmap(assay(r)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df)


###########plot individual genes################

des <- read.csv("/Users/tili/Desktop/Cell_line_bulk_seq_202011/results/figures/des.csv", header = T, sep = ";")
ktz <- read.csv("/Users/tili/Desktop/Cell_line_bulk_seq_202011/results/figures/ktz.csv", header = T, sep = ";")

col <- c("DMSO" = "#f0f0f0",#"KTZ_10-9M" = "#d9f0a3",
         #"KTZ_10-5M" = "#006d2c")
         "DES_10-10M" = "#b3cde3","DES_10-6M" = "#88419d") 

theme <- theme(text=element_text(size=14),
               axis.text.x=element_text(color="black", angle = 90, vjust = 0.3, hjust = 0.3),
               axis.text.y=element_text(color="black"),
               axis.ticks.x.bottom=element_blank(),
               #axis.ticks.length.x=unit(-.25, "cm"),
               strip.background = element_rect(colour=NA, fill=NA),
               plot.background = element_blank(),
               #panel.grid.major = element_blank(),
               #panel.grid.minor = element_blank(),
               panel.spacing.x=unit(.5, "lines"),
               panel.border = element_rect(color = "black", fill = NA, size = 1),
               axis.line.x = element_line(size=0.4),
               axis.line.y = element_line(size=0.4),
               legend.key.size = unit(2, 'lines'),
               panel.background = element_blank()
)

des_plot <- list()
for (i in 1:23) {
  geneid <- des[i, 2]
  gene <- des[i,1]
  d <- plotCounts(dds, gene = geneid, intgroup = c("group","cell_line"), normalized = T, returnData=TRUE)
  des_plot[[gene]] <- ggplot(d, aes(x = group, y = count, fill = group)) + geom_boxplot() +
    geom_point(size = 2, alpha = 0.5) + facet_grid(~ cell_line) + ylab("Normalized counts") + xlab(" ") +
    theme + ggtitle(gene) + scale_x_discrete(limit = c("DMSO", "DES_10-10M", "DES_10-6M"
                                                       #, "KTZ_10-9M", "KTZ_10-5M"
                                                       )) +
    scale_fill_manual(values = col) + theme(legend.position="none")
}

cowplot::plot_grid(nrow = 3, ncol = 3, plotlist = des_plot[1:9], 
                   labels = "AUTO")
cowplot::plot_grid(nrow = 1, ncol = 2, plotlist = des_plot[22:23], 
                   labels = "AUTO")

des_df <- cbind(des_plot[[1]]$data, des_plot[[2]]$data$count)
des_df <- cbind(des_df, des_plot[[3]]$data$count)
des_df <- cbind(des_df, des_plot[[4]]$data$count)
des_df <- cbind(des_df, des_plot[[5]]$data$count)
des_df <- cbind(des_df, des_plot[[6]]$data$count)
des_df <- cbind(des_df, des_plot[[7]]$data$count)
des_df <- cbind(des_df, des_plot[[8]]$data$count)
des_df <- cbind(des_df, des_plot[[9]]$data$count)
des_df <- cbind(des_df, des_plot[[10]]$data$count)
des_df <- cbind(des_df, des_plot[[11]]$data$count)
des_df <- cbind(des_df, des_plot[[12]]$data$count)
des_df <- cbind(des_df, des_plot[[13]]$data$count)
des_df <- cbind(des_df, des_plot[[14]]$data$count)
des_df <- cbind(des_df, des_plot[[15]]$data$count)
des_df <- cbind(des_df, des_plot[[16]]$data$count)
des_df <- cbind(des_df, des_plot[[17]]$data$count)

colnames(des_df) <- c("ERBIN", "group", "cell_line", des$gene[2:17])

p6 <- des_df %>% mutate(sample = rownames(.)) %>% reshape2::melt() %>% 
  group_by(group, variable, cell_line) %>% 
  summarize(average = mean(value)) %>% 
  mutate(gro = paste0(group, "_", cell_line)) %>% 
  ggplot(aes(gro, variable, fill = log2(average))) + geom_tile() + theme +
  scale_fill_gradient2(low = "#075AFF", high = "#FF0000") +
  coord_fixed() + scale_x_discrete(limit = c("DMSO_COV434", "DES_10-10M_COV434", "DES_10-6M_COV434",
                                             "DMSO_KGN", "DES_10-10M_KGN", "DES_10-6M_KGN",
                                             "DMSO_Primary", "DES_10-10M_Primary", "DES_10-6M_Primary")) +
  xlab(" ") + ylab(" ") + ggtitle("DES")

ktz_plot <- list()
for (i in 1:18) {
  geneid <- ktz[i, 2]
  gene <- ktz[i,1]
  d <- plotCounts(dds, gene = geneid, intgroup = c("group","cell_line"), normalized = T, returnData=TRUE)
  ktz_plot[[gene]] <- ggplot(d, aes(x = group, y = count, fill = group)) + geom_boxplot() +
    geom_point(size = 2, alpha = 0.5) + facet_grid(~ cell_line) + ylab("Normalized counts") + xlab(" ") +
    theme + ggtitle(gene) + scale_x_discrete(limit = c("DMSO", 
                                                       #"DES_10-10M", "DES_10-6M",
                                                       "KTZ_10-9M", "KTZ_10-5M")) +
    scale_fill_manual(values = col) + theme(legend.position="none")
}

cowplot::plot_grid(nrow = 3, ncol = 3, plotlist = ktz_plot[1:9], 
                   labels = "AUTO")
cowplot::plot_grid(nrow = 3, ncol = 3, plotlist = ktz_plot[10:18], 
                   labels = "AUTO")

ktz_df <- cbind(ktz_plot[[1]]$data, ktz_plot[[2]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[3]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[4]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[5]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[6]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[7]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[8]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[9]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[10]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[11]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[12]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[13]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[14]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[15]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[16]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[17]]$data$count)
ktz_df <- cbind(ktz_df, ktz_plot[[18]]$data$count)

colnames(ktz_df) <- c("DSTYK", "group", "cell_line", ktz$gene[2:18])

p5 <- ktz_df %>% mutate(sample = rownames(.)) %>% reshape2::melt() %>% 
  group_by(group, variable, cell_line) %>% 
  summarize(average = mean(value)) %>% 
  mutate(gro = paste0(group, "_", cell_line)) %>% 
  ggplot(aes(gro, variable, fill = log2(average))) + geom_tile() + theme +
  scale_fill_gradient2(low = "#075AFF", high = "#FF0000") +
  coord_fixed() + scale_x_discrete(limit = c("DMSO_COV434", "KTZ_10-9M_COV434", "KTZ_10-5M_COV434",
                                             "DMSO_KGN", "KTZ_10-9M_KGN", "KTZ_10-5M_KGN",
                                             "DMSO_Primary", "KTZ_10-9M_Primary", "KTZ_10-5M_Primary")) +
  xlab(" ") + ylab(" ") + ggtitle("KTZ")

p6 + p5 + plot_annotation(tag_levels = "A") 

d <- plotCounts(dds, gene="ENSG00000169083", intgroup=c("group","cell_line"), normalized = T,
                returnData=TRUE) # gene argument can be either row name or numeric index
ggplot(d, aes(x = group, y = count, fill = group)) + geom_boxplot() +
  geom_point(size = 2, alpha = 0.5) + facet_grid(~ cell_line) + ylab("Normalized counts") + xlab(" ") +
  theme + ggtitle("AR") + scale_x_discrete(limit = c("DMSO", 
                                                     "DES_10-10M", "DES_10-6M"#,
                                                     #"KTZ_10-9M", "KTZ_10-5M"
                                                     )) +
  scale_fill_manual(values = col) + theme(legend.position="none")

ggplot(d, aes(x=group, y=count, color=group)) + 
  geom_boxplot() + 
  scale_y_log10(breaks=c(25,100,400)) + theme + ggtitle("PHF7")

ggsave(filename="adult COCs DES_0.003 Golga2 201126.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)



sessionInfo()

