
################################################ DES ##################################################################
#######################################################################################################################

setwd("/Users/tili/Desktop/Results/PCR/bulk seq_primary_KGN_COV/DESeq2/")

info <- read.csv("sample info.csv", header = T, sep = ";")
cells <- read.delim("DES.txt")
cells[is.na(cells)] <- 0 

info <- read.csv("info.csv", header = T, sep = ";")
gene <- cells[,1] 
rownames(cells) <- gene 
cells <- cells[,-1] 
library(dplyr)
info <- info %>% filter(chemical != "KTZ") %>% filter(cell_line == "Primary")

# reorder the count matrx using match() function
genomic_idx <- match(info$sample, colnames(cells))
genomic_idx
cells_ordered  <- cells[,genomic_idx]

all(info$sample == colnames(cells_ordered)) 

############DES dataset construction#################
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=cells_ordered, colData = info, design = ~ cell_line + group) 
dds <- DESeqDataSetFromMatrix(countData=cells_ordered, colData = info, design = ~ group) 
hist(log(counts(dds))) #if it is a skewed dataset, needs transformation
d <- info 
d$read <- colSums(counts(dds)) # how many reads were sequenced for each sample (library size)
ggplot(d, aes(group, read)) + theme + geom_boxplot() + ggtitle("Primary Aligned reads") +
  xlab("") + scale_x_discrete(limit = c("DMSO", "DES_10-10M", "DES_10-6M"))
boxplot(log2(counts(dds) + 1), notch = TRUE) 

#########remove low counts##########
keep <- rowSums(counts(dds)) >= 1 
dds <- dds[keep,]

boxplot(log2(counts(dds) + 1), notch = TRUE) 

############Checking the data cluster and plot PCA########
library(ggplot2)
r <- rlog(dds, blind = F)
theme <- theme(text=element_text(size=14),
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
               legend.position = "right")
pca <- plotPCA(r, intgroup = c("group"), returnData = T) 
pca$read <- d$read
percentVar <- round(100 * attr(pca, "percentVar")) 
ggplot(pca, aes(PC1, PC2, color = group, label = name)) + ggtitle("PCA - COV") + theme +
  geom_point(size = 3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
 geom_text(aes(label = name),hjust=.5, vjust=-1, size=3)

# To explore more PCs since plotPCA only plot PC1 and PC2
r_mat <- assay(r) #extract the vst matrix from the object
PCA <- prcomp(t(r_mat))
df <- cbind(info, PCA$x)
percentVar <- round(100 * attr(PCA, "percentVar")) 
ggplot(df) + geom_point(aes(x=PC3, y = PC4, color = group)) + theme

ggsave(filename="cell bulk seq PCA 201208.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)

##############QC: hierarchical clustering########
r_cor <- cor(r_mat)
library(pheatmap)
df <- as.data.frame(colData(dds)[,c("cell_line","group")])

pheatmap(r_cor, scale = "column", clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan")

#############umap###########
library(umap)
library(dplyr)
set.seed(1)
normalized_counts <- assay(r) %>% t()
umap_results <- umap::umap(normalized_counts, n_neighbors = 3)
out <- data.frame(umap_results$layout)
umap_plot_df <- data.frame(umap_results$layout) %>%
  tibble::rownames_to_column("sample") %>% 
  dplyr::inner_join(info, by = "sample")

ggplot(umap_plot_df, aes(x = X1, y = X2, color = group, shape = cell_line)) + geom_point(size = 2.5) + theme +
  ggtitle("UMAP")

ggsave(filename="cell bulk seq UMAP 201223.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)

########setting reference group##########
dds$group <- relevel(dds$group, ref = "DMSO")
dds$exposure <- relevel(dds$exposure, ref = "control")

###############differential expression analysis######
dds <-DESeq(dds)
plotDispEsts(dds)
resultsNames(dds)
res <- results(dds, contrast = c("group", "DES_10-10M", "DMSO"), independentFiltering = T, pAdjustMethod = "fdr")
res <- results(dds, contrast = c("exposure", "exposed", "control"), independentFiltering = T, pAdjustMethod = "fdr")

res.sig <- res[which(res$padj<0.05),]
res.sig <- res.sig[order(res.sig$padj),]
summary(res.sig)
plotMA(res)
write.csv(final,file="DMSO_vs_DES_10 COV_20210414.csv")

contrast <- c("group", "DES_10-6M", "DMSO") # control group should be placed in the last 

# for the functional analysis such as GO analysis, shrinkage of the fold changes is recommended
res_shrunk <- results(dds, contrast = contrast, alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "group_DES_10.6M_vs_DMSO", type = "apeglm")
plotMA(res_shrunk)

res.sig.shrunk <- data.frame(res_shrunk)
res.sig.shrunk <- res_shrunk[which(res_shrunk$padj<0.05),]
res.sig.shrunk <- res.sig.shrunk[order(res.sig.shrunk$padj),]
summary(res.sig.shrunk)

###########getting common gene symbols and names if you only have the gene id#####
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"), mart = mart)
res.sig <- as.data.frame(res.sig)
res.sig$geneID <- rownames(res.sig)
final <- merge(res.sig, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)

ensembl_genes <- rownames(res.sig.shrunk)
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"), mart = mart)
res.sig.shrunk <- as.data.frame(res.sig.shrunk)
res.sig.shrunk$geneID <- rownames(res.sig.shrunk)
final <- merge(res.sig.shrunk, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)

# listFilters(mart)[1:10,]
# listFilterValues(mart, filter = "chromosome_name")[1:20]
# res.sig$GeneSymbol <- output$external_gene_name

##############venn diagram#############3
library(VennDiagram)
setwd("/Users/tili/Desktop/Results/PCR/bulk seq_primary_KGN_COV/DESeq2/DES Results/")
DES10_primary <- read.delim("DMSO vs DES 10-10M primary.txt")
DES6_primary <- read.delim("DMSO vs DES 10-6M primary.txt")

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(DES10_primary$external_gene_name, DES6_primary$external_gene_name),
  category.names = c("DES_10" , "DES_6"),
  filename = 'Primary_DES_venn_diagramm.png',
  output=TRUE
)


############Functional analysis########

library(org.Hs.eg.db)

# Choose to look at ontology CC(cellular component), MF(molecular function) or BP(biological process)

library(DOSE)
library(pathview)
library(clusterProfiler)

all_genes <- data.frame(res_shrunk)
ensembl_genes <- rownames(all_genes)
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), mart = mart)
all_genes <- as.data.frame(all_genes)
all_genes$geneID <- rownames(all_genes)
all <- merge(all_genes, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)


ego <- enrichGO(gene = final$geneID, universe = all$geneID, keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", 
                qvalueCutoff = 0.1, readable = TRUE)

cluster_summary <- data.frame(ego)
cluster_summary$number <- seq(1,517)
# To save the object as .rda file, load it later with load() function
# save(ego, file="results/ego.rda")
# ego <- load(file="results/ego.rda")

dotplot(ego, showCategory=50) # dotplot showes the number of genes associated with the first 50 terms and p values

emapplot(ego, showCategory = 50) # enrichment GO plot to show relationship between most significantly enriched GO terms

foldchange <- final$log2FoldChange
names(foldchange) <- final$external_gene_name

# category netplot showes relationship between genes associated with most significant GO terms and their fold changes
cnetplot(ego, categorySize="pvalue", showCategory = 5, foldChange=foldchange, vertex.label.font=6)

# Manully choose the interested GO term from the original dataset by subsetting, only choose no.1,2,5 and no.6
ego2 <- ego
ego2@result <- ego@result[c(45,48),]
cnetplot(ego2, categorySize="pvalue", showCategory = 5, foldChange=foldchange, vertex.label.font=6)


########Gene set analysis##############
res.sig.shrunk <- data.frame(res_shrunk)
res.sig.shrunk$geneID <- rownames(res.sig.shrunk)
foldchange <- res.sig.shrunk$log2FoldChange
names(foldchange) <- res.sig.shrunk$geneID
foldchange <- na.omit(foldchange)
foldchange <- sort(foldchange, decreasing = T)

#gene set analysis#
gse <- gseGO(geneList = foldchange, ont = "ALL", keyType = "ENSEMBL", pvalueCutoff = 0.1,  minGSSize = 10,
             maxGSSize = 500, verbose = T, pAdjustMethod = "none", OrgDb = org.Hs.eg.db)
gseaGO_results <- gse@result
write.csv(gseaGO_results, file="gse COV_DES_10_vs_DMSO.csv",row.names=T)

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
res_shrunk <- as.data.frame(res_shrunk)
res_shrunk$geneID <- rownames(res_shrunk)
final <- merge(res_shrunk, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F) 

enrich <- enricher(gene = final$entrezgene_id, TERM2GENE = msigHs_h,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 500)
enrich <- data.frame(enrich)
write.csv(enrich, file = "enricher with hallmark_0_05_BH_with DES6 primary gene_20210414.csv")

genelist <- final %>% 
  dplyr::select(entrezgene_id, log2FoldChange) 

List <- genelist$log2FoldChange
names(List) <- as.character(genelist$entrezgene_id)
List <- sort(List, decreasing = T)

edo <- GSEA(geneList = List, minGSSize = 10, maxGSSize = 500,
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            TERM2GENE = msigHs_h)
top10 <- data.frame(edo)
write.csv(top10, file = "gsea with hallmark_1_BH_with DES6 primary gene_20210414.csv")

gseaplot(edo, geneSetID = 1, by = "runningScore", title = edo$Description[1])
gseaplot(edo, geneSetID = 1, by = "preranked", title = edo$Description[1])
enrichplot::gseaplot2(edo, geneSetID = 1, title = edo$Description[2])

ggsave(filename="gsea plot with hallmark_0_05_BH_with DES6 KGN gene_20210414.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)

###############visualization##############
library(pheatmap)
NR <- read.csv("/Users/tili/Desktop/Results/PCR/bulk seq_primary_KGN_COV/DESeq2/nuclear receptor.csv", sep = ";", header = T)
NR <- merge(output, NR, by.x = "ensembl_gene_id", by.y = "geneID", all.x = F, all.y = T)

NR <- na.omit(NR)
select <- NR$ensembl_gene_id
library(tidyverse)
normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  left_join(output, by=c("gene" = "ensembl_gene_id"))

id <- match(NR$ensembl_gene_id, normalized_counts$gene)
id
normalized_counts <- normalized_counts[id,] 
d1 <- na.omit(normalized_counts)
d3 <- d1[,2:10] 
rownames(d3) <- d1$external_gene_name
df <- as.data.frame(colData(dds)[,c("group")])
pheatmap(d3, cluster_rows = T,
         scale = "none", fontsize_row = 10, 
         height = 20, main = "COV434 nuclear receptors expression")

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("exposure","group")])
pheatmap(assay(r)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df, fontsize_row = 10, 
         height = 20)


###########plot individual genes################

d <- plotCounts(dds, gene="ENSG00000100883", intgroup=c("group","cell_line"), normalized = T,
                returnData=TRUE) # gene argument can be either row name or numeric index
ggplot(d, aes(x=group, y=count)) + geom_boxplot() + facet_grid(~cell_line) +
  geom_point(size=2, alpha = 0.5) + xlab("") + ylab("normalized count") +
  theme + ggtitle("SRP54")

ggplot(d, aes(x=group, y=count)) + 
  geom_boxplot() + xlab("") + theme + ggtitle("ERBIN KGN") + geom_point(position=position_jitter(w=0.1,h=0), alpha = 0.5)
  scale_y_log10(breaks=c(25,100,400)) + theme + ggtitle("ERBIN")

ggsave(filename="adult COCs DES_0.003 Golga2 201126.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)



sessionInfo()
