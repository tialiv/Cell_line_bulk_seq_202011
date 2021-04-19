
##################################################### KTZ ###############################################################
#########################################################################################################################

setwd("/Users/tili/Desktop/Results/PCR/bulk seq_primary_KGN_COV/DESeq2/")

info <- read.csv("sample info.csv", header = T, sep = ";")
cells <- read.delim("KTZ.txt")
cells[is.na(cells)] <- 0 

info <- read.csv("info.csv", header = T, sep = ";")
gene <- cells[,1] 
rownames(cells) <- gene 
cells <- cells[,-1] 
library(dplyr)
info <- info %>% filter(chemical != "DES") %>% filter(cell_line == "Primary")

# reorder the count matrx using match() function
genomic_idx <- match(info$sample, colnames(cells))
genomic_idx
cells_ordered  <- cells[,genomic_idx]

all(info$sample == colnames(cells_ordered)) 

############DES dataset construction#################
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=cells_ordered, colData = info, design = ~ group) 
hist(log(counts(dds))) #if it is a skewed dataset, needs transformation
d <- info 
d$read <- colSums(counts(dds)) # how many reads were sequenced for each sample (library size)
ggplot(d, aes(sample, read, fill=cell_line)) + theme + geom_bar(stat = "identity") + ggtitle("Aligned reads") +
  xlab("") + scale_x_discrete(limit = c("DMSO", "KTZ_10-9M", "KTZ_10-5M"))
boxplot(log2(counts(dds) + 1), notch = TRUE) 

#########remove low counts##########
keep <- rowSums(counts(dds)) >= 1 
dds <- dds[keep,]

boxplot(log2(counts(dds) + 1), notch = TRUE) 

############Checking the data cluster and plot PCA########
library(ggplot2)
r <- vst(dds, blind = F)
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
ggplot(pca, aes(PC1, PC2, color = group, label = name)) + ggtitle("PCA - KGN") + theme +
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

pheatmap(r_cor, scale = "column", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

#############umap###########
library(umap)
library(dplyr)
set.seed(10000)
normalized_counts <- assay(r) %>% t()
umap_results <- umap::umap(normalized_counts)
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
res <- results(dds, contrast = c("group", "KTZ_10-9M", "DMSO"), independentFiltering = T, pAdjustMethod = "fdr")
res <- results(dds, contrast = c("exposure", "exposed", "control"), independentFiltering = T, pAdjustMethod = "fdr")

res.sig <- res[which(res$padj<0.05),]
res.sig <- res.sig[order(res.sig$padj),]
summary(res.sig)
plotMA(res)
write.table(final,file="DMSO vs KTZ 10-9M primary.txt", sep="\t",quote=F,row.names=T)

contrast <- c("group", "KTZ_10-9M", "DMSO") # control group should be placed in the last 

# for the functional analysis such as GO analysis, shrinkage of the fold changes is recommended
res_shrunk <- results(dds, contrast = contrast, alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "group_KTZ_10.9M_vs_DMSO", type = "apeglm")
plotMA(res_shrunk)

res.sig.shrunk <- data.frame(res_shrunk)
res.sig.shrunk <- res_shrunk[which(res_shrunk$padj<0.05),]
res.sig.shrunk <- res.sig.shrunk[order(res.sig.shrunk$padj),]
summary(res.sig.shrunk)

###########getting common gene symbols and names if you only have the gene id#####
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "description"), mart = mart)
res.sig <- as.data.frame(res.sig)
res.sig$geneID <- rownames(res.sig)
final <- merge(res.sig, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)

ensembl_genes <- rownames(res.sig.shrunk)
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "description"), mart = mart)
res.sig.shrunk <- as.data.frame(res.sig.shrunk)
res.sig.shrunk$geneID <- rownames(res.sig.shrunk)
final <- merge(res.sig.shrunk, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F)

# listFilters(mart)[1:10,]
# listFilterValues(mart, filter = "chromosome_name")[1:20]
# res.sig$GeneSymbol <- output$external_gene_name

##############venn diagram#############3
library(VennDiagram)
setwd("/Users/tili/Desktop/Results/PCR/bulk seq_primary_KGN_COV/DESeq2/KTZ Results/")
KTZ9_primary <- read.delim("DMSO vs KTZ 10-9M primary.txt")
KTZ5_primary <- read.delim("DMSO vs KTZ 10-5M primary.txt")
KTZ9_COV <- read.delim("DMSO vs KTZ 10-9M COV.txt")
KTZ5_COV <- read.delim("DMSO vs KTZ 10-5M COV.txt")
KTZ9_KGN <- read.delim("DMSO vs KTZ 10-9M KGN.txt")
KTZ5_KGN <- read.delim("DMSO vs KTZ 10-5M KGN.txt")

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(KTZ9_KGN$external_gene_name, KTZ5_KGN$external_gene_name),
  category.names = c("KTZ_9" , "KTZ_5"),
  filename = 'KGN_KTZ_venn_diagramm.png',
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


########Gene set analysis#################
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
res_shrunk <- as.data.frame(res_shrunk)
res_shrunk$geneID <- rownames(res_shrunk)
final <- merge(res_shrunk, output, by.x = "geneID", by.y = "ensembl_gene_id", all.x = T, all.y = F) 

enrich <- enricher(gene = final$entrezgene_id, TERM2GENE = msigHs_h,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 500)
enrich <- data.frame(enrich)
write.csv(enrich, file = "enricher with hallmark_0_05_BH_with KTZ9 primary gene_20210414.csv")

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
write.csv(top10, file = "gsea with hallmark_1_BH_with KTZ9 primary gene_20210414.csv")

gseaplot(edo, geneSetID = 1, by = "runningScore", title = edo$Description[1])
gseaplot(edo, geneSetID = 1, by = "preranked", title = edo$Description[1])
enrichplot::gseaplot2(edo, geneSetID = 7:13, title = edo$Description[7:13])

ggsave(filename="gsea plot with hallmark_0_05_BH_with DES6 KGN gene_20210414.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)


###############visualization##############
library(pheatmap)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("cell_line","group")])
pheatmap(assay(r)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df)


###########plot individual genes################

d <- plotCounts(dds, gene="ENSG00000188191", intgroup=c("group"), normalized = T,
                returnData=TRUE) # gene argument can be either row name or numeric index
ggplot(d, aes(x=group, y=count)) + geom_boxplot() + # facet_grid(~cell_line) +
  geom_point(size=2, alpha = 0.5) + xlab("") + ylab("normalized count") +
  theme + ggtitle("PRKAR1B")

ggplot(d, aes(x=group, y=count, color=group)) + 
  geom_boxplot() + xlab("") +
  scale_y_log10(breaks=c(25,100,400)) + theme + ggtitle("CASK")

ggsave(filename="adult COCs DES_0.003 Golga2 201126.jpeg", 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, 
       dpi = 300)



sessionInfo()
