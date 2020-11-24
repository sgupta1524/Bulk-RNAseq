#########################Run these on command line/ bash####################################################
#indexing the transcriptome
salmon index -t chr22_transcripts.fa -i chr22_index
#quantifying the reads
for i in *_R1.fastq.gz
do
prefix=$(basename $i _R1.fastq.gz)
salmon quant -i chr22_index --libType A \
-1 ${prefix}_R1.fastq.gz -2 ${prefix}_R2.fastq.gz -o quant/${prefix};
done
############################################################################################

library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(readr)

##Linking the transcript names to gene names
setwd('~/run2/')
txdb <- makeTxDbFromGFF("~/run2/hg38.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
head(tx2gene)

##Importing salmon quantification
samples <- read.table("samples.txt", header = TRUE)
files <- file.path("quant2", samples$sample, "quant2.sf")
names(files) <- paste0(samples$sample)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,  dropInfReps=TRUE, , ignoreTxVersion= TRUE, ignoreAfterBar=TRUE)
tail(txi.salmon$counts)

##Differential expression
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
plotDispEsts(dds, main="Dispersion plot")
rld <- rlogTransformation(dds)
head(assay(rld))

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition", "sample")])
ntd <- normTransform(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

library(ggplot2)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathview)
library(gage)
library(gageData)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$name <- mapIds(org.Hs.eg.db,
                   keys=row.names(res),
                   column="GENENAME",
                   keytype="ENSEMBL",
                   multiVals="first")

head(res)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

foldchanges <- res$log2FoldChange
names(foldchanges) <- res$entrez
keggres <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)

library(dplyr)

# Get the pathways
keggrespathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=5) %>%
  .$id %>%
  as.character()
keggrespathways

# Get the IDs.
keggresids <- substr(keggrespathways, start=1, stop=8)
keggresids

# Define plotting function for applying later
plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# Unload dplyr since it conflicts with the next line
detach("package:dplyr", unload=T)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))


