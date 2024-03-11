wd<- "/home/ies/jonathan.fiorentino/"
setwd(wd)

library("DESeq2")
library("scran")
library("BiocParallel")
register(MulticoreParam(4))

DE_IVF_NT <-function(sce){

  dds <- convertTo(sce, type="DESeq2")
  dds$exp <- sce$exp
  dds$isnt <- sce$isnt
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  print(dds)
  dds <- DESeqDataSet(dds, design = ~ exp + isnt)                     
  dds <- DESeq(dds,parallel=TRUE)
  res <- results(dds, contrast=c("condition","NT","IVF"),parallel=TRUE)
  res <- res[order(res$padj, decreasing = F),]
  res <- res[order(res$log2FoldChange, decreasing = T),]
  res <- res[!is.na(res$padj),]
  res <- res[res$padj< 0.1,]
  res <- as.data.frame.matrix(res)
  res
}

xenopus.combined.raw.sce <- readRDS('xenopus_combined_raw_sce.rds')

res.tot <- DE_IVF_NT(xenopus.combined.raw.sce)

write.csv(res.tot,"res_tot.csv", row.names = TRUE)


