wd<- "/home/ies/jonathan.fiorentino/"
setwd(wd)

library("DESeq2")
library("scran")
library("BiocParallel")
register(MulticoreParam(4))

DE_NT_IVF <-function(sce,i){
  print(i)
  cl.lab <- (sce$seurat_clusters == i)
  
  dds <- convertTo(sce[,cl.lab], type="DESeq2")
  dds$exp <- sce[,cl.lab]$exp
  dds$isnt <- sce[,cl.lab]$isnt
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

res0 <- DE_NT_IVF(xenopus.combined.raw.sce,"0")
res1 <- DE_NT_IVF(xenopus.combined.raw.sce,"1")
res2 <- DE_NT_IVF(xenopus.combined.raw.sce,"2")
res3 <- DE_NT_IVF(xenopus.combined.raw.sce,"3")
res4 <- DE_NT_IVF(xenopus.combined.raw.sce,"4")
res5 <- DE_NT_IVF(xenopus.combined.raw.sce,"5")
res6 <- DE_NT_IVF(xenopus.combined.raw.sce,"6")
res7 <- DE_NT_IVF(xenopus.combined.raw.sce,"7")
res8 <- DE_NT_IVF(xenopus.combined.raw.sce,"8")
res9 <- DE_NT_IVF(xenopus.combined.raw.sce,"9")

write.csv(res0,"res0.csv", row.names = TRUE)
write.csv(res1,"res1.csv", row.names = TRUE)
write.csv(res2,"res2.csv", row.names = TRUE)
write.csv(res3,"res3.csv", row.names = TRUE)
write.csv(res4,"res4.csv", row.names = TRUE)
write.csv(res5,"res5.csv", row.names = TRUE)
write.csv(res6,"res6.csv", row.names = TRUE)
write.csv(res7,"res7.csv", row.names = TRUE)
write.csv(res8,"res8.csv", row.names = TRUE)
write.csv(res9,"res9.csv", row.names = TRUE)



