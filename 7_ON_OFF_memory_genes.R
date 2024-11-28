
###############USED TO MAKE PAPER FIGS AND LISTS OF GENES#################################

setwd("/Users/jonathan/Desktop/MUNCHEN/Project_Xenopus/PAPER_DRAFT_JAN_2024/ANALYSES_FEB_2024/EVA_script_MAY2024/")

getwd()

library(ggplot2)
library(tidyverse)
library(readxl)
library(Seurat)
library(ggplot2)
library(cowplot)
################################################
##################PREPARE DATA (DE ALL CLUSTERS TOGETHER) as used in Fig 2 A and gene sets for B and C ##################
################################################

AllData<- read.csv("AllData.csv") #same data as res_tot.csv but added columns from bulk rna seq (fold changes, pvaluesan d donor cell expression)
colnames(AllData)

#when visualizing the data, there was a problem with very few genes with very low pvalue- they were written as 0. after log10transformation number was infinite, so I removed those for graphical visualization (after talking to Gabriele and Anotnio)
padj<- AllData$padj
AllData$padjlog10 <- -log(padj,10) # log10 transform for visualization
padjnotINF<- subset(AllData, padj> 0) 

#add the donor mean
Donor_padjnotINF<-(padjnotINF[,10:12])
padjnotINF$Donormean = apply(Donor_padjnotINF,1,mean)
padjnotINF$Donormeanlog2 <- log(1+padjnotINF$Donormean,2)

####################################################################################
##################FILTER DATA (DE ALL CLUSTERS TOGETHER)#############################
####################################################################################
UPreg <- subset(padjnotINF, padj < 0.05 & log2FoldChange >= 0)
DOWNreg <- subset(padjnotINF, padj < 0.05 & log2FoldChange <= 0)

###### ON and OFF memory
colnames(padjnotINF)

#Filter ON memory just based on mean RPKM in donors: -> THEY GET A B EXTENSION
FCcutoffNTIVF=log2(2)
FCcutoffIVFDonor = log2(1)

OnMemB2FC<- subset(padjnotINF,  padj < 0.05 & log2FoldChange >= FCcutoffNTIVF & Donormean >1 & logFC.IVF_vs_donor <= -FCcutoffIVFDonor)
nrow(OnMemB2FC)
write.csv(OnMemB2FC$gene_names, file="scRNAseq_OnMemB2FC$gene_names.csv",row.names=F,quote=F) 

OffMemB2FC<- subset(padjnotINF,  padj < 0.05 & log2FoldChange <= -FCcutoffNTIVF & Donormean <=20, logFC.IVF_vs_donor >= FCcutoffIVFDonor)
nrow(OffMemB2FC)
write.csv(OffMemB2FC$gene_names, file="scRNAseq_OffMemB2FC$gene_names.csv",row.names=F,quote=F) 


####MA Plot

plot(padjnotINF$Donormeanlog2,padjnotINF$log2FoldChange , 
     main="",
     col="lightgrey",
     pch=16,
     cex=0.4)

points(UPreg$Donormeanlog2, UPreg$log2FoldChange,
       main="",
       col="#db868c",
       pch=16,
       cex=0.4)

points(DOWNreg$Donormeanlog2, DOWNreg$log2FoldChange,
       main="",
       col="#bab6cc",
       pch=16,
       cex=0.4)

points(OnMemB2FC$Donormeanlog2, OnMemB2FC$log2FoldChange,
       main="",
       col="#dd3f4b",
       pch=16,
       cex=1)

points(OffMemB2FC$Donormeanlog2, OffMemB2FC$log2FoldChange,
       main="",
       col="#301c7a",
       pch=16,
       cex=1)



####################################################################################
##################Prepare and Filter Data DE genes Individual Clusters as used in Fig 2 D and gene sets for E and F #############################
####################################################################################

#DE genes in IVF vs NT  Up in IVF has positive log FC>0



DE0a <- read.csv(file = 'res0.csv')
names(DE0a)[1] <- "gene_names"
DE0= merge(DE0a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE0)[3] <- "log2FoldChange"
names(DE0)[7] <- "padj"

DE1a <- read.csv(file = 'res1.csv')
names(DE1a)[1] <- "gene_names"
DE1= merge(DE1a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE1)[3] <- "log2FoldChange"
names(DE1)[7] <- "padj"

DE2a <- read.csv(file = 'res2.csv')
names(DE2a)[1] <- "gene_names"
DE2= merge(DE2a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE2)[3] <- "log2FoldChange"
names(DE2)[7] <- "padj"

DE3a<- read.csv(file = 'res3.csv')
names(DE3a)[1] <- "gene_names"
DE3= merge(DE3a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE3)[3] <- "log2FoldChange"
names(DE3)[7] <- "padj"

DE4a <- read.csv(file = 'res4.csv')
names(DE4a)[1] <- "gene_names"
DE4= merge(DE4a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE4)[3] <- "log2FoldChange"
names(DE4)[7] <- "padj"

DE5a <- read.csv(file = 'res5.csv')
names(DE5a)[1] <- "gene_names"
DE5= merge(DE5a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE5)[3] <- "log2FoldChange"
names(DE5)[7] <- "padj"

DE6a <- read.csv(file = 'res6.csv')
names(DE6a)[1] <- "gene_names"
DE6= merge(DE6a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE6)[3] <- "log2FoldChange"
names(DE6)[7] <- "padj"

DE7a<- read.csv(file = 'res7.csv')
names(DE7a)[1] <- "gene_names"
DE7= merge(DE7a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE7)[3] <- "log2FoldChange"
names(DE7)[7] <- "padj"

DE8a <- read.csv(file = 'res8.csv')
names(DE8a)[1] <- "gene_names"
DE8= merge(DE8a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE8)[3] <- "log2FoldChange"
names(DE8)[7] <- "padj"

DE9a <- read.csv(file = 'res9.csv')
names(DE9a)[1] <- "gene_names"
DE9= merge(DE9a,padjnotINF , by.x="gene_names", by.y= "gene_names", all.x = TRUE)
names(DE9)[3] <- "log2FoldChange"
names(DE9)[7] <- "padj"
colnames(DE0)



#FILTER FOR FC and padj
#UPin NT

FCcutoffNTIVF=log2(2)


UpNT_0 <- subset(DE0, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_0)
UpNT_1 <- subset(DE1, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_1)
UpNT_2 <- subset(DE2, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_2)
UpNT_3 <- subset(DE3, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_3)
UpNT_4 <- subset(DE4, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_4)
UpNT_5 <- subset(DE5, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_5)
UpNT_6 <- subset(DE6, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_6)
UpNT_7 <- subset(DE7, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_7)
UpNT_8 <- subset(DE8, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_8)
UpNT_9 <- subset(DE9, log2FoldChange > FCcutoffNTIVF & padj<0.05)
nrow(UpNT_9)


#FILTER FOR logFC greater than 1 and padj
#DOWN in NT

DownNT_0 <- subset(DE0, log2FoldChange < -FCcutoffNTIVF & padj<0.05)
nrow(DownNT_0)
DownNT_1 <- subset(DE1, log2FoldChange < -FCcutoffNTIVF & padj<0.05)
nrow(DownNT_1)
DownNT_2 <- subset(DE2, log2FoldChange < -FCcutoffNTIVF & padj<0.05)
nrow(DownNT_2)
DownNT_3 <- subset(DE3, log2FoldChange < -FCcutoffNTIVF & padj<0.05)
nrow(DownNT_3)
DownNT_4 <- subset(DE4, log2FoldChange < -FCcutoffNTIVF & padj<0.05)
nrow(DownNT_4)
DownNT_5 <- subset(DE5, log2FoldChange < -FCcutoffNTIVF & padj<0.05)
nrow(DownNT_5)
DownNT_6 <- subset(DE6, log2FoldChange < -FCcutoffNTIVF & padj<0.05)
nrow(DownNT_6)
DownNT_7 <- subset(DE7, log2FoldChange < -FCcutoffNTIVF & padj<0.05)
nrow(DownNT_7)
DownNT_8 <- subset(DE8, log2FoldChange < -FCcutoffNTIVF& padj<0.05)
nrow(DownNT_8)
DownNT_9 <- subset(DE9, log2FoldChange < -FCcutoffNTIVF& padj<0.05)
nrow(DownNT_9)

DEnumbersFC <- matrix( c(nrow(DownNT_1),
                         nrow(DownNT_2),
                         nrow(DownNT_4),
                         nrow(DownNT_5),
                         nrow(DownNT_6),
                         nrow(DownNT_7),
                         nrow(DownNT_3),
                         nrow(DownNT_8),
                         nrow(DownNT_9),
                         nrow(DownNT_0), 
                      nrow(UpNT_1),
                      nrow(UpNT_2),
                      nrow(UpNT_4),
                      nrow(UpNT_5),
                      nrow(UpNT_6),
                      nrow(UpNT_7),
                      nrow(UpNT_3),
                      nrow(UpNT_8),
                      nrow(UpNT_9),
                      nrow(UpNT_0)  
), 
nrow=2,ncol=10, byrow=TRUE,
dimnames=list(c("Down in NT genes (logFC<1, padj<0.05)","Up in NT genes (logFC<-1)"),
              c("Gob1","Gob2", "CGP","ANPN", "CNPB", "CEP","BSC3","BSC8", "Mixed","NNE") )
) 
rownames(DEnumbersFC)


barplot(DEnumbersFC, main="DE genes per cluster |logFC|>1 padj<0.05", 
        ylab="Number of Genes", beside=TRUE, axis.lty= 0,col = 1:nrow(DEnumbersFC),cex.names=1,space=c(0,0.5),
        legend.text = TRUE, 
        args.legend = list(x = "topleft",
                           cex=0.75,
                           box.lty=0,
                           bg=0,
                           bty= "n",
                           inset = c(0.1, 0)))

#############

#####################################################################################################################

#add donor ivf fdr and logfc filter

colnames(DownNT_0)


FCcutoffIVFDonor = log2(1)

#OFF
scOFF_cluster_0_FDR_FC <- subset(DownNT_0, FDR.IVF_vs_donor <= 0.05 & Donormean < 20 & logFC.IVF_vs_donor > FCcutoffIVFDonor )
nrow(scOFF_cluster_0_FDR_FC)
write.csv(scOFF_cluster_0_FDR_FC, file="scOFF_cluster_0_FDR_FC.csv",row.names=F,quote=F) 

scOFF_cluster_1_FDR_FC <- subset(DownNT_1, FDR.IVF_vs_donor <= 0.05 & Donormean < 20 & logFC.IVF_vs_donor > FCcutoffIVFDonor )
nrow(scOFF_cluster_1_FDR_FC)
write.csv(scOFF_cluster_1_FDR_FC, file="scOFF_cluster_1_FDR_FC.csv",row.names=F,quote=F) 

scOFF_cluster_2_FDR_FC <- subset(DownNT_2, FDR.IVF_vs_donor <= 0.05 & Donormean < 20 & logFC.IVF_vs_donor > FCcutoffIVFDonor )
nrow(scOFF_cluster_2_FDR_FC)
write.csv(scOFF_cluster_2_FDR_FC, file="scOFF_cluster_2_FDR_FC.csv",row.names=F,quote=F)

scOFF_cluster_3_FDR_FC <- subset(DownNT_3, FDR.IVF_vs_donor <= 0.05 & Donormean < 20 & logFC.IVF_vs_donor > FCcutoffIVFDonor )
nrow(scOFF_cluster_3_FDR_FC)
write.csv(scOFF_cluster_3_FDR_FC, file="scOFF_cluster_3_FDR_FC.csv",row.names=F,quote=F)

scOFF_cluster_4_FDR_FC <- subset(DownNT_4, Donormean < 20 & FDR.IVF_vs_donor< 0.05& logFC.IVF_vs_donor > FCcutoffIVFDonor)
nrow(scOFF_cluster_4_FDR_FC)
write.csv(scOFF_cluster_4_FDR_FC, file="scOFF_cluster_4_FDR_FC.csv",row.names=F,quote=F) 

scOFF_cluster_5_FDR_FC <- subset(DownNT_5, Donormean < 20 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor > FCcutoffIVFDonor)
nrow(scOFF_cluster_5_FDR_FC)
write.csv(scOFF_cluster_5_FDR_FC, file="scOFF_cluster_5_FDR_FC.csv",row.names=F,quote=F) 

scOFF_cluster_6_FDR_FC <- subset(DownNT_6, Donormean < 20 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor > FCcutoffIVFDonor)
nrow(scOFF_cluster_6_FDR_FC)
write.csv(scOFF_cluster_6_FDR_FC, file="scOFF_cluster_6_FDR_FC.csv",row.names=F,quote=F) 

scOFF_cluster_7_FDR_FC <- subset(DownNT_7, Donormean < 20 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor > FCcutoffIVFDonor)
nrow(scOFF_cluster_7_FDR_FC)
write.csv(scOFF_cluster_7_FDR_FC, file="scOFF_cluster_7_FDR_FC.csv",row.names=F,quote=F) 

scOFF_cluster_8_FDR_FC <- subset(DownNT_8, FDR.IVF_vs_donor <= 0.05 & Donormean < 20 & logFC.IVF_vs_donor > FCcutoffIVFDonor )
nrow(scOFF_cluster_8_FDR_FC)
write.csv(scOFF_cluster_8_FDR_FC, file="scOFF_cluster_8_FDR_FC.csv",row.names=F,quote=F) 

scOFF_cluster_9_FDR_FC <- subset(DownNT_9, Donormean < 20 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor > FCcutoffIVFDonor)
nrow(scOFF_cluster_9_FDR_FC)
write.csv(scOFF_cluster_9_FDR_FC, file="scOFF_cluster_9.csv_FDR_FC",row.names=F,quote=F) 


#ON
scON_cluster_0_FDR_FC <- subset(UpNT_0, Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor < -FCcutoffIVFDonor)
nrow(scON_cluster_0_FDR_FC)
write.csv(scON_cluster_0_FDR_FC, file="scON_cluster_0_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_1_FDR_FC <- subset(UpNT_1, Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor < -FCcutoffIVFDonor)
nrow(scON_cluster_1_FDR_FC)
write.csv(scON_cluster_1_FDR_FC, file="scON_cluster_1_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_2_FDR_FC <- subset(UpNT_2, Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor < -FCcutoffIVFDonor)
nrow(scON_cluster_2_FDR_FC)
write.csv(scON_cluster_2_FDR_FC, file="scON_cluster_2_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_3_FDR_FC <- subset(UpNT_3, Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor < -FCcutoffIVFDonor)
nrow(scON_cluster_3_FDR_FC)
write.csv(scON_cluster_3_FDR_FC, file="scON_cluster_3_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_4_FDR_FC <- subset(UpNT_4,  Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor < -FCcutoffIVFDonor)
nrow(scON_cluster_4_FDR_FC)
write.csv(scON_cluster_4_FDR_FC, file="scON_cluster_4_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_5_FDR_FC <- subset(UpNT_5, Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor < -FCcutoffIVFDonor)
nrow(scON_cluster_5_FDR_FC)
write.csv(scON_cluster_5_FDR_FC, file="scON_cluster_5_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_6_FDR_FC <- subset(UpNT_6,  Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor < -FCcutoffIVFDonor)
nrow(scON_cluster_6_FDR_FC)
write.csv(scON_cluster_6_FDR_FC, file="scON_cluster_6_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_7_FDR_FC <- subset(UpNT_7,  Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor < -FCcutoffIVFDonor)
nrow(scON_cluster_7_FDR_FC)
write.csv(scON_cluster_7_FDR_FC, file="scON_cluster_7_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_8_FDR_FC <- subset(UpNT_8,  Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor< -FCcutoffIVFDonor)
nrow(scON_cluster_8_FDR_FC)
write.csv(scON_cluster_8_FDR_FC, file="scON_cluster_8_FDR_FC.csv",row.names=F,quote=F) 

scON_cluster_9_FDR_FC <- subset(UpNT_9,  Donormean > 1 & FDR.IVF_vs_donor< 0.05 & logFC.IVF_vs_donor< -FCcutoffIVFDonor)
nrow(scON_cluster_9_FDR_FC)
write.csv(scON_cluster_9_FDR_FC, file="scON_cluster_9_FDR_FC.csv",row.names=F,quote=F) 


######


numbers_FDR_FC <- matrix( c(nrow(scOFF_cluster_1_FDR_FC),
                            nrow(scOFF_cluster_2_FDR_FC),
                            nrow(scOFF_cluster_4_FDR_FC),
                         nrow(scOFF_cluster_5_FDR_FC),
                         nrow(scOFF_cluster_6_FDR_FC),
                         nrow(scOFF_cluster_7_FDR_FC),
                         nrow(scOFF_cluster_3_FDR_FC),
                         nrow(scOFF_cluster_8_FDR_FC),
                         nrow(scOFF_cluster_9_FDR_FC), 
                         nrow(scOFF_cluster_0_FDR_FC),                          
                         nrow(scON_cluster_1_FDR_FC),
                         nrow(scON_cluster_2_FDR_FC),
                         nrow(scON_cluster_4_FDR_FC),      
                         nrow(scON_cluster_5_FDR_FC),
                         nrow(scON_cluster_6_FDR_FC),
                         nrow(scON_cluster_7_FDR_FC),
                         nrow(scON_cluster_3_FDR_FC),
                         nrow(scON_cluster_8_FDR_FC),
                         nrow(scON_cluster_9_FDR_FC),
                         nrow(scON_cluster_0_FDR_FC)
), 
nrow=2,ncol=10, byrow=TRUE,
dimnames=list(c("OFF","ON"),
              c("Gob1","Gob1","CGP","ANPN", "CNPB", "CEP","BSC3","BSC8", "Mixed", "NNE"))
) 
rownames(numbers_FDR_FC)



barplot(numbers_FDR_FC, main="Memory genes per cluster FDR _FC rpkm donor ", 
        ylab="Number of Genes", beside=TRUE, axis.lty= 0,col = 1:nrow(numbers_FDR_FC),cex.names=1,space=c(0,0.5),
        legend.text = TRUE, 
        args.legend = list(x = "topleft",
                           cex=0.75,
                           box.lty=0,
                           bg=0,
                           bty= "n",
                           inset = c(0.1, 0)))



######################################################
##################UPSET plot used in Figure 2 E and F###########################
######################################################
#install.packages("UpSetR")
library(UpSetR)
ONgenes0 <- as.vector(scON_cluster_0_FDR_FC$gene_names)
ONgenes1 <- as.vector(scON_cluster_1_FDR_FC$gene_names)
ONgenes2 <- as.vector(scON_cluster_2_FDR_FC$gene_names)
ONgenes3 <- as.vector(scON_cluster_3_FDR_FC$gene_names)
ONgenes4 <- as.vector(scON_cluster_4_FDR_FC$gene_names)
ONgenes5 <- as.vector(scON_cluster_5_FDR_FC$gene_names)
ONgenes6 <- as.vector(scON_cluster_6_FDR_FC$gene_names)
ONgenes7 <- as.vector(scON_cluster_7_FDR_FC$gene_names)
ONgenes8 <- as.vector(scON_cluster_8_FDR_FC$gene_names)
ONgenes9 <- as.vector(scON_cluster_9_FDR_FC$gene_names)

OFFgenes0 <- as.vector(scOFF_cluster_0_FDR_FC$gene_names)
OFFgenes1 <- as.vector(scOFF_cluster_1_FDR_FC$gene_names)
OFFgenes2 <- as.vector(scOFF_cluster_2_FDR_FC$gene_names)
OFFgenes3 <- as.vector(scOFF_cluster_3_FDR_FC$gene_names)
OFFgenes4 <- as.vector(scOFF_cluster_4_FDR_FC$gene_names)
OFFgenes5 <- as.vector(scOFF_cluster_5_FDR_FC$gene_names)
OFFgenes6 <- as.vector(scOFF_cluster_6_FDR_FC$gene_names)
OFFgenes7 <- as.vector(scOFF_cluster_7_FDR_FC$gene_names)
OFFgenes8 <- as.vector(scOFF_cluster_8_FDR_FC$gene_names)
OFFgenes9 <- as.vector(scOFF_cluster_9_FDR_FC$gene_names)



listInputON <- list(ONgenes1,ONgenes2,ONgenes4,ONgenes5,ONgenes6,ONgenes7,ONgenes3,ONgenes8, ONgenes9,ONgenes0)
names(listInputON) <- c('Goblet cells I','Goblet cells II','Cement gland primordium','Anterior neural plate border','Chordal neural plate border','Ciliated epidermal progenitors','Basal stem cells I','Basal stem cells II','Mixed states', 'Non neural ectoderm')

listInputOFF <- list(OFFgenes1,OFFgenes2,OFFgenes4,OFFgenes5,OFFgenes6,OFFgenes7,OFFgenes3,OFFgenes8,OFFgenes9,OFFgenes0)
names(listInputOFF) <- c('Goblet cells I','Goblet cells II','Cement gland primordium','Anterior neural plate border','Chordal neural plate border','Ciliated epidermal progenitors','Basal stem cells I','Basal stem cells II','Mixed states', 'Non neural ectoderm')


upset(fromList(listInputON), 
      order.by = "freq",
      nintersects = 20, 
      nsets = 10,
      sets =  c('Goblet cells I','Goblet cells II','Cement gland primordium','Anterior neural plate border','Chordal neural plate border','Ciliated epidermal progenitors','Basal stem cells I','Basal stem cells II','Mixed states', 'Non neural ectoderm'),
      decreasing = T, 
      group.by = "degree",
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1,
      )

upset(fromList(listInputOFF), 
      order.by = "freq",
      nintersects = 20, 
      nsets = 10,
      sets =  c('Goblet cells I','Goblet cells II','Cement gland primordium','Anterior neural plate border','Chordal neural plate border','Ciliated epidermal progenitors','Basal stem cells I','Basal stem cells II','Mixed states', 'Non neural ectoderm'),
      decreasing = T, 
      group.by = "degree",
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1,
)


##################################################################
############### UMAP PLOTS OF DATA ################################
##################################################################
##################################################################


############## BEGIN FUNCTIONS ####################################
setwd("/Users/Eva/Library/Mobile Documents/com~apple~CloudDocs/Documents/LAB/PROJECTS/scRNAseq/Xenopus_scRNAseq_Analyses/R_scRNAseqAnalyses_Eva")

library(Seurat)
library(ggplot2)
library(cowplot)

# The following functions take as input the Seurat object and a gene name and plots

#Load the data
xenopus.combined <- readRDS(file = "xenopus_combined_and_clustered.rds")

# Run the scaling of the data for each gene: it is needed for plotting the expression
# of multiple genes together in the UMAP
DefaultAssay(xenopus.combined)<-"RNA"
xenopus.combined <- ScaleData(xenopus.combined, verbose = FALSE)


################################## BEGIN FUNCTIONS ################################################
#################################################################################################

PlotGeneExpression <- function(Sobj,gname){
  #plot gene expression in UMAP:
  p1 <- FeaturePlot(Sobj, features = c(gname),cols =    c("lightgrey", "navy"), order = TRUE, pt.size = 2) 
  #p1 <- FeaturePlot(Sobj, features = c(gname),cols =    c("lightgrey", "navy"), order = TRUE,split.by = "isnt", pt.size = 2) 
  #violin plot of gene expression per cluster split by NT and IVF:
  p2 <- VlnPlot(xenopus.combined, features = c(gname),split.by = "isnt",pt.size = 0, split.plot = TRUE)
   
  #plot UMAP colored according to IVF/NT condition
  p3 <- DimPlot(Sobj, reduction = "umap", group.by = "isnt") 
  #plot cluster ID in UMAP: 
  p4 <- DimPlot(Sobj, reduction = "umap", label = TRUE,  pt.size = 5)
  #Arrange them together in a grid
  plot_grid(p1, p2, p3,p4)
  #plot_grid(p1)
  #plot_grid(p3)
}

# EXAMPLE OF USAGE  - Plotting expression of one gene
gname<- "mab21l3.S"
PlotGeneExpression(xenopus.combined,gname)

################
##################################################################################
#################################################################################
##################################################################################
# The following function takes as input the Seurat object and a list of gene names and 
# plots the averaged scaled expression of the group of genes in the UMAP (compared to the 
# UMAP colored according to cluster labels)
PlotMultipleGeneExpression <- function(Sobj,gnames){
  # Get mean expression of genes of interest per cell
  mean.exp <- colMeans(x = Sobj$RNA@scale.data[gnames,], na.rm = TRUE)
  
  # Add mean expression values in 'object@meta.data$gene.set.score'
  if (all(names(x = mean.exp) == rownames(x = Sobj@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'Sobj@meta.data':\n", 
        "adding gene set mean expression values in 'Sobj@meta.data$gene.set.mean'")
    Sobj@meta.data$gene.set.mean <- mean.exp
  }
  
  # Plot mean expression using Seurat::FeaturePlot()
  p1<- FeaturePlot(object = Sobj, features = "gene.set.mean", cols =  c("grey", "navy"), order = T, split.by = "isnt", pt.size = 2)& theme(legend.position = c(0.8,0.8))
 
  plot_grid(p1)
}

#### EXAMPLE OF USAGE  - Plotting expression of a list  of genes
my.list.genesa<-read.csv("scRNAseq_OnMemB2FC$gene_names.csv") ##dd3f4b
#my.list.genesb<-read.csv("scRNAseq_OffMemB2FC$gene_names.csv") ##301c7a
my.list.genes<- c(my.list.genesa)
#my.list.genes<- c(my.list.genesb)
my.gene.set<- c(my.list.genes$x)

#Plot their averaged (scaled) expression in the UMAP and compare it to the clusters
PlotMultipleGeneExpression(xenopus.combined,my.gene.set)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
# The following function takes as input the Seurat object and a list of gene names and 
# plots the averaged scaled expression of the group of genes in the UMAP (compared to the 
# UMAP colored according to cluster labels)
PlotMultipleGeneExpression <- function(Sobj,gnames){
  # Get mean expression of genes of interest per cell
  mean.exp <- colMeans(x = Sobj$RNA@scale.data[gnames,], na.rm = TRUE)
  
  # Add mean expression values in 'object@meta.data$gene.set.score'
  if (all(names(x = mean.exp) == rownames(x = Sobj@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'Sobj@meta.data':\n", 
        "adding gene set mean expression values in 'Sobj@meta.data$gene.set.mean'")
    Sobj@meta.data$gene.set.mean <- mean.exp
  }
  
  # Plot mean expression using Seurat::FeaturePlot()
  p1<- FeaturePlot(object = Sobj, features = "gene.set.mean", cols =  c("grey", "navy"), order = T, pt.size = 2)& theme(legend.position = c(0.8,0.8))
  
  plot_grid(p1)
}

#### EXAMPLE OF USAGE  - Plotting expression of a list  of genes
my.gene.set<- c("otog.L","otog.S")  # few names of genes in parallel 
#Plot their averaged (scaled) expression in the UMAP and compare it to the clusters
PlotMultipleGeneExpression(xenopus.combined,my.gene.set)

###########
