#Libraries used
library(ggplot2)
library(ggpubr)
library(edgeR)
library(reshape2)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db) #?
library(GSEABase) #?
library(clusterProfiler) #?

#Remove POB subjects from table of counts, cpm and cohort
ix4 <- which(colnames(CpmTblB) %in% c(POB.f0,POB.h0))
CpmTblC <- CpmTblB[, -ix4]
Counts3 <- Counts2[, -ix4]

ix5 <- which(rownames(cohort.female) %in% rownames(cohort.female[cohort.female$Type == "POB",]))
cohort.female2 <- cohort.female[-ix5,]

#Open subsets textfile and extract subset of interest
subset <- read.table('subset.txt',header=F,sep='\t')
S94 <- subset[100,2:27]  #Extract subset of interest
S94 <- as.matrix(S94)
S94 <- c(S94)
S94.OB <- cohort.female2[cohort.female2$Patnr %in% S94,]
S94.OBf <- S94.OB[S94.OB$newcond == 'OB.f0',]$colnames
S94.OBh <- S94.OB[S94.OB$newcond == 'OB.h0',]$colnames
S94Counts <- Counts3[, c(NO.f0,NO.h0,S94.OBf,S94.OBh)]
CpmS94Counts <- cpm(S94Counts)

#Make DGElists w groups and normalize
S94Groups <- S94.OB$newcond
NOGroups <- cohort.female2[cohort.female2$Type == 'NO',]$newcond
NOS94 <- c(NOGroups, S94Groups)
NOS94 <- factor(NOS94)
Groups2 <- cohort.female2$newcond    #?
Groups2 <- factor(Groups2)

load('AnnotationFile_allTC.Rdata')
annoTC <- which(annotation$TC %in% rownames(S94Counts))
y <- DGEList(counts = S94Counts, group = NOS94, genes = annotation[annoTC,]$gene)
z <- DGEList(counts = CpmS94Counts, group = NOS94, genes = annotation[annoTC,]$gene)
Counts_RLE <- calcNormFactors(y, method = "RLE")
TPM_RLE <- calcNormFactors(z, method = "RLE")


#Differential expression analysis
DEgroups <- NOS94
design <- model.matrix(~0+DEgroups, data = Counts_RLE$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(Counts_RLE, design) 
plotBCV(y)

fit <- glmQLFit(y, design)
my.contrasts <- makeContrasts(OBfvsNOf=OB.f0-NO.f0, OBhvsNOh=OB.h0-NO.h0, OBfvsOBh=OB.h0-OB.f0, NOfvsNOh=NO.h0-NO.f0, levels=design)

qlf1 <- glmQLFTest(fit, contrast=my.contrasts[,"OBfvsNOf"])
qlf2 <- glmQLFTest(fit, contrast=my.contrasts[,"OBhvsNOh"])
qlf3 <- glmQLFTest(fit, contrast=my.contrasts[,"OBfvsOBh"])
qlf4 <- glmQLFTest(fit, contrast=my.contrasts[,"NOfvsNOh"])
topTags(qlf4)

summary(decideTests(qlf4, p.value = 0.01))

opar <- par(no.readonly = TRUE)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 3))
## get yer plot on
plotMD(qlf3, cex.axis = 0.6, cex.lab = 0.6, hl.cex = 0.5, pch = 16, legend = F, col = c(rep("black",3), rep("red",3), rep("blue",3)))
legend(par("usr")[2], par("usr")[4], c("non-DE","Up","Down"), pch = 16, col = c("black","red","blue"), bty = "n", cex = 0.7)
par(xpd = F)
abline(h=c(-1, 1), col="black")
## set par back to original
#par(opar)


#------------------Testing-----Am I doing it right?

#Add Entrez gene ID to annotation
egSYMBOL2EG <- toTable(org.Hs.egSYMBOL2EG)
m <- match(y$genes$genes, egSYMBOL2EG$symbol)
y$genes$EZgene_id <- egSYMBOL2EG$gene_id[m]   #Add some gene_id...
qlf$genes$EZgene_id <- egSYMBOL2EG$gene_id[m]
qlf4$genes$EZgene_id <- egSYMBOL2EG$gene_id[m]

egREFSEQ <- toTable(org.Hs.egREFSEQ)   #Necessary?
n <- match(y$genes$EZgene_id, egREFSEQ$gene_id)
y$genes$REFSEQ <- egREFSEQ$accession[n]
qlf$genes$REFSEQ <- egREFSEQ$accession[n]

#---------------------end of Test

#GO and KEGG analysis, (require annotated genes)
got <- goana(qlf4$genes$EZgene_id, FDR = 0.05)
goTop <- topGO(got)

keg <- kegga(qlf4$genes$EZgene_id, FDR = 0.05)
kegTop <- topKEGG(keg)

