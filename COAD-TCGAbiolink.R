##This is an analysis pipeline for RNA-seq data using TCGAbiolinks R-Package
#Author: Neda Ghohabi
#Date: 10.5.2022

##Install Bioconductor
#if(!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
##Queries HT-Seq counts data from COAD project
COADquery <- GDCquery(project = "TCGA-COAD", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts", legacy = F,
                      experimental.strategy = "RNA-Seq") 

##Downloads the query above
GDCdownload(query = COADquery, method = "api",)
##Puts all of the data together into a summarized experiment object, adds in the clinical data as well.
COADprpr <- GDCprepare(query = COADquery, summarizedExperiment = T)

save(COADprpr, file = "COADprepare.RData")
f= load(file = "COADprepare.RData")

clin <- COADprpr@colData@listData
library(SummarizedExperiment)
clinic <- colData(COADprpr)
clinical <- data.frame(clinic)
clindata <- subset(clinical, select = -treatments)
clindata <- unlist(clindata)

write.table(clindata, file = "D:/TCGA/COAD/clindata.txt",
            quote = F, sep = "\t")
#Treatment data from clinical
treatment <- clinical$treatments
treatment <- unlist(subset(clinical, select = treatments))
treatments <- data.frame(treatment)
#getting FPKM-UQ  Data
Exdata <-COADprpr@assays@data$fpkm_uq_unstrand
rawdata <- assay(COADprpr)
write.table(rawdata, file = "D:/TCGA/COAD/rawdata.txt",
            quote = F, sep = "\t")

gene <- SummarizedExperiment::rowRanges(COADprpr)
gene <- data.frame(gene)
boxplot(rawdata[1:50, 50:100], outline = F)

raw <- read.delim(file = "D:/TCGA/COAD/rawdata.txt")

#Distinguish number of normal(N) and cancer(C) data based on their TCGA Barcodes(bar)
bar <- data.frame(colnames(raw))
colnames(bar) <- "barcode"
bar$substr <- substr(bar$barcode, start = 14, stop = 15)
bar <- bar[order(as.numeric(bar$substr), decreasing = T),]
C <- bar[bar$substr < 10,]
N <- bar[bar$substr> 10,]
rawdata <- raw[ , as.character(bar$barcode)]
gr <- c(rep("normal", 41), rep("cancer", 480))
gr <- factor(gr)
design <- model.matrix(~ 0 + gr)
colnames(design)<- levels(gr)
##preprocessing 
#Filtering and Normalisation with limma and edgeR
library(edgeR)
library(limma)
dge <- DGEList(rawdata)
keep <- filterByExpr(dge, design = design)
keep2 <- data.frame(keep)
filt <- dge[keep, ,keep.lib.size = F]
norm <- calcNormFactors(filt, method = "TMM")
head(norm)
head(filt)
v <- voom(norm, design = design, plot = T)
E <- v$E
boxplot(E[ , c(1:20, 50:70)])
ex <- 2^E
Ex <- log2(ex + 1)
boxplot(Ex[ , c(1:20, 50:70)])
write.table(Ex, file = "D:/TCGA/COAD/edgeR_limma/COADnormalized_edgeR_limma.txt",
            quote = F, sep = "\t")
Ex <- read.delim(file ="D:/TCGA/COAD/edgeR_limma/COADnormalized_edgeR_limma.txt" )

#TCGAbiolink Normalization and Filtering
library(TCGAbiolinks)
COADquery <- GDCquery(project = "TCGA-COAD", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts", legacy = F,
                      experimental.strategy = "RNA-Seq") 


COADprpr <- GDCprepare(query = COADquery, summarizedExperiment = T)  


COADprepro <- TCGAanalyze_Preprocessing(object = COADprpr, cor.cut = 0.6)

bar <- data.frame(colnames(COADprepro))
colnames(bar) <- "barcode"                                   
bar$substr <- substr(bar$barcode, 14,15)
bar <- bar[order(bar$substr, decreasing = T),]
preExpr <- COADprepro[ , as.character(bar$barcode)]
gr <- c(rep("normal", 41), rep("cancer", 480))
gr <- factor(gr)
design <- model.matrix(~ 0 + gr)
colnames(design)<- levels(gr)

library(EDASeq)

COADNorm <- TCGAanalyze_Normalization(tabDF = preExpr, 
                                      geneInfo = geneInfoHT, method = "gcContent" )

save(COADNorm, file = "COADNorm-TCGAbiolinks.RData")
f = load(file = "COADNorm-TCGAbiolinks.RData")
#filtering data up to 3 times
COADFilt <- TCGAanalyze_Filtering(COADNorm, method = "quantile",
                                  qnt.cut = 0.25)
filt2 <- TCGAanalyze_Filtering(COADFilt, method = "quantile",
                               qnt.cut = 0.25)
filt3 <- TCGAanalyze_Filtering(filt2,method = "quantile",
                               qnt.cut = 0.25)
write.table(filt3,file = "D:/TCGA/COAD/COADnorm-filt.TCGAbiolinks.txt",
            quote = F, sep = "\t")

library(limma)
v <- voom(filt3, design = design, plot = T)
E <- v$E
boxplot(E[ , c(1:20,400:420)] , outline = F)
write.table(E, file = "D:/TCGA/COAD/TCGAbiolinks.COADNormalized.txt",
            quote = F, sep = "\t")
##Differential Gene Expression Analaysis 
##DEG using edgeR and limma##

fit <- lmFit(Ex, design = design)
contrast <- makeContrasts(cancer-normal, levels = design)
fit2 <- contrasts.fit(fit = fit,contrasts = contrast)
fit2 <- eBayes(fit = fit2)  
diff <- topTable(fit = fit2, number = Inf, sort.by = "P" ) 
diffup <- diff[diff$logFC>2, ]
up <- diffup[diffup$adj.P.Val<0.01, ]

diffup <- diff[c(diff$logFC>2 & diff$adj.P.Val<0.01),]
diffup2 <- diff[c(diff$logFC>2 | diff$adj.P.Val<0.01),]

write.table(diff, file = "E:/Maryam/TCGA/prostate/pradDEG.limma.txt",
            quote = F, sep = "\t")
diff <- read.delim("E:/Maryam/TCGA/prostate/pradDEG.limma.txt")


data <- pradprpr@colData
data <- data.frame(data)
library(SummarizedExperiment)
gene <- SummarizedExperiment::rowRanges(pradprpr)
gene <- data.frame(gene)
symbol <- gene[gene$gene_id %in% as.character(rownames(diff)), ]
diff$GeneName <- symbol$gene_name


##DEG using TCGAbiolinks##

C <- bar[bar$substr < 10,]
N <- bar[bar$substr> 10,]

sampleNT <- TCGAquery_SampleTypes(barcode = N$barcode, typesample = "NT")
sampleTM_TP <- TCGAquery_SampleTypes(barcode = C$barcode,
                                     typesample = c("TP","TM"))
ex <- 2^E
ex <- log2(ex+1)

NTdata<- ex[,as.character(sampleNT)]                                                                  
TMTPdata <- ex[, as.character(sampleTM_TP)]

DEG <- TCGAanalyze_DEA(mat1 = TMTPdata, mat2 = NTdata, Cond1type = "Tumor",
                       Cond2type = "Normal", method = "glmLRT",
                       fdr.cut = 0.05, logFC.cut = 0)
degup<- DEG[DEG$FDR<0.05, ]

DEG2 <- TCGAanalyze_DEA(mat1 = TMTPdata, mat2 = NTdata, Cond1type = "Tumor",
                        Cond2type = "Normal", method = "glmLRT",
                        fdr.cut = 0.01, logFC.cut = 2)

l <- lm(logFC ~ FDR, data = DEG )
plot(l)
anova(l)

write.table(DEG, file = "E:/Maryam/TCGA/prostate/pradDEG.TCGAbiolinks.txt",
            quote = F , sep = "\t")
read.delim("E:/Maryam/TCGA/prostate/pradDEG.TCGAbiolinks")



###converting Gene Id (genecode version) to Ensemble Id and gene symbol
#BiocManager::install("biomaRt")
library(biomaRt)

Ex <- read.delim(file ="D:/TCGA/COAD/edgeR_limma/COADnormalized_edgeR_limma.txt" )

mart <- useMart(biomart = "ensembl", dataset ="hsapiens_gene_ensembl")
gene_ids_version <- rownames(Ex)
getBM(attributes = c('ensembl_gene_id_version',
                     'hgnc_symbol'),
      filters = 'ensembl_gene_id_version', 
      values = gene_ids_version,
      mart = mart)

library(stringr)
gene_ids <- str_replace(gene_ids_version,
                        pattern = ".[0-9]+$",
                        replacement = "")
ensembl <- getBM(attributes = c('ensembl_gene_id',
                                'hgnc_symbol','ensembl_gene_id_version'),
                 filters = 'ensembl_gene_id', 
                 values = gene_ids,
                 mart = mart)

write.table(ensembl, file = "D:/TCGA/COAD/COAD-ensembl-symbol-vesion.txt",
            quote = F, sep = "\t")

#Convert gene identifiers to Entrez ids
library(clusterProfiler)
library(org.Hs.eg.db)

entrez_convert <- clusterProfiler::bitr(rownames(Ex), 
                                        fromType = "ENSEMBL",
                                        toType= "ENTREZID", 
                                        OrgDb =org.Hs.eg.db) 
entrez_convert <- entrez_convert %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::filter(!duplicated(ENSEMBL) & !duplicated(ENTREZID))
Ex <- Ex[rownames(Ex) %in% entrez_convert$ENSEMBL,]
rownames(Ex) <- entrez_convert$ENTREZID
