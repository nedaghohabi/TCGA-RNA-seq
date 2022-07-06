# TCGA-RNA-seq
# TCGA Colorectal Cancer RNA seq Data Analysis Pipeline

## software requirement:
 R version 4.2.0 
The pipeline currently will run three different analyses based on several different R packages 
## Packages:
TCGAbiolinks, EDAseq, Limma, EdgeR, BiomaRt, ClusterProfiler, org.Hs.eg.db, CMScaller

## Description:

Firstly, STAR Counts are retrieved by R-package TCGAbiolinks through the GDC Portal in May 2022. Then, FPKM-UQ (upper quantile) data are preprocessed (Normalization and Filtering) by Limma and EdgeR for further analysis. The normalized result should be checked by plotting Boxplot or Voom. 
Additionally, Corresponding  Clinical data was also retrieved by TCGAbiolinks.
Differential expression analysis can be done to get DEG (differential expressed genes).

Normalized RNAseq data have Genecode IDs ( Ensembl gene id version) as their row names, so for further analysis gene ids should convert to Ensembl ID or gene symbol (gene name). Here, BiomaRt were used to convert them to Ensembl IDs and HGNC symbols. It is important to note that some gene IDs do not have an Ensembl ID because of different gene mapping or they could be noble transcripts, also some Ensembl ids do not have HGNC symbol and name because some of them are pseudogenes or non-protein coding, so in this case, they are removed.
Then, Ensembl IDs were converted to Entrez IDs by ClusterProfiler and org.Hs.eg.db, which is needed as an input for CMSclassifier.

