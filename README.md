# DESeq2GO-Workflow
Workflow from DESeq2 to Gene ontology analysis using TopGO for non-model organisms
# install and load required packages 
install.packages("tidyverse")
install.packages("DESeq2")

library(DESeq2)
library(tidyverse)
library(vsn)
library(RColorBrewer)
library(pheatmap)
library(dplyr)

### Generate a metadata file from the gtf
## Need to read in the GTF that contains the product information
library(data.table)
gtf_file <- fread("/scratch1/fan048/RNA_seq_Fang_2022/featurecount/GCF_002156985.1_Harm_1.0_genomic_FIXED.gtf")
class(gtf_file)

## keep the rows where the gene identified (LOC ID) is and also the text describing that gene (this seems to be the exon or CDS rows)
## Remove duplicate rows so only 1 row is present for reach LOC ID
library(stringr)
gtf_file_edited <- gtf_file %>% 
  filter(V3 == "CDS") %>% ## filter out the columns contains CDS
  separate(col = V9, 
           into = c("1","2","3","4","5","6","7","8","9"),
           sep = ";") %>% ## split the string by delimiter ";"
  select(c(3,9,14,15,16)) ## selected the columns of interest

## select only the LOC ID using stringr
gtf_file_edited$`1` <- str_sub(string = gtf_file_edited$`1`, start = 10, end = 21)

## remove duplicated LOC IDs
gtf_file_edited <- distinct(gtf_file_edited,`1`, .keep_all = TRUE)

## use mutate to add new columns for the product information and protein_id
gtf_file_edited <- gtf_file_edited %>%
  rowwise() %>%
  mutate(products = list(str_subset(c_across(everything()), pattern = "(?<= product ).*"))) %>%
  mutate(protein_id = list(str_subset(c_across(everything()), pattern = "(?<=protein_id ).*")))

## cleaning the products and protein_id columns
gtf_file_edited$products <- str_extract(string = gtf_file_edited$products, pattern = "(?<= product ).*") # grab only the info after product
gtf_file_edited$protein_id <- str_extract(string = gtf_file_edited$protein_id, pattern = "(?<=protein_id ).*") # grab only the info after protein_id

## select the columns of interest and rename the column
gtf_product_proteinid <- select(gtf_file_edited, `1`, products, protein_id) %>%
  rename(LOC_id = `1`)

## export to directory
write.csv(gtf_product_proteinid, "/scratch1/fan048/Rstudio/DataSet/gtf_product_proteinid.csv")


# load featurecount raw output from Linux (from subread), otherwise featurecount could be done in R using Rsubread package from bioconductor
# change the column names
# assign the genes into row names
counts_reverse <- read_tsv("/scratch1/fan048/Rstudio/DEseq2/counts_reverse.txt", skip = 1, col_select = c(1, 7:15), col_names = TRUE) %>%
  rename(, 
         "Gene_id" = "Geneid", 
         "GRuntreat1" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/GRuntreat1._hisat.bam",
         "GRuntreat2" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/GRuntreat2._hisat.bam",
         "GRuntreat3" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/GRuntreat3._hisat.bam",
         "CADuntreat1" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/Caduntreat1._hisat.bam",
         "CADuntreat2" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/Caduntreat2._hisat.bam",
         "CADuntreat3" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/Caduntreat3._hisat.bam",
         "CADtreat1" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/Cadtreat1._hisat.bam",
         "CADtreat2" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/Cadtreat2._hisat.bam",
         "CADtreat3" = "/scratch1/fan048/RNA_seq_Fang_2022/Samtools/Cadtreat3._hisat.bam") %>%
  column_to_rownames("Gene_id")

counts_reverseq <- counts_reverse[, c(7, 8, 9, 4, 5, 6, 1, 2, 3 )]

# read the sample information file in as tsv, make sure the rownames in sample information matches the column names in the countfiles exactly as DEseq2 does not guess and match
metaData <- read_tsv("sample_info", col_names = TRUE) %>%
  column_to_rownames("strain")

# Prior proceeding, check with all(rownames(##sample_info) == colnames(##countfile))
rownames(metaData) == colnames(counts_reverse)

# Set the factor for the DESeq2
metaData$phenotype <- factor(metaData$phenotype)
metaData$treatment <- factor(metaData$treatment)

# To now construct the data object from the count matrix and the metadata table
# design formula test the effect of the second factor by controlling the effect of the first factor
# test the effect of treatment controlling the effect of phenotype
dds_reverse <- DESeqDataSetFromMatrix(countData = counts_reverse,
                                      colData = metaData,
                                      design = ~ phenotype + treatment)

## relevel the reference for phenotype and treatment
dds_reverse$phenotype = relevel(dds_reverse$phenotype, "sus")
dds_reverse$phenotype ## check the reference level

dds_reverse$treatment = relevel(dds_reverse$treatment, "untreated")
dds_reverse$treatment ## check the reference level

# running DESeq2
DESeq_dds_reverse <- DESeq(dds_reverse)

resultsNames(DESeq_dds_reverse) ## check the comparisons made in DESeq2 and derive biological meaningful comparisons

res_reverse_unfilPhe <- results(DESeq_dds_reverse, contrast = c("phenotype", "res", "sus")) %>%
  as.data.frame()

res_reverse_phenotype <- results(DESeq_dds_reverse, contrast = c("phenotype", "res", "sus")) %>% # output file difference between mutant and wildtype without treatment
  as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1)

# check HaCad gene
res_reverse_phenotype[c("LOC110380022"),]

gtf_product_proteinid$LOC_id <- sub("LOC", "", gtf_product_proteinid$LOC_id)

res_reverse_unfilTreatment <- results(DESeq_dds_reverse, contrast = c("treatment", "treated", "untreated")) %>%
  as.data.frame()

res_reverse_unfilTreatment <- left_join(res_reverse_unfilTreatment, gtf_product_proteinid, by = c("rowname" = "LOC_id"))

write.csv(res_reverse_unfilTreatment, "/scratch1/fan048/Rstudio/DataSet/res_reverse_unfilTreatment.csv")

res_reverse_treatment <- res_reverse_unfilTreatment %>% 
  as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1)

write.csv(res_reverse_treatment, "/scratch1/fan048/Rstudio/DataSet/res_reverse_treatment.csv")

# check HaCad gene
res_reverse_treatment[c("LOC110380022"),]

# subset for the upregulated DGE when log2fc > 0 at diff cutoff
res_reverse_treatment_up <- res_reverse_treatment %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange))
res_reverse_treatment_down <- res_reverse_treatment %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(abs(log2FoldChange)))

## Conduct Gene Ontology analysis
setwd("/scratch1/fan048/Rstudio/GO/")
BiocManager::install("topGO")
library(topGO)

## go_mapping contains two columns contains the mapping between genes and categories of interest. A list where the names are genes and each entry is a vector containing GO categories associated with that gene
## The gene_ontologies table was obtained from the matchTheLOCID script, map the LOCID to the HaOG ID which linked to GOterms
gene_ontologies <-
  cleanGOAdded %>%
  dplyr::select(SYMBOL, GOTerm) %>% 
  dplyr::distinct(SYMBOL, GOTerm)  %>% 
  dplyr::filter(SYMBOL != "NA")

## Check if any NA gene IDs under gene_ontologies SYMBOL column
table(is.na(gene_ontologies$SYMBOL))

## this makes a list with SYMBOL = LOCID and columns become list of vectors for GOTerm
gene_ontologies_wide <-
  pivot_wider(gene_ontologies, names_from = GOTerm, values_from = GOTerm) %>%
  unite(GO, c(2:5666), sep = ",", remove = TRUE, na.rm = TRUE) %>%
  as.data.frame()

class(gene_ontologies_wide) # make sure it's a dataframe
write.csv(gene_ontologies_wide, "/scratch1/fan048/Rstudio/GO/gene_ontologiesWide.txt")

## gene2go format gene_ID<TAB>GO_ID1, GO_ID2....
## object returned by readMappings is a named list of character vectors. 
## list names give the genes identifiers
## each element of the list is a character vector and contains the GO identifiers annotated to be specific gene
geneID2GO <- readMappings("/scratch1/fan048/Rstudio/GO/gene_ontologiesWide.txt")
geneID2GO
str(head(geneID2GO))
class(geneID2GO)

## Predefined list of interesting genes
## geneNames become the gene universe, representing set of all genes from the study
geneNames <- names(geneID2GO)
head(geneNames)
class(geneNames)


## THIS IS FOR TREAMTNET DGE UP REGULATION 
## get the list of genes of interest for up and down regulated genes
res_reverse_treatment_up <- res_reverse_treatment_up %>%
  rownames_to_column()

## get the list of genes of interest
myInterestingGenes_treatmentUp <- res_reverse_treatment_up$rowname
head(myInterestingGenes_treatmentUp)
class(myInterestingGenes_treatmentUp)

# create binary vector for topGO
geneList_treUp <- factor(as.integer(geneNames %in% myInterestingGenes_treatmentUp))
class(geneList_treUp)
table(geneList_treUp) ## check if the 1s (for expressed) matches up regulated genes, it's ok if it doesnt match as some genes do not have assigned GO

names(geneList_treUp) <- geneNames
str(geneList_treUp) ## apply the LOC names to the geneList_treUp

## geneList object is a named factor that indicates which genes are interesting and which are not
## It should be straightforward to compute such a named vector, where user has its own pre-defined list of interesting genes (DGE)

## Now, construct the topGOdata object
## To build the topGOdata object
## ontology : character string specifying the ontology of interest (BP, MF or CC)
## allGenes: named vector of type numeric or factor. The names attribute contains the genes identifiers. The genes listed in this object the gene universe
## annotationFUN: function which maps genes identifiers to GO terms. There are a couple of annotation function included in the package. 
## annFUN.gene2GO this function is used when annotations are provided as a gene-to-GOs mapping.
GOdataBP_treUp <- new("topGOdata",
              allGenes = geneList_treUp,
              gene2GO = geneID2GO,
              ontology = "BP",
              nodeSize = 5,
              annot = annFUN.gene2GO)

GOdataMF_treUp <- new("topGOdata",
                      allGenes = geneList_treUp,
                      gene2GO = geneID2GO,
                      ontology = "MF",
                      nodeSize = 5,
                      annot = annFUN.gene2GO)

GOdataCC_treUp <- new("topGOdata",
                      allGenes = geneList_treUp,
                      gene2GO = geneID2GO,
                      ontology = "CC",
                      nodeSize = 5,
                      annot = annFUN.gene2GO)

## One advantage of having the gene scores (or better genes measurements) as well as a way to define which are the interesting genes, in the topGOdata object is that one can apply various group testing procedure.
## see topGO Vignette section 3.5 Using the genes score
GOdataBP_treUp
GOdataMF_treUp
GOdataCC_treUp
## not all the genes that are provided by geneList can be annotated to the GO
## this can be seen by comparing the number of all available genes, the genes present in geneList, with the number of feasible genes
## We are therefore forced at this point to restrict the gene universe to the set of feasible genes for the rest of the analysis
## The summary on the GO graph shows the number of GO terms and the relations between them of the specified GO ontology

## running the enrichment tests
## Three types of enrichment tests are available: 1. tests based on gene counts, 2. tests based on gene scores or gene ranks, 3. tests based on gene expression
## run topGO with elimination test

runTest.treUp.elim <- runTest(GOdataBP_treUp, 
                              algorithm = "elim", 
                              statistic = "Fisher")
res.treUp.BP <- GenTable(GOdataBP_treUp, 
                          elimKS = runTest.treUp.elim,
                          orderBy = "elimKS", topNodes = 20) # the elimKS represents the pvalues for the individual GOterm
class(res.treUp.BP)

runTest.treUp.elim.MF <- runTest(GOdataMF_treUp, 
                              algorithm = "elim", 
                              statistic = "Fisher")
res.treUp.MF <- GenTable(GOdataMF_treUp, 
                         elimKS = runTest.treUp.elim.MF,
                         orderBy = "elimKS", topNodes = 20)

class(res.treUp.MF)
write.csv(res.treUp.MF, "/scratch1/fan048/Rstudio/GO/Treatment_up/res.treUp.MF.csv")

runTest.treUp.elim.CC <- runTest(GOdataCC_treUp, 
                                 algorithm = "elim", 
                                 statistic = "Fisher")
res.treUp.CC <- GenTable(GOdataCC_treUp, 
                         elimKS = runTest.treUp.elim.CC,
                         orderBy = "elimKS", topNodes = 20)

class(res.treUp.CC)
write.csv(res.treUp.CC, "/scratch1/fan048/Rstudio/GO/Treatment_up/res.treUp.CC.csv")

## For accessing the GO terms' p values from a topGOresult object, use score() function
p.res.treUp.BP <- score(runTest.treUp.elim)
p.res.treUp.BP <- as.data.frame(p.res.treUp.BP) %>%
  rownames_to_column()
class(p.res.treUp.BP)

p.res.treUp.MF <- score(runTest.treUp.elim.MF)
p.res.treUp.MF <- as.data.frame(p.res.treUp.MF) %>%
  rownames_to_column()
class(p.res.treUp.MF)

## left join the p-values to the GOterms
colnames(p.res.treUp.BP)[1] <- c("GO")
colnames(res.treUp.BP)[1] <- c("GO")
res.treUp.BP <-
  dplyr::left_join(res.treUp.BP, p.res.treUp.BP, by = "GO")

write.csv(res.treUp.BP, "/scratch1/fan048/Rstudio/GO/Treatment_up/res.treUp.BP.csv")


## THIS IS FOR TREAMTNET DGE DOWN REGULATION 
res_reverse_treatment_down <- res_reverse_treatment_down %>%
  rownames_to_column()

myInterestingGenes_treatmentDown <- res_reverse_treatment_down$rowname
head(myInterestingGenes_treatmentDown)
class(myInterestingGenes_treatmentDown)

geneList_treDown <- factor(as.integer(geneNames %in% myInterestingGenes_treatmentDown))
class(geneList_treDown)
table(geneList_treDown) ## check if the 1s (for expressed) matches Down regulated genes, it's ok if it doesnt match as some genes do not have assigned GO

names(geneList_treDown) <- geneNames
str(geneList_treDown) ## apply the LOC names to the geneList_treDown

GOdataBP_treDown <- new("topGOdata",
                        allGenes = geneList_treDown,
                        gene2GO = geneID2GO,
                        ontology = "BP",
                        nodeSize = 5,
                        annot = annFUN.gene2GO)

GOdataMF_treDown <- new("topGOdata",
                        allGenes = geneList_treDown,
                        gene2GO = geneID2GO,
                        ontology = "MF",
                        nodeSize = 5,
                        annot = annFUN.gene2GO)

GOdataCC_treDown <- new("topGOdata",
                        allGenes = geneList_treDown,
                        gene2GO = geneID2GO,
                        ontology = "CC",
                        nodeSize = 5,
                        annot = annFUN.gene2GO)

GOdataBP_treDown
GOdataMF_treDown
GOdataCC_treDown

## Run topGO
runTest.treDown.elim <- runTest(GOdataBP_treDown, 
                                algorithm = "elim", 
                                statistic = "Fisher")
res.treDown.BP <- GenTable(GOdataBP_treDown, 
                           elimKS = runTest.treDown.elim,
                           orderBy = "elimKS", topNodes = 20) # the elimKS represents the pvalues for the individual GOterm
class(res.treDown.BP)
write.csv(res.treDown.BP, "/scratch1/fan048/Rstudio/GO/Treatment_Down/res.treDown.BP.csv")


runTest.treDown.elim.MF <- runTest(GOdataMF_treDown, 
                                   algorithm = "elim", 
                                   statistic = "Fisher")
res.treDown.MF <- GenTable(GOdataMF_treDown, 
                           elimKS = runTest.treDown.elim.MF,
                           orderBy = "elimKS", topNodes = 20)


class(res.treDown.MF)
write.csv(res.treDown.MF, "/scratch1/fan048/Rstudio/GO/Treatment_Down/res.treDown.MF.csv")

runTest.treDown.elim.CC <- runTest(GOdataCC_treDown, 
                                   algorithm = "elim", 
                                   statistic = "Fisher")
res.treDown.CC <- GenTable(GOdataCC_treDown, 
                           elimKS = runTest.treDown.elim.CC,
                           orderBy = "elimKS", topNodes = 20)


class(res.treDown.CC)
write.csv(res.treDown.CC, "/scratch1/fan048/Rstudio/GO/Treatment_Down/res.treDown.CC.csv")


## THIS IS FOR GENOTYPE DGE UP REGULATION 
res_reverse_phenotype <- res_reverse_phenotype %>%
  rownames_to_column()
write.csv(res_reverse_phenotype, "/scratch1/fan048/Rstudio/DataSet/res_reverse_phenotype.csv")

# subset for the upregulated DGE when log2fc > 0 at diff cutoff
res_reverse_phenotype_up <- res_reverse_phenotype %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange))
res_reverse_phenotype_down <- res_reverse_phenotype %>%
  filter(log2FoldChange < 0) %>%
  arrange(desc(abs(log2FoldChange)))

## get the list of genes of interest for up regulated genes
myInterestingGenes_phenotypeUp <- res_reverse_phenotype_up$rowname
head(myInterestingGenes_phenotypeUp)
class(myInterestingGenes_phenotypeUp)

# create binary vector for topGO
geneList_pheUp <- factor(as.integer(geneNames %in% myInterestingGenes_phenotypeUp))
class(geneList_pheUp)
table(geneList_pheUp) ## check if the 1s (for expressed) matches up regulated genes, it's ok if it doesnt match as some genes do not have assigned GO

names(geneList_pheUp) <- geneNames
str(geneList_pheUp) ## apply the LOC names to the geneList_pheUp

## construct topGO object
GOdataBP_pheUp <- new("topGOdata",
                      allGenes = geneList_pheUp,
                      gene2GO = geneID2GO,
                      ontology = "BP",
                      nodeSize = 5,
                      annot = annFUN.gene2GO)

GOdataMF_pheUp <- new("topGOdata",
                      allGenes = geneList_pheUp,
                      gene2GO = geneID2GO,
                      ontology = "MF",
                      nodeSize = 5,
                      annot = annFUN.gene2GO)

GOdataCC_pheUp <- new("topGOdata",
                      allGenes = geneList_pheUp,
                      gene2GO = geneID2GO,
                      ontology = "CC",
                      nodeSize = 5,
                      annot = annFUN.gene2GO)

GOdataBP_pheUp
GOdataMF_pheUp
GOdataCC_pheUp

## runtest

runTest.pheUp.elim <- runTest(GOdataBP_pheUp, 
                              algorithm = "elim", 
                              statistic = "Fisher")
res.pheUp.BP <- GenTable(GOdataBP_pheUp, 
                         elimKS = runTest.pheUp.elim,
                         orderBy = "elimKS", topNodes = 20) # the elimKS represents the pvalues for the individual GOterm
class(res.pheUp.BP)
write.csv(res.pheUp.BP, "/scratch1/fan048/Rstudio/GO/Phenotype_up/res.pheUp.BP.csv")

runTest.pheUp.elim.MF <- runTest(GOdataMF_pheUp, 
                                 algorithm = "elim", 
                                 statistic = "Fisher")
res.pheUp.MF <- GenTable(GOdataMF_pheUp, 
                         elimKS = runTest.pheUp.elim.MF,
                         orderBy = "elimKS", topNodes = 20)
class(res.pheUp.MF)
write.csv(res.pheUp.MF, "/scratch1/fan048/Rstudio/GO/Phenotype_up/res.pheUp.MF.csv")

runTest.pheUp.elim.CC <- runTest(GOdataCC_pheUp, 
                                 algorithm = "elim", 
                                 statistic = "Fisher")
res.pheUp.CC <- GenTable(GOdataCC_pheUp, 
                         elimKS = runTest.pheUp.elim.CC,
                         orderBy = "elimKS", topNodes = 20)
write.csv(res.pheUp.CC, "/scratch1/fan048/Rstudio/GO/Phenotype_up/res.pheUp.CC.csv")

## ## THIS IS FOR GENOTYPE DGE DOWN REGULATION 
## get the list of genes of interest for up regulated genes

myInterestingGenes_phenotypeDown <- res_reverse_phenotype_down$rowname
head(myInterestingGenes_phenotypeDown)
class(myInterestingGenes_phenotypeDown)

# create binary vector for topGO
geneList_pheDown <- factor(as.integer(geneNames %in% myInterestingGenes_phenotypeDown))
class(geneList_pheDown)
table(geneList_pheDown) 

names(geneList_pheDown) <- geneNames
str(geneList_pheDown) 

## construct topGO object
GOdataBP_pheDown <- new("topGOdata",
                      allGenes = geneList_pheDown,
                      gene2GO = geneID2GO,
                      ontology = "BP",
                      nodeSize = 5,
                      annot = annFUN.gene2GO)

GOdataMF_pheDown <- new("topGOdata",
                      allGenes = geneList_pheDown,
                      gene2GO = geneID2GO,
                      ontology = "MF",
                      nodeSize = 5,
                      annot = annFUN.gene2GO)

GOdataCC_pheDown <- new("topGOdata",
                      allGenes = geneList_pheDown,
                      gene2GO = geneID2GO,
                      ontology = "CC",
                      nodeSize = 5,
                      annot = annFUN.gene2GO)

GOdataBP_pheDown
GOdataMF_pheDown
GOdataCC_pheDown

## runtest

runTest.pheDown.elim <- runTest(GOdataBP_pheDown, 
                              algorithm = "elim", 
                              statistic = "Fisher")

res.pheDown.BP <- GenTable(GOdataBP_pheDown, 
                         elimKS = runTest.pheDown.elim,
                         orderBy = "elimKS", topNodes = 20) # the elimKS represents the pvalues for the individual GOterm
class(res.pheDown.BP)

res.pheDown.BP <- res.pheDown.BP %>%
  dplyr::filter(as.numeric(elimKS) < 0.05)

runTest.pheDown.elim.MF <- runTest(GOdataMF_pheDown, 
                                 algorithm = "elim", 
                                 statistic = "Fisher")
res.pheDown.MF <- GenTable(GOdataMF_pheDown, 
                         elimKS = runTest.pheDown.elim.MF,
                         orderBy = "elimKS")
class(res.pheDown.MF)


runTest.pheDown.elim.CC<- runTest(GOdataCC_pheDown, 
                                 algorithm = "elim", 
                                 statistic = "Fisher")
res.pheDown.CC <- GenTable(GOdataCC_pheDown, 
                         elimKS = runTest.pheDown.elim.CC,
                         orderBy = "elimKS")


write.csv(res.pheDown.BP, "/scratch1/fan048/Rstudio/GO/Phenotype_down/res.pheDown.BP.csv")
write.csv(res.pheDown.MF, "/scratch1/fan048/Rstudio/GO/Phenotype_down/res.pheDown.MF.csv")
write.csv(res.pheDown.CC, "/scratch1/fan048/Rstudio/GO/Phenotype_down/res.pheDown.CC.csv")


  
