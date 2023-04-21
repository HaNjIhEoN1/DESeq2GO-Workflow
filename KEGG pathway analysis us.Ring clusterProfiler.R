## ClusterProfile
BiocManager::install("clusterProfiler")
library(clusterProfiler)

hawKegg <- search_kegg_organism('Helicoverpa armigera', by = 'scientific_name') # 250 haw is helicoverpa armigera
dim(hawKegg)
hawKegg



## enrichKEGG() find overrepresented KEGG categories (pathways) in a set of significantly regulated genes
## It presents a representations of most biological metabolic or signalling transduction pathways, the genes are conducted through reactions
res_reverse_treatment <- res_reverse_treatment %>%
  rownames_to_column()

res_reverse_unfilTreatment <- res_reverse_unfilTreatment %>%
  rownames_to_column()

res_reverse_unfilPhe <- res_reverse_unfilPhe %>%
  rownames_to_column()

res_reverse_treatment$rowname <- sub("LOC", "", res_reverse_treatment$rowname) # remove LOC
res_reverse_unfilTreatment$rowname <- sub("LOC", "", res_reverse_unfilTreatment$rowname)

res_reverse_unfilPhe$rowname <- sub("LOC", "", res_reverse_unfilPhe$rowname)
res_reverse_phenotype$rowname <- sub("LOC", "", res_reverse_phenotype$rowname)

res_reverse_treatment_up$rowname <- sub("LOC", "", res_reverse_treatment_up$rowname)
res_reverse_treatment_down$rowname <- sub("LOC", "", res_reverse_treatment_down$rowname)

# create a vector for the gene universe
kegg_gene_list <- res_reverse_unfilTreatment$log2FoldChange
kegg_gene_list

kegg_gene_list_phe <- res_reverse_unfilPhe$log2FoldChange
# name vector
names(kegg_gene_list) <- res_reverse_unfilTreatment$rowname
kegg_gene_list

names(kegg_gene_list_phe) <- res_reverse_unfilPhe$rowname
kegg_gene_list_phe

# omis NAs
kegg_gene_list <- na.omit(kegg_gene_list)
kegg_gene_list_phe <- na.omit(kegg_gene_list_phe)

# sort the list in descreasing order need for clusterprofiler
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_gene_list_phe = sort(kegg_gene_list_phe, decreasing = T)

# extract significant result
kegg_genes <- res_reverse_treatment$log2FoldChange
kegg_gene_phe <- res_reverse_phenotype$log2FoldChange

kegg_treUp <- res_reverse_treatment_up$log2FoldChange
kegg_treDown <- res_reverse_treatment_down$log2FoldChange

# name the vector with ncbi ID
names(kegg_treUp) <- res_reverse_treatment_up$rowname 
names(kegg_treDown) <- res_reverse_treatment_down$rowname

# omit NA
kegg_genes <- na.omit(kegg_genes)
kegg_gene_phe <- na.omit(kegg_gene_phe)

kegg_treUp <- na.omit(kegg_treUp)
kegg_treDown <- na.omit(kegg_treDown)
kegg_treUp
kegg_treDown

# create enrichKEGG object
# use bitr_kegg to check if the genes are included in the database
kegg_organism = "haw"

kk <- enrichKEGG(gene = names(kegg_genes),
                 organism = "haw",
                 universe = names(kegg_gene_list),
                 pvalueCutoff = 0.05,
                 keyType = "kegg")
head(kk)

kk.tre.up <- enrichKEGG(gene = names(kegg_treUp),
                        organism = "haw",
                        universe = names(kegg_gene_list),
                        pvalueCutoff = 0.05,
                        keyType = "kegg")
kk.tre.down <- enrichKEGG(gene = names(kegg_treDown),
                        organism = "haw",
                        universe = names(kegg_gene_list),
                        pvalueCutoff = 0.05,
                        keyType = "kegg")
kk.tre.up.bar <- barplot(kk.tre.up,
                         title = "Enriched Pathways upregulation")
kk.tre.down.bar <- barplot(kk.tre.down,
                           showCategory = 20,
                         title = "Enriched Pathways downregulation")
bar.kk.tre.upANDdown <- ggarrange(kk.tre.up.bar,
                                  kk.tre.down.bar,
                                  labels = c("A", "B"),
                                  ncol = 1, nrow = 2,
                                  heights = c(0.3, 1.3))
                         
                         
kk.phe <- enrichKEGG(gene = names(kegg_gene_phe),
                     organism = "haw",
                     universe = names(kegg_gene_list_phe),
                     pvalueCutoff = 0.05,
                     keyType = "kegg")
head(kk.phe)

KEGGBar.treatment <- barplot(kk,
        showCategory = 20,
        title = "Enriched Pathways")

KEGGBar.phenotype <- barplot(kk.phe,
                             title = "Enriched Pathways")

dotplot(kk,
        showCategory = 13,
        title = "Enriched Pathways")

clusterProfiler::browseKEGG(kk,
                            'haw01200')

clusterProfiler::browseKEGG(kk,
                            'haw00010')



# cetegory net plot
BiocManager::install("ggnewscale")
library(ggnewscale)
cnetplot(kk, 
         categorySize = "pvalue",
         foldChange = kegg_gene_list)

# pathway view
BiocManager::install("pathview")
library(pathview)

## produce a native KEGG plot
Oxi.Pho<- pathview(gene.data = kegg_genes,
         pathway.id = "haw00190",
         species = "haw",
         gene.idtype = "KEGG")
