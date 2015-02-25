# Usage: manually execuate on local R, since lonestar cannot install biomaRt successfully.
# Content: Execuate the biomaRt to extract gene annotations we need in ARIC aSPU project.
# Date:"Wed Jan 29 22:12:48 2014"
# Author: Yang Yang
##########################################################################################
library(biomaRt)
ensembl=useMart("ensembl")
ensembl_Ds = useDataset("hsapiens_gene_ensembl",mart=ensembl)
listFilters(ensembl_Ds)
genes_annotation <- getBM(attributes=c('hgnc_symbol', 'chromosome_name','start_position','end_position', 'band','strand'), mart = ensembl_Ds)
save(genes_annotation,file = "genes_annotation_genomewide.Rdata")