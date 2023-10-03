# conda activate r-mofa

library(Seurat)
library(SeuratDisk)


setwd('/storage/htc/joshilab/Su_Li/Alg_development/scRNA_TCRBCR_surfaceProtein/data_collection/anti-PD1_human_Tcell_NSCLC/')
raw = readRDS("GSE179994_all.Tcell.rawCounts.rds")

# Read the first column as a list from your table (replace 'your_table.csv' with your actual file path)
first_column_list <- read.csv('GSE179994_all.scTCR.tsv', header = TRUE,sep='\t')$CellName


intersection <- intersect(first_column_list, colnames(raw))


# Assuming your_matrix is a matrix with column names
subsetted_raw <- raw[, intersection]

library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
#[11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
#[16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
#[21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
#[26] "UNIPROT"  

ENSEMBL_id <- mapIds(org.Hs.eg.db, keys=rownames(subsetted_raw), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
head(ENSEMBL_id)

rownames(subsetted_raw) = ENSEMBL_id

seurat_object <- CreateSeuratObject(counts = subsetted_raw, 
                                    project = "NSCLC",
                                    min.cells = 1, min.features = 1)

# 17745 features across 77030 samples within 1 assay 


SaveH5Seurat(seurat_object, filename = "NSCLC_subsetted_raw.h5Seurat",overwrite = TRUE)
Convert("NSCLC_subsetted_raw.h5Seurat", dest = "h5ad", overwrite = TRUE)



