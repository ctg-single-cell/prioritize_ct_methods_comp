#!/usr/bin/env Rscript

# Process Allen Brain dataset according to EPIC
# Example Usage: Rscript process_scrna_epic_aba.R --matrix ${matrix} --meta ${meta} --project_name Allen_M1_2020 --out_rpkm Allen_M1_2020_epic_rpkm.csv --out_loc Allen_M1_2020_epic_loc.csv
# This script should run out of the box for Allen Brain Multiple cortical region (2019) and M1 (2020).
# For running Allen Brain MTG 2022 (or any other dataset), there needs to be some slight modification that I addressed below:
# NOTE 1: before running this script, you probably need to update the key word "sample_name" in order to get the column name of the column that has the cells. For eaxmple, in the MTG Allen Brain 2022, this is called specimen_name. However, in the multiple cortical areas Allen Brain 2019 and M1 Allen Brain 2020, it's sample_name. You can probably modify the script to take this as an argument but I just have not done so yet.
# NOTE 2: I copied the function process_scRNA from the R EPIC program to this script, You can also source this script. Just be aware that you need to correctly set the path of the file gene.noMHC.loc (line 30 here)
# Date: 14-03-2023
# Author: Tanya Phung (t.n.phung@vu.nl)
suppressPackageStartupMessages(library("argparse"))
library(data.table)
library(tibble)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--matrix", help ="Path to the matrix csv. Row: cell. Column: gene")
parser$add_argument("--meta", help = "Path to the meta file")
parser$add_argument("--project_name", help = "Name of the project")
parser$add_argument("--out_rpkm", help = "Path to the rpkm output")
parser$add_argument("--out_loc", help = "Path to the loc output")

args <- parser$parse_args()

# epic process_scRNA function (I'm including the function here for convenience but you could also source the script)
process_scRNA <- function(SeuratObject, meta_ct_col) {
  gene.loc = read.table(file = "/home/tphung/codes/EPIC/package/inst/extdata/gene.noMHC.loc", sep = "\t",
                        header = FALSE, stringsAsFactors = FALSE)
  
  grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice")
  grch37 = useDataset("hsapiens_gene_ensembl", mart=grch37)
  
  HGNC_to_ENSE=getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', "strand"),
                     filters = 'hgnc_symbol',
                     values = rownames(SeuratObject),
                     mart = grch37)
  
  dictionary = as.data.frame(HGNC_to_ENSE)
  dictionary <- dictionary[dictionary$ensembl_gene_id %in% gene.loc$V1, ]
  
  dictionary$chromosome_name <- as.numeric(dictionary$chromosome_name)
  o <- order(dictionary$chromosome_name, dictionary$start_position, dictionary$end_position)
  dictionary <- dictionary[o,]
  
  rawCount <- GetAssayData(object = SeuratObject, assay = "RNA", slot = "counts")
  keep.idx <- which(!is.na(match(rownames(rawCount), dictionary$hgnc_symbol)))
  rawCount <- rawCount[keep.idx, ]
  # rawCount <- rawCount[match(dictionary$hgnc_symbol, rownames(rawCount)), ]
  
  # RPKM
  gene.length <- (dictionary$end_position - dictionary$start_position)/1000
  library.size <- apply(rawCount, 2, sum)
  rpkm <- (10^6)*t(t(rawCount)/library.size)/gene.length
  cts <- sort(unique(SeuratObject@meta.data[[meta_ct_col]]))
  
  log2_rpkm_ct <- matrix(NA, nrow = nrow(rpkm), ncol = length(cts))
  rownames(log2_rpkm_ct) <- rownames(rpkm)
  colnames(log2_rpkm_ct) <- cts
  for (ct in cts) {
    cell.idx = rownames(SeuratObject@meta.data)[which(SeuratObject@meta.data[[meta_ct_col]] == ct)]
    log2_rpkm_ct[,ct] <- log2(1 + apply(rpkm[, cell.idx, drop = FALSE], 1, mean))
  }
  
  SeuratObject <- subset(SeuratObject, features = rownames(rpkm))
  
  ngene = 8000
  var.features <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = ngene)
  
  common.genes.ngene <- var.features@assays$RNA@var.features
  dictionary.ngene <- dictionary[dictionary$hgnc_symbol %in% common.genes.ngene, ]
  log2_rpkm_ct.ngene <- log2_rpkm_ct[rownames(log2_rpkm_ct) %in% common.genes.ngene, ]
  colnames(log2_rpkm_ct.ngene) <- gsub(" ", "_", colnames(log2_rpkm_ct.ngene))
  
  log2_rpkm_ct.ngene = cbind(log2_rpkm_ct.ngene, apply(log2_rpkm_ct.ngene, 1, mean))
  colnames(log2_rpkm_ct.ngene)[ncol(log2_rpkm_ct.ngene)] = "Average"
  rownames(log2_rpkm_ct.ngene) <- dictionary$ensembl_gene_id[match(rownames(log2_rpkm_ct.ngene), dictionary$hgnc_symbol)]
  
  log2_rpkm_ct.ngene = as.data.frame(log2_rpkm_ct.ngene)
  
  dictionary.ngene <- dictionary.ngene[,-1]
  colnames(dictionary.ngene) <- c("V1","V2","V3","V4","V5")
  
  return(list(scRNA.rpkm = log2_rpkm_ct.ngene, scRNA.loc = dictionary.ngene))
}


scrna_data = fread(args$matrix) #read in matrix
# sanity check
dim(scrna_data)
scrna_data[1:5, 1:5]

meta = fread(args$meta) #read in metadata
# sanity check
dim(meta)
meta[1:5, 1:5]

if (dim(meta)[1] != dim(scrna_data)[1]) {
  print("The number of rows in the metadata is not equal to the number of rows in the scRNAseq data. Will subset")
  # make a copy of the original
  meta_orig = meta
  meta = meta_orig %>% filter(sample_name %in% scrna_data$sample_name)
  }

# after subsetting sanity check
print("Checking the dimension of the meta dataframe after subsetting")
dim(meta)

# transpose
scrna_data_mx = t(scrna_data[,-1])

colnames(scrna_data_mx) = scrna_data$sample_name #assign column names #TODO: make the column name into a variable (in the MTG 2022 Allen brain dataset, the colname of the samples is called sample_name)
dim(scrna_data_mx)
scrna_data_mx[1:5, 1:5]

rownames(meta) = scrna_data$sample_name #assign rownames to meta dataframe #TODO: make the column name into a variable
dim(meta)

# loading libraries
library(Seurat)
library(Matrix)
library(biomaRt)

# convert matrix to sparse
print("Convert matrix to sparse")
scrna_data_mx_sparse <- as(object = scrna_data_mx, Class = "dgCMatrix")

# create a seurat object
print("Create Seurat Object")
scrna_seurat <- CreateSeuratObject(counts = scrna_data_mx_sparse, project = args$project, min.cells = 5, min.features = 200,  meta.data = meta)

# remove the cells that don't have annotations
scrna_seurat@meta.data = scrna_seurat@meta.data[scrna_seurat@meta.data$subclass_label!="", ]

# process scRNAseq according to EPIC
scrna_seurat_processed <- process_scRNA(SeuratObject = scrna_seurat, meta_ct_col = "subclass_label")

scRNA.rpkm <- scrna_seurat_processed$scRNA.rpkm
scRNA.loc <- scrna_seurat_processed$scRNA.loc

# convert the rowname into a column
scRNA.rpkm_updated = rownames_to_column(scRNA.rpkm, var="Gene")
scRNA.loc_updated = rownames_to_column(scRNA.loc, var="Index")

write.table(scRNA.rpkm_updated, args$out_rpkm, row.names = FALSE, sep = ",", quote = FALSE)
write.table(scRNA.loc_updated, args$out_loc, row.names = FALSE, sep = ",", quote = FALSE)
