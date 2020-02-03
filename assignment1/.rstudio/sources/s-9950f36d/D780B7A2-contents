if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GEOmetadb", quietly = TRUE))
  BiocManager::install("GEOmetadb")

if (!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library("GEOmetadb")
library("edgeR")
library("biomaRt")

#METADATA -- NEEDED AT ALL?
gse <- getGEO("GSE77108",GSEMatrix=FALSE)
colnames(gse)
current_gpl <- names(GPLList(gse))[1]
current_gpl_info <- Meta(getGEO(current_gpl))
current_gpl_info$title
current_gpl_info$last_update_date
current_gpl_info$organism


#GET THE EXPRESSION DATA
sfiles = getGEOSuppFiles('GSE77108')
fnames = rownames(sfiles)
# there is only one supplemental file
diabetic.HDAC = read.delim(fnames[1],header=TRUE, check.names = FALSE)
#SUBSET ONLY COLUMNS WHICH have SAHA treatment
diabetic.HDAC.SAHA <- diabetic.HDAC[grep("SAHA|Geneid", colnames(diabetic.HDAC))]
#dim(diabetic.HDAC.SAHA)


samples <- data.frame(lapply(colnames(diabetic.HDAC.SAHA)[2:7], 
                             FUN=function(x){unlist(strsplit(x, split = "\\."))[c(2)]}))
colnames(samples) <- colnames(diabetic.HDAC.SAHA)[2:7]
rownames(samples) <- c("cell_type")
samples <- data.frame(t(samples), stringsAsFactors=FALSE)

#LABEL CELL TYPES AS "DIABETIC" OR "NORMAL"
diabetic <- grep("D", rownames(samples))
normal <- grep("N", rownames(samples))
samples$cell_type[diabetic] = "diabetic"
samples$cell_type[normal] = "normal"

#FILTER OUT GENES THAT HAVE LOW COUNTS
#use edgeR to calculate counts per million (cpms)
cpms = cpm(diabetic.HDAC.SAHA[2:7])
rownames(cpms) <- diabetic.HDAC.SAHA[,1]
#use cpms to determine which of our ENSGs to throw out
keep = rowSums(cpms >1) >=3
diabetic.HDAC.SAHA.filtered = diabetic.HDAC.SAHA[keep,]

#APPLYING TMM TO OUR DATASET
filtered_data_matrix <- as.matrix(diabetic.HDAC.SAHA.filtered[,2:7])
rownames(filtered_data_matrix) <- diabetic.HDAC.SAHA.filtered$Geneid
d = DGEList(counts=filtered_data_matrix, group=samples$cell_type)
d = calcNormFactors(d)
normalized_counts <- cpm(d)
#estimate common and tagwise dispersion, model_design not used since there is only 1 dimension (i.e. cell type)
d <- estimateDisp(d)

#get mart object from useast mirror (this mirror was working at the time of running)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

conversion_stash <- "SAHA_id_conversion.rds"
if(file.exists(conversion_stash)){
  HDAC_id_conversion <- readRDS(conversion_stash)
} else {
  HDAC_id_conversion <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                               filters = c("ensembl_gene_id"),
                               values = diabetic.HDAC.SAHA.filtered$Geneid,
                               mart = ensembl)
  saveRDS(HDAC_id_conversion, conversion_stash)
}

#number of ENSG IDs not mapped to HGNC symbols
nrow(diabetic.HDAC.SAHA.filtered) - nrow(HDAC_id_conversion)

