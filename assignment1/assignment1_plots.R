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

#METADATA
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
diabetic.HDAC = read.delim(fnames[1],header=TRUE,
                       check.names = FALSE)
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

#Distribution of our data - Boxplot
data2plot <- log2(cpm(diabetic.HDAC.SAHA.filtered[,2:7]))
boxplot(data2plot, xlab = "Samples", ylab = "log2 CPM", 
        las = 2, cex = 0.5, cex.lab = 0.5,
        cex.axis = 0.5, main = "normal and diabetic cells treated with SAHA")

#Distribution of our data - Density plot
counts_density <- apply(log2(cpm(diabetic.HDAC.SAHA.filtered[,2:7])), 2, density)
#calculate the limits across all the samples
xlim <- 0; ylim <- 0
for (i in 1:length(counts_density)) {
  xlim <- range(c(xlim, counts_density[[i]]$x)); 
  ylim <- range(c(ylim, counts_density[[i]]$y))
}
cols <- rainbow(length(counts_density))
ltys <- rep(1, length(counts_density))
#plot the first density plot to initialize the plot
plot(counts_density[[1]], xlim=xlim, ylim=ylim, type="n", 
     ylab="Smoothing density of log2-CPM", main="", cex.lab = 0.85)
#plot each line
for (i in 1:length(counts_density)) lines(counts_density[[i]], col=cols[i], lty=ltys[i])
#create legend
legend("topright", colnames(data2plot),  
       col=cols, lty=ltys, cex=0.75, 
       border ="blue",  text.col = "green4", 
       merge = TRUE, bg = "gray90")         


#APPLYING TMM TO OUR DATASET
filtered_data_matrix <- as.matrix(diabetic.HDAC.SAHA.filtered[,2:7])
rownames(filtered_data_matrix) <- diabetic.HDAC.SAHA.filtered$Geneid
d = DGEList(counts=filtered_data_matrix, group=samples$cell_type)
#d = DGEList(counts=filtered_data_matrix)
d = calcNormFactors(d)
normalized_counts <- cpm(d)

plotMDS(d, labels=rownames(samples),
        col = c("darkgreen","blue")[factor(samples$cell_type)])

#estimate common and tagwise dispersion
#model_design <- model.matrix(samples$cell_type+0)
#cell_types <- samples$cell_type
#model_design <- model.matrix(cell_types)
#d <- estimateDisp(d, model_design)
d <- estimateDisp(d)

#Graphing the BCV
plotBCV(d,col.tagwise = "black",col.common = "red")

#Create a visual representation of the mean-variance relationship
plotMeanVar(d, show.raw.vars = TRUE, 
            show.tagwise.vars=TRUE, 
            show.ave.raw.vars = TRUE,  
            NBline=TRUE,
            show.binned.common.disp.vars = TRUE)

#try(ensembl <- useMart("ensembl"))
#
ensembl <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl",
                   mirror = "useast")

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
