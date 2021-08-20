#This is the code for DEG analysis (Differential gene expression) with DESeq2.
#DESEQ2

#Read the metadata (either downloaded metadata from ArrayExpress database or send to you)
metadata <- read.table(file.path("/PATH_TO_FILE/SraRunTable.txt"), header = TRUE, stringsAsFactors=FALSE, sep = ",")
metadata
#Read the kallisto inputs (abundance.h5)
sample_id <- dir(file.path("/PATH_TO_FOLDER/kallisto/"))
sample_id
kal_files <- file.path("/PATH_TO_FOLDER/kallisto/", sample_id, "abundance.h5")
kal_files
#Vectorize transcripts to genes from ENSEMBL through bioMart (check most updated versions !-NOT ALWAYS WORKING!-)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = "feb2021.archive.ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description","transcript_biotype"),mart = mart)
head(t2g)
t2g$ensembl_transcript_id_version = paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep=".")
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id_version,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
t2g <- dplyr::select(t2g, c('target_id', 'ens_gene', 'ext_gene','transcript_biotype'))
head(t2g)
#Manipulate metadata to present a new dataframe with the desired information variables.
#Hushed comments refer to another analysis that took into account tumour glial content as an example. 
#The current unhashed one refer to ACP vs Fetal, ACP vs NFPA).
#Experiment with the dataframe to get a sense of the dataset manipulations
metadata <-dplyr::select(metadata, Run, disease)
metadata
metadata$disease[c(15,16,17)] = "Fetal"
metadata <- dplyr::rename(metadata, Type = disease)
metadata <- dplyr::rename(metadata, Run_ID = Run)
metadata
metadata$new_col <- NA
metadata <- dplyr::rename(metadata, Control = new_col)
#metadata <- dplyr::rename(metadata, Glial = new_col)
#for glial
#metadata$Glial[c(3,4,6,11,14,15,16,17,23,24)] = "No"
#metadata$Glial[c(1,2,5,10,12,18,21,22)] = "Yes"
metadata$Control[c(1,2,3,4,5,6,10,11,12,13,14,18,19,20,21,22,23,24)] = "No"
metadata$Control[c(7,8,9,15,16,17)] = "Yes"
metadata$Type[c(7,8,9)] = "NFPA"
metadata$Type[c(15,16,17)] = "Fetal"
#Import DESeq2 requiremnets to initiate the analysis
library('tximport')
library('DESeq2')
#[Optional: IF YOU WANT TO ANALYZE SPESIFIC KALLISTO FILES FROM THE DIRECTORY] 
#I will have to ommit some entries form kal_files cause the metadata and txi.kallisto need to have same number of rows [unhash on demand]
txi_kal_files<-kal_files#[c(1,2,3,4,5,6,10,11,12,14,15,16,17,18,21,22,23,24)]
txi_kal_files
#tximmport 
txi.kallisto <- tximport(txi_kal_files, type = "kallisto", tx2gene = t2g, ignoreAfterBar = TRUE,countsFromAbundance = "lengthScaledTPM",txOut = FALSE)
names(txi.kallisto)
head(txi.kallisto$counts)
#construct DESeq object
dds <- DESeqDataSetFromTximport(txi.kallisto, colData = ACP_metadata, design = ~Type)
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
dds <- DESeq(dds)

#Use the contrast command to isolate individual comparisons e.g., ACP vs Fetal and ACP vs NFPA.
#!!!THE ORDER IS IMPORTANT!! SEE LINE BELOW.
#contrast <- c("condition", "level_to_compare", "base_level")
res1 <- results(dds, contrast=c("Type", "Adamantinomatous Craniopharyngioma", "Fetal"))
res2 <- results(dds, contrast=c("Type", "Adamantinomatous Craniopharyngioma", "NFPA"))

#convert to df
res_data1 <- as.data.frame(res1)
res_data2 <- as.data.frame(res2)

#library(data.table)
#setDT(res_data1, keep.rownames = "ID")[]


#data manipulation to add gene annotation to transcripts in Fetal comparison.(lol is a temporary name)
#only unique entries in t2g were accounted for.
lol <- t2g[!duplicated(t2g[,3]),]
row.names(res_data1[,])
res_ens <- res_data1
res_ens$ens_gene <- rownames(res_ens)
x_fetal<- mergedTable <- merge(res_ens, lol, by = "ens_gene")
#get significant results from fetal
x_fetal
head(x_fetal)
sig_fetal <- subset(x_fetal, padj < 0.05, 
                    select=c(ens_gene, ext_gene,log2FoldChange,padj))

#for nfpa comparison
row.names(res_data2[,])
res_ens2 <- res_data2
res_ens2$ens_gene <- rownames(res_ens)
x_nfpa<- mergedTable <- merge(res_ens2, lol, by = "ens_gene")
x_nfpa
sig_nfpa <- subset(x_nfpa, padj < 0.05, 
                   select=c(ens_gene, ext_gene,log2FoldChange,padj))
sig_nfpa

#Ordering 

sig_fetal_ordered <- sig_fetal[order(sig_fetal$padj),]
sig_nfpa_ordered <- sig_nfpa[order(sig_nfpa$padj),]


# Getting downregulates and upregulated genes for each comparison (change the thresholds to whatever you want to)
sig_nfpa_ordered_downregulated <- subset(sig_nfpa_ordered,log2FoldChange < 0)
sig_nfpa_ordered_upregulated <- subset(sig_nfpa_ordered,log2FoldChange > 0)


sig_fetal_ordered_downregulated <- subset(sig_fetal_ordered,log2FoldChange < 0)
sig_fetal_ordered_upregulated <- subset(sig_fetal_ordered,log2FoldChange > 0)


#Now I will extract common ens_gene entries from both expressions from the significant results 
#for downregulated
common_downreg <- subset(sig_nfpa_ordered_downregulated,ens_gene %in% sig_fetal_ordered_downregulated$ens_gene)
common_downreg
common_upreg <- subset(sig_nfpa_ordered_upregulated,ens_gene %in% sig_fetal_ordered_upregulated$ens_gene)
common_upreg 

#save if you want
#write.table(common_upreg,file = "/Users/manthospitoulias/Desktop/common_upreg")

#log2foldchange is also known as effect size meaning the size of the deregulation. You can do something like:
effect_size_downreg <- common_downreg[order(common_downreg$log2FoldChange),]
effect_size_downreg

#I will import 404 splicing factors and check if any of those is present in our case 
#(file downloaded from a paper, see my thesis)
library(readxl)
X404_Splice_factors <- read_excel("PATH/404_Splice_factors.xlsx")
sf_list <- X404_Splice_factors
#some downregulated splicing factors
sf_downreg <- subset(common_downreg,ext_gene %in% sf_list$ext_gene)
sf_downreg
#some upregulated splicing factors
sf_upreg <- subset(common_upreg,ext_gene %in% sf_list$ext_gene)
sf_upreg
#BRIO - Manually added a dataframe after using BRIO and downloading the results
#essentially we are checking for commom entries
sf_upreg_brio <- subset(common_upreg,ext_gene %in% tab_enriched_motifs$Protein)
sf_upreg_brio
sf_downreg_brio <- subset(common_downreg,ext_gene %in% tab_enriched_motifs$Protein)
sf_downreg_brio
#save any result if you want. Here I am saving sf_upreg as a table.
#devtools::install_github("ropensci/writexl")
#library("writexl")
#write_xlsx(sf_downreg,"/Users/manthospitoulias/Desktop/downreg_splicing.xlsx")
#write_xlsx(sf_upreg,"/Users/manthospitoulias/Desktop/upreg_splicing.xlsx")

#Some usefull plotting. plotCounts only plots a single gene. We will combine plotCount with ggplot to create
#summary plots(e.g., count plots of the whole CPLX family)
library("ggplot2")
plotCounts(dds,"ENSG00000213578",intgroup = "Type")
cplx1<- plotCounts(dds, gene="ENSG00000168993", intgroup="Type", 
                   returnData=TRUE)

ggplot(cplx1, aes(x=Type, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


cplx2<- plotCounts(dds, gene="ENSG00000145920", intgroup="Type", 
                   returnData=TRUE)
library("ggplot2")
ggplot(cplx2, aes(x=Type, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

cplx3<- plotCounts(dds, gene="ENSG00000213578", intgroup="Type", 
                   returnData=TRUE)
library("ggplot2")
ggplot(cplx3, aes(x=Type, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

cplx4<- plotCounts(dds, gene="ENSG00000166569", intgroup="Type", 
                   returnData=TRUE)
library("ggplot2")
ggplot(cplx4, aes(x=Type, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#counts for cplx family
p <- ggplot() +
  geom_point(data = cplx1, position=position_jitter(w=0.1,h=0),aes(x=Type, y=count, colour ="CPLX1"), shape = 15) + 
  geom_point(data = cplx2, position=position_jitter(w=0.1,h=0),aes(x=Type, y=count, colour ="CPLX2"), shape =16) +
  geom_point(data = cplx3, position=position_jitter(w=0.1,h=0),aes(x=Type, y=count, colour ="CPLX3"), shape = 17) +
  geom_point(data = cplx4, position=position_jitter(w=0.1,h=0),aes(x=Type, y=count, colour ="CPLX4"), shape =18) +
  scale_y_log10(breaks=c(5,10,25,100,250,400,800)) +
  scale_colour_manual(name="Legend",
                      values=c(CPLX1="red", CPLX2="blue", CPLX3="green", CPLX4 = "purple"))
p

#Another usefull graph is the PCA plot. We use rld which are normalized and log transformed counts
plotPCA(rld, intgroup= "Type", ntop = 500)

