# Produce a custom Fasta file from RNAseq data
library(customProDB)
library(devtools)
## install_github("mffrank/RNAseqToCustomFASTA")
## devtools::install("~/code/RNAseqToCustomFASTA")
# load_all(paste0(root_dir, "Master_Project/src/RNAseqToCustomFASTA"))
library(RNAseqToCustomFASTA)
# load_all("Y:/Master_Project/src/RNAseqToCustomFASTA")

#--------------------------
## Execution Parameters

root_dir <- "~/projects/PRPF8/"

## Download annotations

## ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
##                             dataset="hsapiens_gene_ensembl",
##                             ## host="http://jan2019.archive.ensembl.org",
##                             archive=FALSE)
## customProDB::PrepareAnnotationEnsembl(mart=ensembl,
##                                       annotation_path=paste0(root_dir, "data/customFasta/Annotation"),
##                                       splice_matrix=TRUE,
##                                       dbsnp=NULL,
##                                       transcript_ids=NULL,
##                                       COSMIC=FALSE)

annotation_path <- paste0(root_dir,"data/customFasta/Annotation/")
output_directory <-  paste0(root_dir,"data/customFasta/") # dir for output fastas
rpkm_directory <-  paste0(root_dir,"data/RNAseq/cufflinks/")
## vcf_directory <- paste0(root_dir,"Internship/PRPF8_SEC/data/RNAseq/Variant_calls/")
bed_directory <-  paste0(root_dir,"data/RNAseq/Junction_tables/")
sample_name <- "PRPF8"
wd <-  paste0(root_dir,"results/customFasta/") #dir for plots etc.
# wd <-  "~/Desktop/tmp/" #dir for plots etc.

## Download all nescessary annotations from ENSEMBL - only has to be carried out once
# downloadAnnotations("Y:/Master_Project/data/Human_genome/GRCh38/CustomProDB_annotation")

# Load the annotation objects into the working env.

load(paste(annotation_path, "/exon_anno.RData", sep = ""))
load(paste(annotation_path,"/proteinseq.RData", sep = ""))
load(paste(annotation_path,"/ids.RData", sep = ""))
load(paste(annotation_path, "/procodingseq.RData", sep = ""))
load(paste(annotation_path,"/splicemax.RData", sep = ""))
txdb <- AnnotationDbi::loadDb(paste(annotation_path,"/txdb.sqlite", sep = ""))

## Step one: Generate Fasta according to expression levels
# -------------------------------------------------------

# Load FPKM values into data frame
setwd(rpkm_directory)
files <- list.files(pattern = ".*isoforms.*")

fpkm <- list()
for(file in files){
  sample_name <- paste0(strsplit(file, "_")[[1]][1:2],collapse = "_")
  fpkm[[sample_name]] <- getCufflinkstranscriptFPKM(isoforms.fpkm_tracking_file_path = file)
  fpkm[[sample_name]] <- fpkm[[sample_name]][sort(names(fpkm[[sample_name]]))]
}
fpkm_table <- as.data.frame(do.call(cbind, fpkm))
## fpkm_table <- fpkm_table[-which(duplicated(rownames(fpkm_table))),]
## clfpkm_table <- getCufflinksFPKM(files = files, filetype = "isoforms.fpkm_tracking")
setwd(wd)
saveRDS(fpkm_table, paste0(wd,"PRPF8_fpkm_table.rda"))

# Plot the distribution of FPKMs to get a feeling for the cutoff to be applied
allFPKMs <- as.vector(t(fpkm_table))
quantiles1 <- mean(allFPKMs[allFPKMs > 0] <= 1)
quantiles05 <- mean(allFPKMs[allFPKMs > 0] <= 0.5)
quantiles01 <- mean(allFPKMs[allFPKMs > 0] <= 0.1)

FPKMsummaryReport(fpkm_table = fpkm_table, id_table = ids, c(0.6,0.5, 0.3, 0.2))

# Look for protein coding transcripts and assign Protein id

fpkm_table_prot <- getProteinCoding(fpkm_table = fpkm_table,id_table = ids,
                                    transcript_id_name = "tx_name", protein_id_name = "pro_name")

## library(data.table)
## fpkm_table_prot_pl <- data.table::melt(fpkm_table_prot, variable.name="Sample", value.name="FPKM")

## library(ggplot2)
## ggplot(fpkm_table_prot_pl) +
##   geom_density(aes(color=Sample, x=log10(FPKM)))
## hist(log10(fpkm[[2]][fpkm[[2]] > 0.001]), breaks=50)
## ## fpkm_table_prot_measured <- filterDetectedTranscripts(fpkm_table = fpkm_table_prot)
## # cutoff <- quantile(apply(fpkm_table_prot_measured,1, mean),0.25)
## cutoff <- 1
## fpkm_table_prot_filt <- fp

# Build the filtered Fasta library
setwd(output_directory)

OutputsharedPro(fpkm_table_prot[,names(fpkm_table_prot) != "protein_id"],
                cutoff=1, share_sample=3,
                outfile = "expressed_canonical_proteins1FPKM.fasta", ids = ids, proteinseq =proteinseq)
finishFasta("expressed_canonical_proteins1FPKM.fasta", outfile = "expressed_canonical_proteins1FPKMfinished.fasta")

## Step three: Generate Fasta with alternative Splicing Proteins
# -------------------------------------------------------
setwd(wd)

pdf("Remaining_Junctions_different_filters.pdf")
filterSpacePlot(covfilter_unique = c(3,5,7), share_num = c(2,3), bed_directory, splicemax, txdb, ids)
dev.off()


juncs <- getStarJunctions(bed_directory, splicemax = splicemax, txdb = txdb, ids = ids,
                          pattern = "\\.tab", share_num = 3, extend = 100, skip = 0,
                          covfilter_unique = 5, covfilter_multi = 0, overhang = "max_overhang")

saveRDS(juncs,"Junctions_cov5_share3.rda")
junction_type <- JunctionType(juncs, splicemax, txdb, ids)
saveRDS(junction_type,"Junction_type_cov5_share3.rda")


## The genome is not actually used if only the alternative splicing junctions are analysed
## Still needs to be supplied as a parameter
library("BSgenome.Hsapiens.UCSC.hg38")
setwd(output_directory)

## Only output novel alternative splicing junctions
OutputJunctions(juncs = junction_type, procodingseq = procodingseq, proteinseq = proteinseq,
                genome = Hsapiens, exons = exon, all_junc_peptides= F,
                outfile="Splice_Junction_containing_proteins.fasta")

## When including splice proteins we also want to include the reference isoform
library(seqinr)
splice_fasta <- read.fasta("Splice_Junction_containing_proteins.fasta", as.string = T, seqtype = "AA")
splice_proteins <- unique(gsub("_.*", "", names(splice_fasta)))

fpkm_table_splice <- fpkm_table_prot[rownames(fpkm_table_prot) %in% splice_proteins,]

OutputsharedPro(fpkm_table_splice[,names(fpkm_table_splice) != "protein_id"],
                cutoff=0, share_sample=0,
                outfile = "expressed_canonical_proteinsSpliceRef.fasta", ids = ids, proteinseq =proteinseq)

finishFasta("expressed_canonical_proteinsSpliceRef.fasta", outfile = "expressed_canonical_proteinsSpliceReffinished.fasta", removeStars = T, addcRAP = F, addiRT = F)

## Step four: Concatenate all fasta files
# -------------------------------------------------------

concatenateFastas(path = output_directory, pattern = "Splice_Junction_containing_proteins.fasta|expressed_canonical_proteins1FPKMfinished\\.fa|expressed_canonical_proteinsSpliceReffinished\\.fa",
                  outfile = "PRPF8_protgen_combined_1FPKM.fasta", addcRAP = F, addiRT = F, saveMappingTable = T)
