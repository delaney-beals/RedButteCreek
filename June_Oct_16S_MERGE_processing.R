#### 16S AMPLICON PROCESSING ####
getwd()
setwd("C:/Users/Puri Lab/Desktop/Delaney/June_Oct_merge/")
path <- "C:/Users/Puri Lab/Desktop/Delaney/June_Oct_merge"
list.files(path)

(.packages())
library(pacman)

p_load("dada2", "ShortRead", "Biostrings")

# Read in seq tables from runs 1, 2, 3, 4; Note: Sequence table format is sample (rows) x sequence variant (columns)
run1_bact_seq_table <- readRDS("C:/Users/Puri Lab/Desktop/Delaney/June_Oct_merge/bact_seq_table_June")

run2_bact_seq_table <- readRDS("C:/Users/Puri Lab/Desktop/Delaney/June_Oct_merge/bact_seq_table_Oct_rip_cDNA")

run3_bact_seq_table <- readRDS("C:/Users/Puri Lab/Desktop/Delaney/June_Oct_merge/bact_seq_table_Oct_rip_gDNA")

run4_bact_seq_table <- readRDS("C:/Users/Puri Lab/Desktop/Delaney/June_Oct_merge/bact_seq_table_Oct_up")

# Now merge all runs
bact_all_seq_tables <- mergeSequenceTables(run1_bact_seq_table, run2_bact_seq_table, run3_bact_seq_table, run4_bact_seq_table)


#### 11) Remove chimeras ####
seq_table_nochim <- removeBimeraDenovo(bact_all_seq_tables, method = "consensus", multithread = TRUE, verbose = TRUE)
# 23% identified as chimeric. 

# We don't know if these chimeras held a lot in terms of abundance. To find out, we can calculate the proportion of sequences retained after chimeras removed.
sum(seq_table_nochim)/sum(bact_all_seq_tables)
# We retained 95% of sequence abundance. Great!


#### 12) Track reads through the pipeline to verify everything worked as expected ####
getN <- function(x) sum(getUniques(x))
(summary_table <- data.frame(row.names = sample.names, 
                             input = out[, 1],
                             filtered = out[, 2], 
                             denoised_FWD = sapply(dada_FWD, getN),
                             denoised_REV = sapply(dada_REV, getN), 
                             merged = sapply(mergers, getN),
                             non_chim = rowSums(seq_table_nochim), 
                             final_perc_reads_retained = round(rowSums(seq_table_nochim)/out[, 1]*100, 1)))

# Looks good! The most amount of reads were removed during filtering, as expected. On average, retained 60% of reads. 


#### 13) Make summary tables that can be used to assign taxonomy

# giving seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seq_table_nochim)
asv_headers <- vector(dim(seq_table_nochim)[2], mode = "character")

for (i in 1:dim(seq_table_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# making a fasta file of the final ASV seqs. Once this is made, upload to rdp.cme.msu.edu/classifier/classifier.jsp to get taxonomy 
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "C:/Users/Puri Lab/Desktop/Delaney/June_Oct_merge/June_Oct_16S_MERGE.fa")

# make count table
asv_table <- t(seq_table_nochim)
row.names(asv_table) <- sub(">", "", asv_headers)
write.table(asv_table, "June_Oct_MERGE_ASV_counts.tsv", sep = "\t", quote = F, col.names = NA)


#### 14) Assign bacteria and archaea taxonomy with the Ribosomal Database Projects (RDP) training set database from https://zenodo.org/record/801828#.XNHerqZ7lSM ####
taxa <- assignTaxonomy(seq_table_nochim, "rdp_train_set_18.fa.gz", multithread = TRUE) # if filled with NAs; try including tryRC = TRUE
taxa <- addSpecies(taxa, "rdp_species_assignment_18.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# make taxa table
asv_taxa <- taxa
row.names(asv_taxa) <- sub(">", "", asv_headers)
write.table(asv_taxa, "June_Oct_MERGE_ASV_taxa_table.tsv", sep = "\t", quote = F, col.names = NA)
