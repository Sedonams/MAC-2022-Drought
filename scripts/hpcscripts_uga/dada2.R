# PBJ Edits Aug2023
# Added 'pool = TRUE' condition to the dada2() step
# This is useful if you want to report alpha diversity estimates- without pooling ALL single-sample singletons will be removed rather than whole-dataset singletons
# In the 'mergePairs' function I changed the approach to merge paired ends instead of concatenating them

# PBJ Edits July 2024
# Changed how filtered and trimmed file names are given... I am removing anything that changes the base name because they are already in a different directory

####### Load required packages #######

suppressMessages(library(dada2)); print(paste0("loading dada 2 version ", packageVersion("dada2")))
suppressMessages(library(ShortRead)); print(paste0("loading ShortRead version ", packageVersion("ShortRead")))
suppressMessages(library(Biostrings)); print(paste0("loading Biostrings 2 version ", packageVersion("Biostrings")))
suppressMessages(library(stringr)); print(paste0("loading stringr version ", packageVersion("stringr")))
library(dplyr)

####### File parsing #######

# set working directory:
wdPath <- getwd()
list.files(wdPath)

print("Defining dummy variable inputs...")

# Define the dataset dummy variables
dataset.value <- "dataset.value.input"
print(paste0("dataset.value = ", dataset.value))
run.value <- "run.value.input"
print(paste0("run.value = ", run.value))
rc.value <- "rc.value.input"
print(paste0("rc.value = ", rc.value))
concat.value <- "concat.value.input"
print(paste0("concat.value = ", concat.value))
pool.value <- "pool.value.input"
print(paste0("pool.value = ", pool.value))
outdir.value <- "outdir.value.input"
print(paste0("outdir.value = ", outdir.value))
R1trunclen <- R1trunclen.value.input
print(paste0("R1trunclen = ", R1trunclen))
R2trunclen <- R2trunclen.value.input
print(paste0("R2trunclen = ", R2trunclen))
minlen <- minlen.value.input
print(paste0("minlen = ", minlen))


if (concat.value == "true"){

    print("dada2 denoised reads will be concatenated, not merged")

} else {

    print("dada2 denoised reads will be merged, not concatenated")

}

print("Organizing file locations and sample names...")

# specify subdirectories
trimpath <- paste0(wdPath, paste0("/trimmed/",dataset.value, run.value, rc.value)) # directory containing fwd and rev primer-trimmed reads
print('Trimmed input files are named  \n')
list.files(trimpath)[1:20] # Make sure the trimmed files are in there

filtpath <- paste0(wdPath, paste0("/filtered/",dataset.value, run.value)) # directory containing fwd and rev quality-filtered reads
print('Filtered files are named \n')
list.files(filtpath)[1:20] # Should be empty at beginning of script

# This creates a list of file paths for the primer-trimmed input files
# It assumes forward and reverse fastq filenames have format:
# SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fastqFs.trimmed <- sort(list.files(trimpath, pattern="_R1.fastq.gz", full.names = TRUE)) # get forward reads (R1)
fastqRs.trimmed <- sort(list.files(trimpath, pattern="_R2.fastq.gz", full.names = TRUE)) # get reverse reads (R2)
if(length(fastqFs.trimmed) != length(fastqRs.trimmed)) stop("Forward and reverse files do not match.")
print('fastqFs.trimmed are in:\n')
list(fastqFs.trimmed[1:5])
list(fastqFs.trimmed[1:5])

# Rewrite file names for filtered reads:
# It replaces 'trimmed' with 'filtered' in the full path (subdirectory name) without changing the base file name
fastqFs.filtered <- fastqFs.trimmed
fastqRs.filtered <- fastqRs.trimmed
# Replace the file path name
fastqFs.filtered <- str_replace_all(fastqFs.filtered,'trimmed/','filtered/')
fastqRs.filtered <- str_replace_all(fastqRs.filtered,'trimmed/','filtered/')
print('fastqFs.filtered are in:\n')
list(fastqFs.filtered[1:5])
list(fastqFs.filtered[1:5])

# Extract simplified sample names
sample.names <- str_remove_all(basename(fastqFs.filtered),'_R1.fastq.gz')
names(fastqFs.filtered) <- sample.names
names(fastqRs.filtered) <- sample.names
print('Sample names are:\n')
print(sample.names[1:20])

####### Quality filtering + read truncation #######

print("Beginning pre-dada2 sequence filtering...")

# We have set out trunclen and minlen values earlier in the script

filterOutput <- filterAndTrim(fwd=fastqFs.trimmed, filt=fastqFs.filtered,
              rev=fastqRs.trimmed, filt.rev=fastqRs.filtered,truncLen=c(R1trunclen,R2trunclen),
              compress=TRUE, verbose=TRUE, multithread=TRUE, minLen = minlen, truncQ = 2, maxEE = c(2, 4))

saveRDS(filterOutput ,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "filterOutput", dataset.value,  run.value, ".RDS")))

filterOutput


# Check which files actually exist
existing_files <- file.exists(fastqFs.filtered)
# Subset fastqFs.filtered to include only the existing files
fastqFs.filtered <- fastqFs.filtered[existing_files]
# Check which files actually exist
existing_files <- file.exists(fastqRs.filtered)
# Subset fastqRs.filtered to include only the existing files
fastqRs.filtered <- fastqRs.filtered[existing_files]
sample.names <- str_remove_all(basename(fastqFs.filtered),'_R1.fastq.gz')

print("Carrying out error rate determination...")

print("Read 1...")

####### Learn the Error Rates and denoise #######
errF <- learnErrors(fastqFs.filtered, multithread=TRUE)
saveRDS(errF,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "errF", dataset.value,  run.value, ".RDS"))) # save intermediate file

print(" Read 2...")

errR <- learnErrors(fastqRs.filtered, multithread=TRUE)
saveRDS(errR,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "errR", dataset.value,  run.value, ".RDS"))) # save intermediate file

print("Carrying out sequence dereplication...")

# Dereplicate:
derepFs <- derepFastq(fastqFs.filtered, verbose=TRUE)
saveRDS(derepFs,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "derepFs", dataset.value,  run.value, ".RDS"))) # save intermediate file
derepRs <- derepFastq(fastqRs.filtered, verbose=TRUE)
saveRDS(derepFs,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "derepRs", dataset.value,  run.value, ".RDS"))) # save intermediate file

# Rename derep-class objects:
names(derepFs) <- sample.names
names(derepRs) <- sample.names

print("Denoising dereplicated sequences...")

print("Read 1...")

# Apply the core sample inference algorithm to the dereplicated data.
if (pool.value == "true"){

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
saveRDS(dadaFs,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "dadaFs", dataset.value,  run.value, ".RDS"))) # save intermediate file

print("Read 2...")

dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)
saveRDS(dadaRs,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "dadaRs", dataset.value,  run.value, ".RDS"))) # save intermediate file

} else {

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = FALSE)
saveRDS(dadaFs,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "dadaFs", dataset.value,  run.value, ".RDS"))) # save intermediate file

print("Read 2...")

dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = FALSE)
saveRDS(dadaRs,file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "dadaRs", dataset.value,  run.value, ".RDS"))) # save intermediate file

}

print("Merging / concatenating read pairs...")


####### Merge paired-end reads by concatenating #######
# Merge the paired end reads, just concatenate
# instead of looking for overlapping regions.
# Will just stick both ends together with NNNNNNNNNN between.

if (concat.value == "true"){

    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 12, justConcatenate = TRUE, verbose=TRUE)
    saveRDS(mergers,file=paste0(wdPath, paste("/dada2output/", outdir.value,"/", "mergers", dataset.value,  run.value, ".RDS",sep="")))

} else {

   mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 12, justConcatenate = FALSE, verbose=TRUE)
   saveRDS(mergers,file=paste0(wdPath, paste("/dada2output/", outdir.value,"/", "mergers", dataset.value,  run.value, ".RDS",sep="")))

}


rm(derepFs,derepRs)   # clean up to save memory

print("Constructing sequence table...")

#######  Construct sequence table #######
seqtab <- makeSequenceTable(mergers)

print("Performing de novo chimera removal...")

#######  Remove Chimeras #######
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread=TRUE, verbose=TRUE
)

####### Export nochim sequence table for subsequent merging
saveRDS(seqtab.nochim,file=paste0(wdPath, paste("/dada2output/", outdir.value,"/", "ASVtab.nochim", dataset.value,  run.value, ".RDS",sep="")))

print("Creating sequence tracking table...")

####### Report the number of reads that made it through each step in the pipeline: #######

#getN <- function(x) sum(getUniques(x))

# Make a helper to turn any named vector into a data.frame
#make_df <- function(x, colname) {
#  data.frame(Sample = names(x), value = unname(x), stringsAsFactors = FALSE) |>
#    setNames(c("Sample", colname))
#}

#track <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), list(
#  make_df(filterOutput, "Filtered"),
#  make_df(sapply(dadaFs, getN), "DADA_F"),
#  make_df(sapply(dadaRs, getN), "DADA_R"),
#  make_df(sapply(mergers, getN), "Merged"),
#  make_df(rowSums(seqtab.nochim), "Non_chim")
#))

# Optional: replace NA with 0
#track[is.na(track)] <- 0

# See the result
#track


getN <- function(x) {
  sum(getUniques(x))
}

getN(dadaFs)
getN(mergers)
rowSums(seqtab.nochim)

#track <- cbind(
#  filterOutput,
#  sapply(dadaFs, getN),
#  sapply(dadaRs, getN),
#  sapply(mergers, getN),
#  rowSums(seqtab.nochim)
#)

#colnames(track) <- c(
#  "input",
#  "filtered",
#  "denoisedF",
#  "denoisedR",
#  "merged",
#  "nonchim")
#rownames(track) <- sample.names
#head(track)

####### write outputs #######

print("Writing final dada2 asv table and fasta file outputs...")

# Give sequence headers more manageable names (ASV_1, ASV_2...)
# Sike we're not doing that
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste(">",colnames(seqtab.nochim), sep = "") # No sensible name just ASV codes
#asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# Usually we would carry out this code to give such manageable names, but since we are combining seqtab tables together, we do NOT want to do this until$
# for (i in 1:dim(seqtab.nochim)[2]) {
# asv_headers[i] <- paste(">ASV", i, sep="_")
# }

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "ASVs",dataset.value,  run.value,".fasta")))

#count table
asv_tab <- t(seqtab.nochim)
#row.names(asv_tab) <- sub(">", "", asv_headers)
head(asv_tab)

#add placeholder column, then row.names=FALSE
Taxa <- rownames(asv_tab)
head(Taxa)
data <- cbind(Taxa,asv_tab)
head(data)
write.table(
  data,
  file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "ASVtable", dataset.value,  run.value, ".csv")),
              sep=",", quote=F, row.names = FALSE)

write.table(
  data,
  file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "ASVtable", dataset.value,  run.value, ".tsv")),
              sep="\t", quote=F, row.names = FALSE)

# Tracking table
#write.table(
#track,
#                file=paste0(wdPath, paste0("/dada2output/", outdir.value,"/", "tracking_table", dataset.value,  run.value, ".csv")),
#                            sep=",", quote=F, col.names=NA)

print("dada2.R complete")
