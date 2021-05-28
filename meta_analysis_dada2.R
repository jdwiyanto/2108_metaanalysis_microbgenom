path = ''
setwd(path)

library(dada2)
library(mbiome)
library(tidyverse)
library(phyloseq)
library(Biostrings)

# function 1: dada2 error learning and inference
dada2_step1 = function(batch, trunc) {
      library(dada2)
      library(mbiome)
      library(tidyverse)
      library(phyloseq)
      library(Biostrings)
  
      path2 = paste0(path, batch)
      dir.create(paste0(path2, '/out/'))
      
      # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
      fnFs <- sort(list.files(path2, pattern="_1.fastq", full.names = TRUE))
      fnRs <- sort(list.files(path2, pattern="_2.fastq", full.names = TRUE))
      
      # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
      sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
      
      # Place filtered files in filtered/ subdirectory
      filtFs <- file.path(path2, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
      filtRs <- file.path(path2, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
      
      names(filtFs) <- sample.names
      names(filtRs) <- sample.names
      
      out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen= trunc,
                           maxN=0, maxEE=c(2,6), rm.phix=TRUE,
                           compress=TRUE, multithread=F) # On Windows set multithread=FALSE
      
      # Learn error
      errF <- learnErrors(filtFs, multithread=TRUE)
      errR <- learnErrors(filtRs, multithread=TRUE)
      
      # Inference
      dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
      dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
      
      # Merged reads
      mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
      
      # Create sequence table
      seqtab <- makeSequenceTable(mergers)
      sequence_length = table(nchar(getSequences(seqtab)))
      write.table(sequence_length, file = paste0(path2, '/out/seqlength.tsv'), sep = '\t', quote = F, row.names = F)
      saveRDS(seqtab, file = paste0(path2, '/out/seqtab.rds'))
}

# function 2: compiling studies, removing chimera, and annotating sequences
dada2_step2 = function() {
  
  library(dada2)
  library(mbiome)
  library(tidyverse)
  library(phyloseq)
  library(Biostrings)
  
  seqtab.list = list.files(pattern = 'seqtab.rds', recursive = T, full.names = T)
  
  seqtab.list2 = lapply(seqtab.list, function(x) {
    df = readRDS(x)
    return(df)
  })
  
  seqtab = mergeSequenceTables(tables = seqtab.list2)
  seqtab.nochim = removeBimeraDenovo(seqtab, method = 'consensus', multithread = T, verbose = T)
  
  # set sequence length to keep
  seqlength_min = 248; seqlength_max = 260
  
  # filter sequences to expected length
  seqtab.nochim2 <- subset(seqtab.nochim, select = dplyr::between(nchar(getSequences(seqtab.nochim)), seqlength_min, seqlength_max))
  
  # Assign taxa
  taxa <- assignTaxonomy(seqtab.nochim2, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
  
  # combine abundance and taxonomy to ps object
  ps <- phyloseq::phyloseq(phyloseq::otu_table(seqtab.nochim2, taxa_are_rows=FALSE), phyloseq::tax_table(taxa))
  
  # merge ASV sequence to PS object
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- phyloseq::taxa_names(ps)
  ps <- phyloseq::merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(phyloseq::ntaxa(ps)), "_", phyloseq::tax_table(ps)[,6])
  
  # save PS object
  write.csv(phyloseq::tax_table(ps), 'phyloseq_tax.csv')
  write.csv(phyloseq::otu_table(ps), 'phyloseq_otu.csv')
  
  saveRDS(seqtab.nochim2, file = 'seqtab_nochim.rds')
  save.image(file = 'dada2_step2.RData')

}

# main

# list of truncLen criteria
trunc_00 = c('PRJEB34323', 'PRJNA513244', 'PRJNA574565', 'PRJNA525566', 'PRJNA553183', 'PRJNA338148', 'PRJNA353065', 'PRJNA631204')
trunc_140_130 = 'PRJEB20773'

# list of study batch
batches = list.dirs(full.names = F)
batches = batches[-1]

# running dada2_step1()

for (x in batches) {
  if (x %in% trunc_00) {
    dada2_step1(batch = x, trunc = c(0,0))
  } else if (x %in% trunc_140_130) {
    dada2_step1(batch = x, trunc = c(140,130))
    } else {
      dada2_step1(batch = x, trunc = c(150,150))
  }
}

# run dada2_step2