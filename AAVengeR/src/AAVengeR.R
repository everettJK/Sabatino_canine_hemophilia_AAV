library(yaml)
library(stringdist)
library(ShortRead)
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(gintools)
library(lubridate)
options(stringsAsFactors = FALSE)


# Read the config file.
configFile <- commandArgs(trailingOnly = TRUE)
if(! file.exists(configFile)) stop('Error -- configuration file not found.')
config  <- read_yaml(configFile)

#config <- read_yaml('/home/everett/releases/Canine_hemophilia_AAV/AAVengeR/configs/Sabatino.config')
#config <- read_yaml('/home/everett/releases/Sabatino_canine_hemophilia_AAV/AAVengeR/configs/Sabatino.config')


# Add a function to check over config file and indentify conflicting parameters...
if(dir.exists(config$outputDir)) stop('Error -- output directory already exists.')
source(file.path(config$softwareDir, 'AAVengeR.lib.R'))


# Create required directory structure
dir.create(config$outputDir)
invisible(sapply(c('tmp', 'readIDs', 'readsRemoved', 'seqChunks', 'sampleReads', 'logs', 'fragReads'), 
                 function(x) dir.create(file.path(config$outputDir, x))))


# Record session info.
write(capture.output(sessionInfo()), file = file.path(config$outputDir, 'sessionInfo.txt'))


# Store the start time for the run.
# This time will be used by the log file to calculate time elapsed since start.
config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))


# Start log file
logFile <- file.path(config$outputDir, 'logs', 'log')
write(date(), file = logFile)


# Read in the sample data and add default column values if not present.
samples <- read_delim(config$sampleConfigFile, delim  = ',', col_names = TRUE, col_types = cols())


checkConfigFilePaths(samples$refGenomeBLATdb)
if('vectorSeqFile' %in% names(samples)) checkConfigFilePaths(samples$vectorSeqFile)


if(! 'subject' %in% names(samples))   samples$subject <- 'subject'
if(! 'replicate' %in% names(samples)) samples$replicate <- 1


# Create unique sample identifiers, eg. subject~sample~replicate
if(any(grepl('~|\\|', paste(samples$subject, samples$sample, samples$replicate)))) 
  stop('Error -- tildas (~) are reserved characters and can not be used in the subject, sample, or replicate sample configuration columns.')

samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)
if(any(duplicated(samples$uniqueSample))) stop('Error -- all subject, sample, replicate id combinations are not unique.')


# Reverse compliment the provided index1 barcodes if requested.
# Expand to incude RC for a second bar code.
if(config$indexReads.rc) samples$index1Seq <- as.character(reverseComplement(DNAStringSet(samples$index1Seq)))


# Chunk read files. 
cluster <- makeCluster(config$demultiplexing.CPUs)
clusterExport(cluster, c('config', 'samples'))

invisible(parLapply(cluster, 
                    list(c(config$index1ReadsFile, config$sequence.chunk.size, 'index1Reads', file.path(config$outputDir, 'seqChunks')),
                         c(config$breakReadsFile,  config$sequence.chunk.size, 'breakReads',  file.path(config$outputDir, 'seqChunks')),
                         c(config$virusReadsFile,  config$sequence.chunk.size, 'virusReads',  file.path(config$outputDir, 'seqChunks'))), 
                 function(x){
                   library(ShortRead)
                   source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
                   createSeqChunks(x[[1]], x[[2]], x[[3]], x[[4]])
                }))


# Reset run timer.
logMsg(config, 'Sample chunking completed.', logFile)



# Correct the index reads with golayCorrection() if the Golay correction option is set to TRUE in the configuration file.
# This function replaces the existing index reads files and relies on a python script packaged with this software.

if(config$correctGolayIndexReads){
  logMsg(config, 'Starting Gloay correction, resetting timer.', logFile)
  config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))
  
  invisible(parLapply(cluster, list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE, pattern = 'index'), function(x){
    library(ShortRead)
    library(dplyr)
    source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
    golayCorrection(x)
  }))
  
  
  logMsg(config, 'Gloay correction completed.', logFile)
}


logMsg(config, 'Starting sample chunk threads, resetting timer.', logFile)
config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))


invisible(parLapply(cluster, list.files(file.path(config$outputDir, 'seqChunks'), pattern = 'virusReads', full.names = TRUE), function(f){
  library(ShortRead)
  library(tidyverse)
  source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
  
  # Capture the chunk identifier.
  chunk.n <- unlist(str_match_all(f, '\\.(\\d+)$'))[2]
  
  # Create a chunk specific log file.
  logFile <-  file.path(config$outputDir, 'logs', paste0('seqChunk_', chunk.n, '.log'))
  
  # Read in sequence chunks.
  # The virus read file name is provided by f, here we switch out 'virusReads' to read in other files.
  virusReads  <- readFastq(f)
  breakReads  <- readFastq(sub('virusReads', 'breakReads',  f))
  index1Reads <- readFastq(sub('virusReads', 'index1Reads', f))
  
  # Trim reads by base call quality scores using a sliding window approach, trim lane info off of read ids.
  preTrimReadIDs <-  sub('\\s+.+$', '', as.character(virusReads@id))
  writeLines(preTrimReadIDs, file.path(config$outputDir, 'readIDs', paste0(chunk.n, '.readIDs')))
  
  virusReads <- trimTailw(virusReads, 2, config$sequence.trim.qualCode, 5)
  breakReads <- trimTailw(breakReads, 2, config$sequence.trim.qualCode, 5)
  
  # One or more read lists are depleted - return.
  if(length(virusReads) == 0 | length(breakReads) == 0){
    logMsg(config, paste0('Chunk ', chunk.n, ': No reads remaining after quality trimming.'), logFile)
    return()
  }

  # Convert reads to DNAstring sets in order to conserve memory and more easily access read sequences.
  index1Reads <- shortRead2DNAstringSet(index1Reads)
  virusReads  <- shortRead2DNAstringSet(virusReads)
  breakReads  <- shortRead2DNAstringSet(breakReads)
  
  # Sync becuase reads may of been lost in the previous base call quality filter.
  # Better solution for unpacking returned list?
  reads <- syncReads(index1Reads, virusReads, breakReads)
  index1Reads <- reads[[1]];  virusReads  <- reads[[2]];  breakReads  <- reads[[3]]
  
  # Report reads removed from quality trimming.
  ids <- preTrimReadIDs[! preTrimReadIDs %in% names(virusReads)]
  if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('qualTrim.', chunk.n)))
  
  if(length(virusReads) == 0){
    logMsg(config, paste0('Chunk ', chunk.n, ': No reads remaining after quality trimming and read synchronization.'), logFile)
    return()
  }
  
  # Loop through samples in sample data file.
  invisible(lapply(1:nrow(samples), function(r){
    r <- samples[r,]
    
    # Trim over-read sequences.
    virusReads <- trimOverReadSeq(virusReads, r$virusRead.overReadSeq)
    breakReads <- trimOverReadSeq(breakReads, r$breakRead.overReadSeq)
    
    # Ensure that reads are long enough for upcoming sample specific tests.
    preReadSizeCheck1 <- names(virusReads)
    virusReads <- virusReads[width(virusReads) >= (config$virusReads.minLTRseqLength + config$trimmedRead.minLength)]
    breakReads <- breakReads[width(breakReads) >= (nchar(r$breakReadLinkerSeq) + config$trimmedRead.minLength)]
    
    # Sync becuase reads may of been lost in the previous length filter.
    reads <- syncReads(index1Reads, virusReads, breakReads)
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]];  breakReads  <- reads[[3]]
    
    # Report reads removed from min read size filter1.
    ids <- preReadSizeCheck1[! preReadSizeCheck1 %in% names(virusReads)]
    if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('readSizeCheck1.', r$uniqueSample, '.', chunk.n)))
    
    if(length(virusReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after first read size filter.'), logFile)
      return()
    }
    
    # Create barcode demultiplexing vectors.
    v1 <- rep(TRUE, length(virusReads))
    if('index1Reads.maxMismatch' %in% names(config)){
      v1 <- vcountPattern(r$index1Seq, index1Reads, max.mismatch = config$index1Reads.maxMismatch) > 0
    }
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(breakReads))
    if('breakReads.linkerBarcode.maxMismatch' %in% names(config)){
      testSeq <- substr(r$breakReadLinkerSeq, r$breakReadLinkerBarcode.start, r$breakReadLinkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(breakReads, r$breakReadLinkerBarcode.start, r$breakReadLinkerBarcode.end), max.mismatch = config$breakReads.linkerBarcode.maxMismatch) > 0
    }
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- base::intersect(which(v1), which(v2))
    if(length(i) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': No reads demultiplexed for sample ', r$uniqueSample, '.'), logFile)
      return()
    }
    
    reads <- syncReads(index1Reads[i], virusReads[i], breakReads[i])
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
    if(length(virusReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after post demultiplexing read synchronization.'), logFile)
      return()
    }
    
    # Test for break point reads that align to the vector plasmid.
    if(config$filter.removeVectorReadPairs){
      vectorReadIDs <- getVectorReadIDs(breakReads, config, r$vectorSeqFile)
    
      if(length(vectorReadIDs) > 0){
        virusReadsVector <- virusReads[names(virusReads) %in% vectorReadIDs]
        breakReadsVector <- breakReads[names(breakReads) %in% vectorReadIDs]
    
        names(virusReadsVector) <- paste0(names(virusReadsVector), '|', r$uniqueSample)
        names(breakReadsVector) <- paste0(names(breakReadsVector), '|', r$uniqueSample)
      
        writeFasta(virusReadsVector, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.virusReadsVector.', chunk.n, '.fasta')))
        writeFasta(breakReadsVector, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.breakReadsVector.', chunk.n, '.fasta')))
      }
    
      # Filter break point reads against read ids returned by vectorReadIDs().
      breakReads <- breakReads[! names(breakReads) %in% vectorReadIDs]
    
      # Sync reads. 
      reads <- syncReads(index1Reads, virusReads, breakReads)
      index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
      if(length(virusReads) == 0){
        logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after vector read removal and read synchronization.'), logFile)
        return()
      }
    }

    # Test the start of virus reads.
    v1 <- rep(TRUE, length(virusReads))
    if('virusReads.startTest.maxMismatch' %in% names(config)){
      
      v1 <- Reduce('|', lapply(unlist(strsplit(r$virusLTRseq, ';')), function(x){ 
           testSeq <- substr(unlist(strsplit(x, ','))[2], 1,  config$virusReads.startTest.length)
           vcountPattern(testSeq, subseq(virusReads, 1, config$virusReads.startTest.length), max.mismatch = config$virusReads.startTest.maxMismatch) > 0
        }))
    }
    
    # Test the entire LTR sequence.
    v2 <- rep(TRUE, length(virusReads))
    if('virusReads.fullTest.maxMismatch' %in% names(config)){
      # v2 <- vcountPattern(r$virusLTRseq, subseq(virusReads, 1, nchar(r$virusLTRseq)), max.mismatch = config$virusReads.fullTest.maxMismatch) > 0
      
      v2 <- Reduce('|', lapply(unlist(strsplit(r$virusLTRseq, ';')), function(x){ 
        vcountPattern(unlist(strsplit(x, ','))[2], subseq(virusReads, 1, nchar(unlist(strsplit(x, ','))[2])), max.mismatch = config$virusReads.fullTest.maxMismatch) > 0
      }))
      
    }
    
    # Test break read common linker.
    v3 <- rep(TRUE, length(virusReads))
    if('breakReads.linkerCommon.maxMismatch' %in% names(config)){
      testSeq <- substr(r$breakReadLinkerSeq, r$breakReadLinkerCommon.start, r$breakReadLinkerCommon.end)
      v3 <- vcountPattern(testSeq, subseq(breakReads, r$breakReadLinkerCommon.start,r$breakReadLinkerCommon.end), max.mismatch = config$breakReads.linkerCommon.maxMismatch) > 0
    }
    
    # Test to see which reads pass the last three filters.
    i <- base::Reduce(base::intersect, list(which(v1), which(v2), which(v3)))
    if(length(i) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after user specified read tests and read synchronization.'), logFile)
      return()
    }
    
    # Subset and sync reads.
    reads <- syncReads(index1Reads[i], virusReads[i], breakReads[i])
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
  
    # Trim leading adapter sequences.
    if(config$virusReads.captureLTRseqs){
       o <- captureLTRseqs(virusReads, r$virusLTRseq)
       
       o$reads <- o$reads[match(o$LTRs$id, names(o$reads))]                      # Order the reads to match LTR table.
       o$reads <- subseq(o$reads, start = nchar(o$LTRs$LTRseq)+1)                # Remove matched LTR sequences from the reads.
       o$reads <- o$reads[width(o$reads) >= config$trimmedRead.minLength]        # Remove reads which are now too short.
       o$LTRs  <- subset(o$LTRs, id %in% names(o$reads))                         # Now that we have removed some reads, trim the LTR table.
       o$LTRs$readSeq <- as.character(o$reads[match(o$LTRs$id, names(o$reads))]) # Add the trimmed sequences to the LTR table so that we can latter find 
                                                                                 #  additional NTs between LTRs and genomic alignments.
       
       saveRDS(o$LTRs, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.LTRseqs.', chunk.n, '.rds')))
       
       ids <- names(virusReads)[! names(virusReads) %in% names(o$reads)]
       if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('LTRcapture.', chunk.n)))
       
       virusReads <- o$reads
       
       reads <- syncReads(index1Reads, virusReads, breakReads)
       index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
       if(length(virusReads) == 0){
         logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after LTR capture and trimming.'), logFile)
         return()
       }
    } else {
      # This route does not support multiple LTR sequences -- add check.
      s <- unlist(strsplit(r$virusLTRseq, ','))[2]
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') Trimming leading seq: ', s), logFile)
      virusReads <- trimLeadingSeq(virusReads, s)
    }
    
    # Capture random ids.
    randomIDs <- data.frame()
    if('breakReadLinkerBarcode.start' %in% names(r) & 'breakReadLinkerBarcode.end' %in% names(r)){
      randomIDs <- as.character(subseq(breakReads, r$breakReadLinkerRandomID.start, r$breakReadLinkerRandomID.end))
      randomIDs <- data.frame(randomSeqID = unname(randomIDs), readID = names(randomIDs))
    }
    
    # Remove leading linker from break point reads.
    breakReads <- trimLeadingSeq(breakReads, r$breakReadLinkerSeq)
    
    # Select reads which have the minimum read lengths post trimming.
    preFilterReads <- names(virusReads)
    virusReads <- virusReads[width(virusReads) >= config$trimmedRead.minLength]
    breakReads <- breakReads[width(breakReads) >= config$trimmedRead.minLength]
    
    # Report reads removed from min read size filter2.
    ids <- preFilterReads[! preFilterReads %in% names(virusReads)]
    if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('readSizeCheck2.', r$uniqueSample, '.', chunk.n)))
    
    # Sync reads after selecting for reads with min. length post-trimming.
    reads <- syncReads(index1Reads, virusReads, breakReads)
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
    if(length(virusReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after second read size filter.'), logFile)
      return()
    }
    
    # Write out final reads. Add sample names to read IDs.
    names(breakReads) <- paste0(names(breakReads), '|', r$uniqueSample)
    names(virusReads) <- paste0(names(virusReads), '|', r$uniqueSample)
    
    save(randomIDs, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.randomIDs.', chunk.n, '.RData')))
    writeFasta(breakReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.breakReads.', chunk.n, '.fasta')))
    writeFasta(virusReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.virusReads.', chunk.n, '.fasta')))
    
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') completed with ', length(virusReads), ' reads.'), logFile)
  }))
  
  logMsg(config, paste0('Chunk ', chunk.n, ' completed.'), logFile)
  logMsg(config, paste0('Chunk ', chunk.n, ' completed.'), file.path(config$outputDir, 'logs', 'log'))
}))

stopCluster(cluster)


# Collate demultiplexed reads using file name snibets.
collateSampleReads('virusReads')
collateSampleReads('breakReads')


# Collate read pairs where the break point read aligns well to the vector for downstream analysis.
if(config$filter.removeVectorReadPairs){
  logMsg(config, 'Collating read pairs which were removed because the break read aligned to the vector.', logFile)
  collateSampleReads('virusReadsVector')
  collateSampleReads('breakReadsVector')
}


# Organize read files by reference genome and read source (virusRead or breakRead).
# This will allow samples to be aligned to different genomes defined in the samples table.
d <- tibble(file = list.files(file.path(config$outputDir, 'sampleReads'), pattern = '\\.virusReads\\.|\\.breakReads\\.'),
            uniqueSample = unlist(lapply(strsplit(file, '\\.'), function(x) paste0(x[1:(length(x) - 2)], collapse = '.'))),
            source = ifelse(grepl('virusRead', file), 'virusReads', 'breakReads')) %>%
     left_join(select(samples, uniqueSample, refGenomeBLATdb), by = 'uniqueSample')



# Align the sample FASTA files defined in the previous data frame to their respective genomes. 
# This is done by grouping reads by type (virusReads / breakReads) and reference geneome 
# and using parLapply() to disrubte read chunks across a number of CPUs. Unique reads 
# are nested and then unnnested after alignment.

alignments <- bind_rows(lapply(split(d, paste(d$source, d$refGenomeBLATdb)), function(x){
  f <- tmpFile()
  
  # Concatenat all the source / genome fasta files into a single file.
  system(paste('cat ', paste0(file.path(config$outputDir, 'sampleReads', x$file), collapse = ' '), ' > ', 
               paste0(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))))
  
  fasta <- readFasta(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  invisible(file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.fasta'))))
  
  d <- tibble(id  = as.character(fasta@id),
              seq = as.character(fasta@sread)) %>%
       group_by(seq) %>%
        summarise(ids = list(id)) %>%
       ungroup() %>%
       mutate(seqID = paste0('s', 1:n()))
  
  
  fasta <- DNAStringSet(d$seq)
  names(fasta) <- d$seqID
  
  db <- x$refGenomeBLATdb[1]
  cluster <- makeCluster(config$genomAlignment.CPUs)
  clusterExport(cluster, c('config', 'db'), envir = environment())
  
  message(paste0('Starting ', x$source[1], ' alignments'))
  r <- bind_rows(parLapply(cluster, split(fasta, ceiling(seq_along(fasta)/config$alignment.chunk.size)), function(y){
         library(ShortRead)
         library(tidyverse)
         source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
         alignReads.BLAT(y, db)
       }))
  
  stopCluster(cluster)
  
  d <- left_join(select(d, -seq), r,  by = c("seqID" = "qName"))
  d <- unnest(d, ids)
  
  d$source <- x$source[1]
  select(d, -seqID)
}))

save(alignments, file = file.path(config$outputDir, 'alignments.RData'))


# Apply generic alignment filters.
alignments <- 
  filter(alignments, 
         alignmentPercentID >= config$alignment.genome.minPercentID,
         tNumInsert  <= 1, 
         qNumInsert  <= 1,
         tBaseInsert <= 2,
         qBaseInsert <= 2) 


# Apply direction specific filters.
virusReads.aln <- alignments[which(alignments$source == 'virusReads'),]
breakReads.aln <- alignments[which(alignments$source == 'breakReads'),]


# Break point reads have static 5' adapters which are removed with cutadapt.
# Only consider reads which immediately align to the genome.
breakReads.aln <- filter(breakReads.aln, qStart <= 3, matches >= config$trimmedRead.minLength)


# Virus reads have had over-read sequences removed.
# Only consider virus reads where the end of the read aligns to the genome.
virusReads.aln <- filter(virusReads.aln, (qSize - qEnd) <= 3, matches >= config$trimmedRead.minLength)


# If we are not capturing LTR sequences then a static LTR was removed from viral reads
# and we can apply a filter to the beginning of the trimmed read.
if(! config$virusReads.captureLTRseqs){
  logMsg(config, 'Requiring virus reads to begin alignments near their start because the ITR/LTR should of been fully trimmed off.', logFile)
  virusReads.aln <- filter(virusReads.aln, qStart <= 3)
}


alignments <- bind_rows(virusReads.aln, breakReads.aln)


# Rename alignment column headers by appending read sources.
alignments <- lapply(split(alignments, alignments$source), function(x){
    names(x) <- paste0(names(x), '.', x$source[1])
    x
})


# Combine read alignments into rows and extract sample names from read ids. 
frags <- 
  left_join(alignments[["virusReads"]], alignments[["breakReads"]], by = c('ids.virusReads' = 'ids.breakReads')) %>%
  select(ids.virusReads, qStart.virusReads, strand.virusReads, strand.breakReads, tName.virusReads, tName.breakReads,
         tStart.virusReads, tStart.breakReads, tEnd.virusReads, tEnd.breakReads) %>%
  drop_na() %>%
  mutate(uniqueSample = unlist(lapply(str_split(ids.virusReads, '\\|'), '[[', 2 )),
         readID = sub('\\|.+$', '', ids.virusReads)) 


# Remove alignment pairings which can not be real fragments.
# Store indexes in a variable first because indexing on a null returns nothing.
i <- which(frags$tName.virusReads != frags$tName.breakReads)
if(length(i) > 0) frags <- frags[-i,]

i <- which(frags$strand.virusReads == frags$strand.breakReads)
if(length(i) > 0) frags <- frags[-i,]


# Determine the fragment start and end positions and apply additional sanity checks. 
frags <- mutate(frags, 
                fragStart = ifelse(strand.virusReads == '+', tStart.virusReads + 1, tStart.breakReads + 1),
                fragEnd   = ifelse(strand.virusReads == '+', tEnd.breakReads + 1,   tEnd.virusReads + 1),
                fragTest  = ifelse(strand.virusReads == '+', tStart.virusReads < tEnd.breakReads, tStart.breakReads < tEnd.virusReads),  
                fragWidth = (fragEnd - fragStart) + 1) %>%
         filter(fragTest == TRUE, 
                fragWidth <= config$fragments.maxLength,
                fragWidth >= config$fragments.minLength)


# Create a list of reads whoes alignments could not be reduced to a single sane pairing,
# remove conflicting alignments, and create fragment IDs which cab be used to define fragments.

conflictedReads <- unique(frags$ids.virusReads[duplicated(frags$ids.virusReads)])

if(config$filter.removeMultiHitReadPairs){
  logMsg(config, paste0('Removin multi-hit read pairs (', length(conflictedReads), ' read pairs).'), logFile)
  
  if (length(conflictedReads) > 0){
    frags <- filter(frags, ! ids.virusReads %in% conflictedReads) 
  }
}

frags <- mutate(frags, fragID = paste0(uniqueSample, ';', tName.virusReads, ';', strand.virusReads, ';', fragStart, ';', fragEnd))
  
if(config$removeReadFragsWithSameRandomID){
  
  logMsg(config, 'Removing read pairs which share random linker IDs between samples.', logFile)
  
  randomIDs <- bind_rows(lapply(list.files(file.path(config$outputDir, 'tmp'), 
                                           pattern = 'randomIDs', 
                                           full.names = TRUE), function(x){ load(x); randomIDs }))
  frags <- left_join(frags, randomIDs, by = 'readID')


  # Unpack the unique ids.
  frags$subject  <- unlist(lapply(strsplit(frags$uniqueSample, '~'), '[[', 1))
  frags$sample   <- unlist(lapply(strsplit(frags$uniqueSample, '~'), '[[', 2))
  frags$replicate <- unlist(lapply(strsplit(frags$uniqueSample, '~'), '[[', 3))

  frags <- group_by(frags, randomSeqID) %>%
           mutate(samplesPerRandomSeqID = n_distinct(sample)) %>%
           ungroup()

  percentReadsRemoved <- sprintf("%.2f%%", (nrow(subset(frags, samplesPerRandomSeqID > 1)) / nrow(frags))*100) 
  logMsg(config, paste0('Removed ', percentReadsRemoved, ' of fragment reads because one or more other reads shared the same random linker sequence.'), logFile)
  
  frags <- subset(frags, samplesPerRandomSeqID == 1)
}

createFragReadAlignments(config, samples, frags)

save(frags, file = file.path(config$outputDir, 'readFrags.RData'))

cluster <- makeCluster(config$demultiplexing.CPUs)
clusterExport(cluster, c('config'), envir = environment())

# Add caputred LTR sequence information if requested.
if(config$virusReads.captureLTRseqs){
  
  logMsg(config, 'Starting fragment creation in LTR/ITR capture mode.', logFile)
  
  # Add previous mapped LTR sequences to fragments via read ids.
  # Add additional NTs, NTs found between the recognized LTR sequences and the start of genomic alignments. 
  # These addition NTs may be the result of DNA repair mechanisms or alignemnt errors. 
  
  LTR.table <- bind_rows(lapply(list.files(file.path(config$outputDir, 'tmp'), pattern = 'LTRseqs', full.names = TRUE), readRDS))
  
  frags$LTRseq  <- LTR.table[match(frags$readID, LTR.table$id),]$LTRseq
  frags$LTRName <- LTR.table[match(frags$readID, LTR.table$id),]$LTRname
  frags$LTRread <- LTR.table[match(frags$readID, LTR.table$id),]$read
  frags$additionalLTRnts <- substr(frags$LTRread, 1, frags$qStart.virusReads)
  frags$LTRseq2 <- paste0(frags$LTRseq, frags$additionalLTRnts)

  
  # Group reads into fragments and create lists of LTRsequences. 
  frags <- group_by(frags, fragID) %>%
    summarise(reads = n(), readIDs = list(readID), LTRseqs = list(LTRseq), LTRseq2s = list(LTRseq2)) %>%
    ungroup() %>%
    separate(fragID, c('uniqueSample', 'seqnames', 'strand', 'start', 'end'), sep = ';') 
  
  
  # Identify representative LTR sequences for each fragment.
  frags$s <- ntile(1:nrow(frags), config$demultiplexing.CPUs)
  
  frags <- bind_rows(parLapply(cluster, split(frags, frags$s), function(x){
    source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
    library(dplyr)
    x$LTRseqsRep  <- unlist(lapply(x$LTRseqs,  representativeSeq))
    x$LTRseq2sRep <- unlist(lapply(x$LTRseq2s, representativeSeq))
    x
  }))
  
  
  # Standardize the fragments with the gintools package.
  # The standardization has the potential of combining fragments with notably differnt LTR sequences 
  # which may be from separate insertion event. For each fragment, we identify a representive 
  # sequence and remove fragments within each group which is dissimiliar to the representative sequence.
  frags$subject <- unlist(lapply(strsplit(frags$uniqueSample, '~'), '[[', 1))
  frags$sample  <- unlist(lapply(strsplit(frags$uniqueSample, '~'), '[[', 2))
  g <- makeGRangesFromDataFrame(frags, keep.extra.columns = TRUE)
  
  if(config$standardizeBreakPointsBy == 'replicate'){
    g$s <- g$uniqueSample
  } else if (config$standardizeBreakPointsBy == 'sample'){
    g$s <- g$sample
  } else {
    g$s <- g$subject
  }
  
  g <- unlist(GRangesList(parLapply(cluster, split(g, g$s), function(x) gintools::refine_breakpoints(x, counts.col = 'reads')))) 
  
  
  if(config$standardizeSitesBy == 'replicate'){
    g$s <- g$uniqueSample
  } else if (config$standardizeSitesBy == 'sample'){
    g$s <- g$sample
  } else {
    g$s <- g$subject
  }
  
  # There appears to be a ~ 10,000 site limit for standardize sites.
  # Create a 10K splitting vector to use if sites exceed 10K. 
  # This is not an ideal patch because all sites from a subject should be normalized together.
  
  #g <- unlist(GRangesList(parLapply(cluster, split(g, g$s), function(x){
  g <- unlist(GRangesList(lapply(split(g, g$s), function(x){
         if(length(x) <= 10000){
           return(gintools::standardize_sites(x, counts.col = 'reads'))
         } else {
           message('   Applying standardize_sites() patch for large number of sites.')
           i <- unlist(lapply(1:10000, function(i) rep(i, 10000)))
           x$s <- i[1:length(x)]
           return(unlist(GRangesList(lapply(split(x, x$s), function(x2) gintools::standardize_sites(x2, counts.col = 'reads')))))
         }
       }))) 
  
  f <- select(data.frame(g), -LTRseqs, -LTRseq2s, -width, -uniqueSample, -s)

  frags <- bind_rows(parLapply(cluster, split(f, paste(f$subject, f$sample, f$seqnames, f$start, f$end, f$strand)), function(x){
             source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
             library(stringdist)
             library(dplyr)
         
             r <- representativeSeq(x$LTRseq2sRep)
             n <- ceiling(nchar(r)*(config$virusReads.captureLTRseqs.repSeqEditDistMaxPercentage/100))
             
             x <- x[which(stringdist(r, x$LTRseq2sRep) <= n),]
             if(nrow(x) == 0) return(tibble())
             
             x$LTRseqsRep  <- representativeSeq(x$LTRseqsRep)
             x$LTRseq2sRep <- r
             x$reads       <- sum(x$reads)
             x$readIDs     <- list(unique(unlist(x$readIDs)))
             
             x <- x[1,]
             
             x$fragStart <- ifelse(x$strand == '+', x$start, x$end)
             x$posid     <- paste0(x$seqnames, x$strand, x$fragStart)
             x$posid2    <- paste0(x$subject, '_', x$posid)
             x$estAbund  <- 1
             x
           }))
  
  
  createIntUCSCTrack(frags,  title = 'AAVengeR_fragments', 
                     outputFile = file.path(config$outputDir,'AAVengeR_fragments.ucsc'), siteLabel = 'posid2')
  
  
  sites <- bind_rows(parLapply(cluster, split(frags, paste(frags$subject, frags$sample, frags$posid)), function(x){
             source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
             library(stringdist)
             r <- representativeSeq(x$LTRseq2sRep)
             n <- ceiling(nchar(r)*(config$virusReads.captureLTRseqs.repSeqEditDistMaxPercentage/100))
             o <- x[which(stringdist(r, x$LTRseq2sRep) <= n),]
             o$diffRepSeqFragsRemoved = nrow(x) - nrow(o)
             o$LTRseqsRep  <- representativeSeq(o$LTRseqsRep)
             o$LTRseq2sRep <- r
             o$reads       <- sum(o$reads)
             o$estAbund    <- nrow(o)
             
             o$start <- o$fragStart 
             o$end <- o$start
             
             o[1,]
          })) %>% 
         group_by(subject, sample) %>%
         mutate(relAbund = estAbund / sum(estAbund)) %>%
         ungroup()
  
  sites$LTRseqName <- LTR.table[match(sites$LTRseqsRep, LTR.table$LTRseq),]$LTRname
  
  createIntUCSCTrack(sites,  title = 'AAVengeR_sites', 
                     outputFile = file.path(config$outputDir,'AAVengeR_sites.ucsc'), siteLabel = 'posid2')
  
  system(paste('cat', file.path(config$outputDir,'AAVengeR_fragments.ucsc'), file.path(config$outputDir,'AAVengeR_sites.ucsc'), 
               '>', file.path(config$outputDir,'AAVengeR.ucsc')))
  
} else {
  
  logMsg(config, 'Starting fragment creation in fixed LTR/ITR mode.', logFile)
  
  # Group reads into fragments and create lists of LTRsequences. 
  frags <- group_by(frags, fragID) %>%
    summarise(reads = n(), readIDs = list(readID)) %>%
    ungroup() %>%
    separate(fragID, c('uniqueSample', 'seqnames', 'strand', 'start', 'end'), sep = ';') 
  
  
  # Unpack the unique sample ids and standardize fragments.
  frags$subject  <- unlist(lapply(strsplit(frags$uniqueSample, '~'),  '[[', 1))
  frags$sample   <- unlist(lapply(strsplit(frags$uniqueSample, '~'),  '[[', 2))
  frags$replicate <- unlist(lapply(strsplit(frags$uniqueSample, '~'), '[[', 3))
  
  g <- makeGRangesFromDataFrame(frags, keep.extra.columns = TRUE)
  
  if(config$standardizeBreakPointsBy == 'replicate'){
    g$s <- g$uniqueSample
  } else if (config$standardizeBreakPointsBy == 'sample'){
    g$s <- g$sample
  } else {
    g$s <- g$subject
  }
  
  g <- unlist(GRangesList(parLapply(cluster, split(g, g$s), function(x) gintools::refine_breakpoints(x, counts.col = 'reads')))) 
  
  if(config$standardizeSitesBy == 'replicate'){
    g$s <- g$uniqueSample
  } else if (config$standardizeSitesBy == 'sample'){
    g$s <- g$sample
  } else {
    g$s <- g$subject
  }
  
  g <- unlist(GRangesList(parLapply(cluster, split(g, g$s), function(x) gintools::standardize_sites(x, counts.col = 'reads')))) 
  
  # Replicate level merging of standardized fragments.

  repFrags <- data.frame(g)
  repFrags$seqnames <- as.character(repFrags$seqnames)
  repFrags$strand   <- as.character(repFrags$strand)

  # Drop read columns for speed up.
  repFrags <- select(repFrags, -width, -s, -readIDs)
  
  repFrags$s <- paste(repFrags$uniqueSample, repFrags$seqnames, repFrags$start, repFrags$end, repFrags$strand)

  repFrags <- group_by(repFrags, s) %>%
              mutate(reads = sum(reads),
                     posid = paste0(seqnames, strand, ifelse(strand == '+', start, end))) %>%
              dplyr::slice(1) %>%
              ungroup()
      
  
  # Add replicate collapse feature here -- remove replicate value from uniqueSample           
  
  sites <- group_by(repFrags, uniqueSample, posid) %>%
           mutate(fragWidth = end - start + 1,
                  estAbund = n_distinct(fragWidth),
                  reads = sum(reads),
                  start = ifelse(strand == '+', start, end),
                  end = start) %>%
           dplyr::slice(1) %>%
           ungroup() %>%
           group_by(subject, sample) %>%
           mutate(relAbund = estAbund / sum(estAbund)) %>%
           ungroup()
}

stopCluster(cluster)

save(sites, file = file.path(config$outputDir, 'sites.RData'))

logMsg(config, 'done.', logFile)
