logMsg <- function(config, msg, logFile, append = TRUE){
  library(lubridate)
  # Create log time elapsed stamp using the start time in the config object, eg.  [15.5 minutes] log message.
  te <- paste0('[', 
               sprintf("%.1f", 
                       as.numeric(difftime(ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S")), 
                                           config$startTime, 
                                           units = 'min'))), 
               ' minutes]')
  write(paste0(te, '\t', msg), file = logFile, append = append)
}


checkConfigFilePaths <- function(files){
  invisible(sapply(unique(samples$refGenomeBLATdb), function(file){ 
    if(! file.exists(file)){
      logMsg(config, paste0('Error. "', file, '" does not exist.'), logFile)
      stop('Stopping AVVengeR -- file not found error.')
    }
  }))
}


syncReads <-function(...){
  arguments <- list(...)
  
  # Create a list of read IDs common to all read arguments.
  n <- Reduce(base::intersect, lapply(arguments, names))

  lapply(arguments, function(x){ 
    x <- x[names(x) %in% n]; 
    x[order(names(x))]
  })
}



waitForFile <- function(f, seconds = 1){
  repeat
  {
    if(file.exists(f)) break
    Sys.sleep(seconds)
  }
  return(TRUE)
}




createSeqChunks <- function(f, chunkSize, label, ouputDir){
  
  # Create a pointer like object to the read file.
  strm <- FastqStreamer(f, n = as.integer(chunkSize))
  
  # Extract chunks from file and write out chunks with numeric suffixes, eg. virusReads.5
  n <- 1
  repeat {
    fq <- yield(strm)
    if(length(fq) == 0) break
    writeFastq(fq, file = file.path(ouputDir, paste0(label, '.', n)), compress = FALSE)
    n <- n + 1
  }
}


tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }


shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
}




representativeSeq <- function(s){
  if(length(s) == 1 | n_distinct(s) == 1) return(s[1])
  m <- as.matrix(stringdist::stringdistmatrix(s))
  d <- apply(m, 1, sum) 
  s[which(d == min(d))[1]]
}



captureLTRseqs <- function(reads, seqs){
  
  d <- group_by(tibble(ids = names(reads), 
                       read = as.character(reads)), read) %>%
    summarise(id.list = paste0(ids, collapse=',')) %>%
    ungroup() %>%
    mutate(id = paste0('s', 1:n()),
           qlength = nchar(read))
  
  f <- tmpFile()
  write(paste0('>', d$id, '\n', d$read), file = file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  
  seqsNames <- unlist(lapply(unlist(strsplit(seqs, ';')), function(x){
    x <- unlist(strsplit(x, ','))
    x[1]
  }))
  
  invisible(lapply(unlist(strsplit(seqs, ';')), function(x){
    x <- unlist(strsplit(x, ','))
    write(paste0('>', x[1], '\n', x[2]), file =  file.path(config$outputDir, 'tmp', paste0(f, '.db')), append = TRUE)
  }))
  
  system(paste0(config$command.makeblastdb, ' -in ', file.path(config$outputDir, 'tmp', paste0(f, '.db')), ' -dbtype nucl -out ', file.path(config$outputDir, 'tmp', f)), ignore.stderr = TRUE)
  
  waitForFile(file.path(file.path(config$outputDir, 'tmp', paste0(f, '.nin'))))
  
  system(paste0(config$command.blastn, ' -word_size 7 -evalue 50 -outfmt 6 -query ',
                file.path(config$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                file.path(config$outputDir, 'tmp', f),
                ' -num_threads 5 -out ', file.path(config$outputDir, 'tmp', paste0(f, '.blast'))), 
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  waitForFile(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), seconds = 1)
  
  
  if(file.info(file.path(config$outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(character(length = 0))
  
  b <- read.table(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
  names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
  b$alignmentLength <- b$qend - b$qstart + 1
  b <- left_join(b, d, by = c('qname' = 'id'))
  
  b <- subset(b, pident >= config$virusReads.captureLTRseqs.minPercentSeqID & 
                 qstart <= config$virusReads.captureLTRseqs.maxAlignmentStart & 
                 gapopen <= config$virusReads.captureLTRseqs.maxGapOpen)
  
  # For each read alignment, find the longest alignment and use the ordering of the LTR sequences as a tie breaker.
 b <- bind_rows(lapply(split(b, b$read), function(x){
    x <- subset(x, alignmentLength ==  max(x$alignmentLength))
    x$n <- match(x$sseqid, seqsNames)
    arrange(x, n) %>% slice(1)
  }))
 
 b$id.list <- strsplit(b$id.list, ',')
 b <- unnest(b, id.list)  
 
 b$LTRseq <- substr(b$read, b$qstart, b$qend)
 
 reads <- reads[names(reads) %in% b$id.list]
 
 b <- select(b, id.list, sseqid, LTRseq)
 names(b) <- c('id', 'LTRname', 'LTRseq')
 
 return(list(reads = reads, LTRs = b))
}

golayCorrection <- function(x){
  fq <- readFastq(x)
  fq.ids <- sub('\\s+.+$', '', as.character(fq@id))
  
  writeFasta(fq, file = file.path(paste0(x, '.fasta')))
  
  system(paste(config$command.python2, file.path(config$softwareDir, 'bin', 'golayCorrection', 'processGolay.py'), file.path(paste0(x, '.fasta'))))
  invisible(file.remove(c(paste0(x, '.fasta'), x)))
  
  corrected <- readFasta(file.path(paste0(x, '.fasta.corrected')))
  
  a <- left_join(tibble(id = fq.ids, seq = as.character(sread(fq))),
                 tibble(id = as.character(corrected@id), seq2 = as.character(sread(corrected))), by = 'id')
  a <- rowwise(a) %>% mutate(editDist = as.integer(adist(seq, seq2))) %>% ungroup()
  
  i <- which(a$editDist <= 2)
  a[i,]$seq <- a[i,]$seq2
  
  if(! all(fq.ids == a$id)) stop('There was an ordering error during the Golay correction step.')
  
  o <- paste0('@', as.character(fq@id), '\n', a$seq, '\n+\n', as.character(quality(fq@quality)))
  fileConn <- file(x)
  writeLines(o, fileConn)
  close(fileConn)
}





collateSampleReads <- function(label){
  v <- tibble(file = list.files(file.path(config$outputDir, 'tmp'), pattern = paste0('\\.', label, '\\.')),
              sample = unlist(lapply(strsplit(file, paste0('\\.', label, '\\.')), '[[', 1)))
  invisible(lapply(split(v, v$sample), function(x){
    system(paste('cat ', paste0(file.path(config$outputDir, 'tmp', x$file), collapse = ' '), ' > ', 
                 paste0(file.path(config$outputDir, 'sampleReads', paste0(x$sample[1], '.', label, '.fasta')))))
  }))
}



alignReads.BLAT <- function(x, db){
  f <- tmpFile()
  writeFasta(x, file = file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  
  comm <- paste0(config$command.blat, ' ', db, ' ', 
                 file.path(config$outputDir, 'tmp', paste0(f, '.fasta')), ' ', 
                 file.path(config$outputDir, 'tmp', paste0(f, '.psl')),  
                 ' -tileSize=11 -stepSize=9 -minIdentity=90 -out=psl -noHead')
  system(comm)
  file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  
  if(file.exists(file.path(config$outputDir, 'tmp', paste0(f, '.psl')))){
    b <- parseBLAToutput(file.path(config$outputDir, 'tmp', paste0(f, '.psl')))
    file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.psl')))
    b
  } else {
    return(tibble())
  }
}

# BLAT output parser.
parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble())
  b <- read_delim(f, delim = '\t', col_names = FALSE, col_types = cols())
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$queryPercentID       <- (b$matches/b$qSize)*100
  b$tAlignmentWidth      <- (b$tEnd - b$tStart) + 1
  b$queryWidth           <- (b$qEnd - b$qStart) + 1
  b$alignmentPercentID   <- (b$matches/b$tAlignmentWidth)*100
  b$percentQueryCoverage <- (b$queryWidth/b$qSize)*100
  b$qStarts              <- as.character(b$qStarts)
  b$tStarts              <- as.character(b$tStarts)
  b
}






trimLeadingSeq <- function(x, seq){
  f <- tmpFile()
  
  writeFasta(x,  file = file.path(config$outputDir, 'tmp', f))
  
  system(paste0(config$command.cutadapt3, ' -f fasta  -e 0.15 -g ', seq, ' --overlap 2 ', 
                file.path(config$outputDir, 'tmp', f), ' > ', 
                file.path(config$outputDir, 'tmp', paste0(f, '.out'))),
         ignore.stderr = TRUE)
  
  s <- readFasta(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
  x <- s@sread
  names(x) <- as.character(s@id)
  file.remove(file.path(config$outputDir, 'tmp', f))
  file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
  x
}


trimOverReadSeq <- function(x, seq){
  f <- tmpFile()
  
  writeFasta(x,  file = file.path(config$outputDir, 'tmp', f))
  
  system(paste0(config$command.cutadapt3, ' -f fasta  -e 0.15 -a ', seq, ' --overlap 2 ', 
                file.path(config$outputDir, 'tmp', f), ' > ', 
                file.path(config$outputDir, 'tmp', paste0(f, '.out'))),
         ignore.stderr = TRUE)
  
  s <- readFasta(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
  x <- s@sread
  names(x) <- as.character(s@id)
  file.remove(file.path(config$outputDir, 'tmp', f))
  file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.out')))
  x
}


getVectorReadIDs <- function(reads, config, vectorFile){
  d <- group_by(tibble(ids = names(reads), 
                       read = as.character(reads)), read) %>%
    summarise(id.list = paste0(ids, collapse=',')) %>%
    ungroup() %>%
    mutate(id = paste0('s', 1:n()),
           qlength = nchar(read))
  
    f <- tmpFile()
    write(paste0('>', d$id, '\n', d$read), file = file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))

    system(paste0(config$command.makeblastdb, ' -in ', vectorFile, ' -dbtype nucl -out ', file.path(config$outputDir, 'tmp', f)), ignore.stderr = TRUE)
   
    waitForFile(file.path(file.path(config$outputDir, 'tmp', paste0(f, '.nin'))))
    
    system(paste0(config$command.blastn, ' -word_size 10 -evalue 10 -outfmt 6 -query ',
                  file.path(config$outputDir, 'tmp', paste0(f, '.fasta')), ' -db ',
                  file.path(config$outputDir, 'tmp', f),
                  ' -num_threads 5 -out ', file.path(config$outputDir, 'tmp', paste0(f, '.blast'))), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
  
    waitForFile(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), seconds = 1)
  
    if(file.info(file.path(config$outputDir, 'tmp', paste0(f, '.blast')))$size == 0) return(character(length = 0))
    
    b <- read.table(file.path(config$outputDir, 'tmp', paste0(f, '.blast')), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    b <- left_join(b, d, by = c('qname' = 'id'))
    b$pcoverage <- ((b$qend - b$qstart) / b$qlength)*100
  
    b <- subset(b, pident >= config$alignment.vector.minPercentID & pcoverage >= config$alignment.vector.minPercentQueryCoverage & gapopen <= 1)

    invisible(file.remove(list.files(file.path(config$outputDir, 'tmp'), pattern = f, full.names = TRUE)))
    unlist(strsplit(b$id.list, ','))
}


createFragReadAlignments <- function(config, samples, frags){
  invisible(lapply(split(samples, samples$subject), function(x){
  
    f <- list.files(file.path(config$outputDir, 'sampleReads'), pattern = paste0(unique(x$sample), collapse = '|'), full.names = TRUE)
  
    virusReads <- Reduce('append', lapply(f[grep('virus', f)], readFasta))
    breakReads <- Reduce('append', lapply(f[grep('break', f)], readFasta))
  
    breakReads <- breakReads[as.character(breakReads@id) %in% frags$ids.virusReads]
    virusReads <- virusReads[as.character(virusReads@id) %in% frags$ids.virusReads]
  
    breakReads@id <- BStringSet(sub('\\|.+$', '', as.character(breakReads@id)))
    virusReads@id <- BStringSet(sub('\\|.+$', '', as.character(virusReads@id)))
  
    if(! all(as.character(breakReads@id) == as.character(virusReads@id))) stop('Read id mismatch error.')
  
    writeFasta(breakReads, file = file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.breakReads.fasta')))
    writeFasta(virusReads, file = file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.virusReads.fasta')))
  
    system(paste0(config$command.bwa, ' mem -M ', x$refGenomeBWAdb[1], ' ',
                  file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.breakReads.fasta')), ' ',
                  file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.virusReads.fasta')), ' > ',
                  file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sam'))))
  
    system(paste0(config$command.samtools, ' view -S -b ',  
                  file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sam')), ' > ',
                  file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.bam'))))
  
    invisible(file.remove(file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sam'))))
  
    system(paste0(config$command.samtools, ' sort -o ', 
                  file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sorted.bam')), ' ',
                  file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.bam'))))
  
    system(paste0(config$command.samtools, ' index ', 
                  file.path(config$outputDir, 'fragReads', paste0(x$subject[1], '.sorted.bam'))))
  }))
}



createIntUCSCTrack <- function(d, abundCuts = c(5,10,50), 
                               posColors = c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB"),
                               negColors = c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000"),
                               title = 'intSites', outputFile = 'track.ucsc', visibility = 1, 
                               padSite = 0,
                               siteLabel = NA){
  
  # Check function inputs.
  if(length(posColors) != length(negColors)) 
    stop('The pos and neg color vectors are not the same length.')
  
  if(length(abundCuts) != length(posColors) - 1) 
    stop('The number of aundance cut offs must be one less than the number of provided colors.')
  
  if(! all(c('start', 'end', 'strand', 'seqnames', 'estAbund') %in% names(d))) 
    stop("The expected column names 'start', 'end', 'strand', 'seqnames', 'estAbund' were not found.") 
  
  if(is.na(siteLabel) | ! siteLabel %in% names(d)) 
    stop('The siteLabel parameter is not defined or can not be found in your data.')
  
  
  # Cut the abundance data. Abundance bins will be used to look up color codes.
  # We flank the provided cut break points with 0 and Inf in order to bin all values outside of breaks.
  cuts <- cut(d$estAbund, breaks = c(0, abundCuts, Inf), labels = FALSE)
  
  posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
  negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')
  
  
  # Create data fields needed for track table.
  d$score <- 0
  d$color <- ifelse(d$strand == '+', posColors[cuts], negColors[cuts])
  
  # Pad the site n NTs to increase visibility.
  if(padSite > 0){
    d$start <- floor(d$start - padSite/2)
    d$end   <- ceiling(d$end + padSite/2)
  }
  
  # Define track header.
  trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s",
                       title, title, visibility)
  
  # Write out track table.
  write(trackHead, file = outputFile, append = FALSE)
  write.table(d[, c('seqnames', 'start', 'end', siteLabel, 'score', 'strand', 'start', 'end', 'color')], 
              sep = '\t', col.names = FALSE, row.names = FALSE, file = outputFile, append = TRUE, quote = FALSE)
}
