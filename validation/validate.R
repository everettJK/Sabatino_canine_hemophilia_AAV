library(ShortRead)
library(dplyr)
library(readr)
library(parallel)
library(yaml)
library(ggplot2)
options(stringsAsFactors = FALSE)

invisible(lapply(rev(list.dirs()), function(workDir){
if(workDir == '.') return()
setwd(workDir)
message('\n\nStarting ', workDir, '\n\n')
if(file.exists('done')) return()
#--------------------------------------------------------------------------------------------------

config <- read_yaml('config.yml')
config$megahit.path  <- '/home/everett/software/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit'
config$blastBin.path <- '/home/everett/ext/blast+/bin'


shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
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

R1 <- readFastq(config$R1.path)
R2 <- readFastq(config$R2.path)

R1 <- trimTailw(R1, 2, '5', 5)
R2 <- trimTailw(R2, 2, '5', 5)

R1 <- shortRead2DNAstringSet(R1)
R2 <- shortRead2DNAstringSet(R2)

r <- syncReads(R1, R2)
R1 <- r[[1]]; R2 <- r[[2]]

cluster <- makeCluster(config$CPUs)
clusterExport(cluster, 'config', envir = environment())

if(! dir.exists('tmp')) dir.create('tmp')
invisible(file.remove(list.files('tmp', full.names = TRUE))) # Clear out previous runs.

parallelBlat <- function(x, CPUs, dbConfigKey, label = 'x'){
  invisible(parLapply(cluster, split(x, ntile(1:length(x), CPUs)), function(x){
    library(ShortRead)
    f <- paste0(tempfile(pattern = 'tmp', tmpdir = './tmp'), '.', label)
    writeFasta(x, file = paste0(f, '.fasta'))
    system(paste0(config$blat.path, ' ', config[[dbConfigKey]], ' ', paste0(f, '.fasta'), ' ', paste0(f, '.psl'), ' -tileSize=11 -stepSize=9 -minIdentity=90 -out=psl -noHead'))
    file.remove(paste0(f, '.fasta'))
   }))
}

# Blast all unique reads against the vector sequences.
parallelBlat(R1, config$CPUs, 'vectorDB.blat.path', label = 'R1')
parallelBlat(R2, config$CPUs, 'vectorDB.blat.path', label = 'R2')

system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'R1.psl', full.names = TRUE), collapse = ' '), ' > R1.vector.psl'))
system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'R2.psl', full.names = TRUE), collapse = ' '), ' > R2.vector.psl'))

R1.vector.blat <- parseBLAToutput('R1.vector.psl')
R2.vector.blat <- parseBLAToutput('R2.vector.psl')

# Reduce the reads to only reads that partially align to the vectors.
R1 <- R1[names(R1) %in% R1.vector.blat$qName]
R2 <- R2[names(R2) %in% R2.vector.blat$qName]
writeFasta(R1, file = 'R1.vectorAligned.fasta')
writeFasta(R2, file = 'R2.vectorAligned.fasta')

r <- syncReads(R1, R2)
writeFasta(r[[1]], file = 'R1.vectorAligned.paired.fasta')
writeFasta(r[[2]], file = 'R2.vectorAligned.paired.fasta')

# Build contigs with different min-count values using both paired-end and single read approaches.
contigs <- lapply(c(2, 3, 4), function(n){
  dir <- paste0('MEGAHIT_pairedEnd_output', n)
  
  system(paste0(config$megahit.path, ' -1 R1.vectorAligned.paired.fasta -2 R2.vectorAligned.paired.fasta -o ', dir, 
                ' --min-count ', n, ' --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))

  c1 <- readFasta(file.path(dir, 'final.contigs.fa'))
  c1@id <- BStringSet(paste0(as.character(c1@id), '_pairedEnd_minCount_', n))
  unlink(dir, recursive = TRUE)

  dir <- paste0('MEGAHIT_singleEnd_output', n)
  system(paste0(config$megahit.path, ' -r R1.vectorAligned.fasta,R2.vectorAligned.fasta -o ', dir, 
                ' --min-count ', n, ' --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))

  c2 <- readFasta(file.path(dir, 'final.contigs.fa'))
  c2@id <- BStringSet(paste0(as.character(c2@id), '_singleEnd_minCount_', n))
  unlink(dir, recursive = TRUE)

  list('pairedEnd' = c1, 'singleEnd' = c2)
})


se <- lapply(contigs, '[[', 'singleEnd')
se <- se[which(unlist(lapply(se, length)) > 0)]

pe <- lapply(contigs, '[[', 'pairedEnd')
pe <- pe[which(unlist(lapply(pe, length)) > 0)]

se <- Reduce('append', se)
pe <- Reduce('append', pe)
r  <- Reduce('append', (list(se, pe)))
r  <- r[width(r) >= 50]
ids <- gsub('\\s', '_', as.character(r@id))
r   <- sread(r)
names(r) <- ids
r  <- unique(r)
r  <- r[order(width(r), decreasing = TRUE)]
writeFasta(r, file = 'contigs.fasta')

system(paste0(config$blastBin.path, '/blastn -db ../../data/sequenceData/vectorSeqs_ITR_to_ITR -query contigs.fasta -word_size 7 -evalue 100 -out contigs.blast -outfmt 6'))

tbl <- read_delim('contigs.blast', '\t', col_names = FALSE)
names(tbl) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

d <- unlist(strsplit(getwd(), '/'))
d <- unlist(strsplit(d[length(d)], '_'))
if(d[1] %in% c('M50', 'M06', 'M66')){
  tbl <- subset(tbl, sseqid == 'SingleChainVector')
  vectorFeatures <- read.csv('../../data/singleChainVectorFeatures.csv', header = TRUE)
  featureColors <- c('red1', 'gray25', 'gold2', 'gray50', 'green2', 'green4', 'blue1', 'gray90', 'red4')
} else {
  tbl <- subset(tbl, sseqid != 'SingleChainVector')
  vectorFeatures <- read.csv('../../data/splitChainVectorFeatures.csv', header = TRUE)
  featureColors <- c('red1', 'gray25', 'gold2', 'gray50', 'green2', 'gray70', 'blue1', 'gray90', 'red4', 'green4', 'gray10', 'gray20')
}

tbl$contigLength = as.integer(unlist(lapply(stringr::str_match_all(tbl$qname, '_len=(\\d+)'), '[', 2)))
tbl <- subset(tbl, pident >= 95 & gapopen == 0)
tbl$strand <- ifelse(tbl$sstart > tbl$send, '-', '+')
tbl$contigStart <- ifelse(tbl$qstart > tbl$qend, tbl$qend, tbl$qstart)
tbl$contigEnd   <- ifelse(tbl$qstart > tbl$qend, tbl$qstart, tbl$qend)
tbl$vectorStart <- ifelse(tbl$sstart > tbl$send, tbl$send, tbl$sstart)
tbl$vectorEnd   <- ifelse(tbl$sstart > tbl$send, tbl$sstart, tbl$send)
tbl$width <- tbl$contigEnd - tbl$contigStart
tbl <- select(tbl, qname, contigStart, contigEnd, vectorStart, vectorEnd, strand, sseqid, width, contigLength)
tbl$n <- 1:nrow(tbl)
tbl$contigID <- paste('Contig', group_indices(tbl, qname))

vectorFeatures$start <- vectorFeatures$start - vectorFeatures$correction
vectorFeatures$end   <- vectorFeatures$end - vectorFeatures$correction
vectorFeatures$n <- 1:nrow(vectorFeatures)

features <- bind_rows(lapply(split(tbl, paste0(tbl$contigID)), function(x){
            y <- 0
            contig <- bind_rows(lapply(split(x, x$n), function(x2){
              if((x2$vectorEnd - x2$vectorStart + 1) >= 25){
                v <- subset(vectorFeatures, sseqid == x2$sseqid)
                bind_rows(lapply(split(v, v$n), function(x3){
                  # Vector coords
                  featureStartTest <- x3$start >= x2$vectorStart & x3$start <= x2$vectorEnd
                  featureEndTest   <- x3$end >= x2$vectorStart & x3$end <= x2$vectorEnd
                  if(featureStartTest == FALSE & featureEndTest == FALSE) return(tribble()) # Skip if neither end of the feature in within the alignment.
                  
                  vectorToContigShift <- x2$vectorStart - x2$contigStart
                  
                  #if(x2$contigID == 'Contig 1' & x3$label == 'Promoter') browser()
                  
                  # Truncate features if needed.
                  if(featureStartTest == FALSE) x3$start <- x2$vectorStart
                  if(featureEndTest == FALSE)   x3$end <- x2$vectorEnd
                  
                  y <<- y + 1
                  ### if(x2$contigStart > x2$vectorStart) vectorToContigShift <- vectorToContigShift * -1
                  return(tibble(contigID = x2$contigID, 
                                start = x3$start - vectorToContigShift, 
                                end =  x3$end - vectorToContigShift,
                                strand = x2$strand, 
                                feature = x3$label,
                                contigLength = x$contigLength[1],
                                y = y))
               
                }))
              } else {
                return(tibble())
              }
           }))
            
          contig
        }))

if(nrow(features) > 0){
  features$feature <- factor(as.character(features$feature), levels = unique(vectorFeatures$label))

  p <- ggplot(features, aes(color = feature)) +
    theme_bw() +
    scale_color_manual(name = 'Vector feature', values = featureColors, drop = FALSE) +
    scale_linetype_manual(name = 'Strand', values = c('dashed', 'solid')) +
    geom_segment(aes(x = start, y = y, xend = end, yend = y, linetype = strand)) +
    geom_vline(aes(xintercept = contigLength)) +
    facet_grid(contigID~., scales="free_y") +
    labs(x = 'Contig position') +
    theme( panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
         axis.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0),
         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         panel.border = element_rect(colour = "black")) 

  ggsave(p, file = 'contig_features.pdf')
}


# Blast all unique reads against the vector sequences.
invisible(file.remove(list.files('tmp', full.names = TRUE))) 
parallelBlat(R1, config$CPUs, 'genomeDB.blat.path', label = 'R1')
parallelBlat(R2, config$CPUs, 'genomeDB.blat.path', label = 'R2')

system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'R1.psl', full.names = TRUE), collapse = ' '), ' > R1.genome.psl'))
system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'R2.psl', full.names = TRUE), collapse = ' '), ' > R2.genome.psl'))

stopCluster(cluster)

R1.blat <- parseBLAToutput('R1.genome.psl')
R2.blat <- parseBLAToutput('R2.genome.psl')

R1.blat.flankingDownStream <- subset(R1.blat, tName == config$siteChromosome & tStart >= (config$sitePosition - 500) & tStart <= (config$sitePosition - 25))
R2.blat.flankingDownStream <- subset(R2.blat, tName == config$siteChromosome & tStart >= (config$sitePosition - 500) & tStart <= (config$sitePosition - 25))
R2.blat.flankingUpStream   <- subset(R2.blat, tName == config$siteChromosome & tEnd >= (config$sitePosition + 25) & tEnd <= (config$sitePosition + 500))
R1.blat.flankingUpStream   <- subset(R1.blat, tName == config$siteChromosome & tEnd >= (config$sitePosition + 25) & tEnd <= (config$sitePosition + 500))

#-----
# blatLocationNeighborhood <- function(x, posDigits = 4, n = 25){
#   x <- x[! (x$tName == 'chrX' & x$tStart >=  122897024 & x$tStart <= 123043178),]
#   x$neighborhood <- paste(x$tName, substr(as.character(x$tStart), 1, posDigits))
#   sort(table(x$neighborhood), decreasing = TRUE)[1:n]
# }
# 
# blatLocationNeighborhood(R1.blat)
# blatLocationNeighborhood(R2.blat)
#-----

downStream.bridgingReads <- c(R1[names(R1) %in% R1.blat.flankingDownStream$qName],
                              R2[names(R2) %in% R2.blat.flankingDownStream$qName])

upStream.bridgingReads <- c(R1[names(R1) %in% R1.blat.flankingUpStream$qName],
                               R2[names(R2) %in% R2.blat.flankingUpStream$qName])

invisible(file.remove(list.files('tmp', full.names = TRUE))) 

writeFasta(downStream.bridgingReads, file = 'downStream.bridgingReads.fasta')
writeFasta(upStream.bridgingReads,   file = 'upStream.bridgingReads.fasta')


system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' ', config$R1.path, ' ', config$R2.path, ' > R1R2.vector.sam'))
system('samtools view -S -b R1R2.vector.sam > R1R2.vector.bam')     
system('samtools sort -f R1R2.vector.bam R1R2.vector.sorted.bam')                          
system('samtools index R1R2.vector.sorted.bam') 
file.remove(c('R1R2.vector.sam', 'R1R2.vector.bam'))


if(length(downStream.bridgingReads) > 0){
  system(paste0(config$bwa.path, ' mem -M ', config$genomeDB.bwa.path, ' downStream.bridgingReads.fasta > downStream.bridgingReads.sam'))
  system('samtools view -S -b downStream.bridgingReads.sam > downStream.bridgingReads.bam')     
  system('samtools sort -f downStream.bridgingReads.bam downStream.bridgingReads.sorted.bam')                          
  system('samtools index downStream.bridgingReads.sorted.bam')                    
  file.remove(c('downStream.bridgingReads.sam', 'downStream.bridgingReads.bam'))

  system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' downStream.bridgingReads.fasta > downStream.bridgingReads.vector.sam'))
  system('samtools view -S -b downStream.bridgingReads.vector.sam > downStream.bridgingReads.vector.bam')     
  system('samtools sort -f downStream.bridgingReads.vector.bam downStream.bridgingReads.vector.sorted.bam')                          
  system('samtools index downStream.bridgingReads.vector.sorted.bam') 
  file.remove(c('downStream.bridgingReads.vector.sam', 'downStream.bridgingReads.vector.bam'))
}


if(length(upStream.bridgingReads) > 0){
  system(paste0(config$bwa.path, ' mem -M ', config$genomeDB.bwa.path, ' upStream.bridgingReads.fasta > upStream.bridgingReads.sam'))
  system('samtools view -S -b upStream.bridgingReads.sam > upStream.bridgingReads.bam')     
  system('samtools sort -f upStream.bridgingReads.bam upStream.bridgingReads.sorted.bam')                          
  system('samtools index upStream.bridgingReads.sorted.bam')                    
  file.remove(c('upStream.bridgingReads.sam', 'upStream.bridgingReads.bam'))
  
  system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' upStream.bridgingReads.fasta > upStream.bridgingReads.vector.sam'))
  system('samtools view -S -b upStream.bridgingReads.vector.sam > upStream.bridgingReads.vector.bam')     
  system('samtools sort -f upStream.bridgingReads.vector.bam upStream.bridgingReads.vector.sorted.bam')                          
  system('samtools index upStream.bridgingReads.vector.sorted.bam') 
  file.remove(c('upStream.bridgingReads.vector.sam', 'upStream.bridgingReads.vector.bam'))
}

unlink('tmp', recursive = TRUE)
dir.create('alignments')
system('mv *.bam alignments/')
system('mv *.bai alignments/')
invisible(file.remove(list.files(pattern = '*.psl')))

write(c('new',
      'genome canFam3',
      paste0('load  ', getwd(), '/alignments/downStream.bridgingReads.sorted.bam'),
      paste0('load  ', getwd(), '/alignments/upStream.bridgingReads.sorted.bam'),
      paste0('snapshotDirectory ', getwd()),
      paste0('goto ', config$siteChromosome, ':', config$sitePosition-1000, '-',  config$sitePosition+1000),         
      'sort position',
      'collapse',
      'snapshot genome.bridgingReads.png',
      'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/upStream.bridgingReads.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto HeavyChainVector:1-4042',         
        'sort position',
        'collapse',
        'snapshot vector.heavy.bridgingReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/upStream.bridgingReads.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto LightChainVector:1-3904',         
        'sort position',
        'collapse',
        'snapshot vector.light.bridgingReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/upStream.bridgingReads.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto SingleChainVector:1-5411',         
        'sort position',
        'collapse',
        'snapshot vector.single.bridgingReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')




write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/R1R2.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto HeavyChainVector:1-4042',         
        'sort position',
        'collapse',
        'snapshot vector.heavy.allReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/R1R2.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto LightChainVector:1-3904',         
        'sort position',
        'collapse',
        'snapshot vector.light.allReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')



write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/R1R2.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto SingleChainVector:1-5411',         
        'sort position',
        'collapse',
        'snapshot vector.single.allReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


invisible(file.remove(list.files('tmp', full.names = TRUE))) 
dir.create('workFiles')
system('mv batchFile alignments *.blast *.fasta workFiles')

write('done', file = 'done')
#--------------------------------------------------------------------------------------------------
setwd('..')
}))





