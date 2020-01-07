library(ShortRead)
library(dplyr)
library(readr)
library(parallel)
library(yaml)

config <- read_yaml('config.yml')
config$megahit.path <- '/home/everett/software/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit'

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

R1.ids <- sub('\\s+.+$', '', as.character(R1@id))
R2.ids <- sub('\\s+.+$', '', as.character(R2@id))

R1 <- sread(R1)
R2 <- sread(R2)

names(R1) <- R1.ids
names(R2) <- R2.ids

R1.unique <- unique(R1)
R2.unique <- unique(R2)

cluster <- makeCluster(config$CPUs)
clusterExport(cluster, 'config')

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

# Blast all unique readsa against the vector sequences.
parallelBlat(R1.unique, config$CPUs, 'vectorDB.blat.path', label = 'unique_R1')
parallelBlat(R2.unique, config$CPUs, 'vectorDB.blat.path', label = 'unique_R2')

system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'unique_R1.psl', full.names = TRUE), collapse = ' '), ' > unique_R1.vector.psl'))
system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'unique_R2.psl', full.names = TRUE), collapse = ' '), ' > unique_R2.vector.psl'))

R1.unique.vector.blat <- parseBLAToutput('unique_R1.vector.psl')
R2.unique.vector.blat <- parseBLAToutput('unique_R2.vector.psl')

# Reduce the reads to only reads that partially align to the vectors.
R1.unique <- R1.unique[names(R1.unique) %in% R1.unique.vector.blat$qName]
R2.unique <- R2.unique[names(R2.unique) %in% R2.unique.vector.blat$qName]

writeFasta(R1.unique, file = 'unique_R1.vectorAligned.fasta')
writeFasta(R2.unique, file = 'unique_R2.vectorAligned.fasta')

if(dir.exists('MEGAHIT_singleEnd_output')) unlink('MEGAHIT_singleEnd_output1', recursive = TRUE)
system(paste0(config$megahit.path, ' -r unique_R1.vectorAligned.fasta,unique_R2.vectorAligned.fasta -o MEGAHIT_singleEnd_output1 ',
              '--min-count 1 --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))
c1 <- readFasta('MEGAHIT_singleEnd_output1/final.contigs.fa')
c1@id <- BStringSet(paste(as.character(c1@id), ' min-count: 1'))

if(dir.exists('MEGAHIT_singleEnd_output')) unlink('MEGAHIT_singleEnd_output2', recursive = TRUE)
system(paste0(config$megahit.path, ' -r unique_R1.vectorAligned.fasta,unique_R2.vectorAligned.fasta -o MEGAHIT_singleEnd_output2 ',
              '--min-count 2 --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))
c2 <- readFasta('MEGAHIT_singleEnd_output2/final.contigs.fa')
c2@id <- BStringSet(paste(as.character(c2@id), ' min-count: 2'))


if(dir.exists('MEGAHIT_singleEnd_output')) unlink('MEGAHIT_singleEnd_output3', recursive = TRUE)
system(paste0(config$megahit.path, ' -r unique_R1.vectorAligned.fasta,unique_R2.vectorAligned.fasta -o MEGAHIT_singleEnd_output3 ',
              '--min-count 3 --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))
c3 <- readFasta('MEGAHIT_singleEnd_output3/final.contigs.fa')
c3@id <- BStringSet(paste(as.character(c3@id), ' min-count: 3'))

writeFasta(Reduce('append', list(c1, c2, c3)), file = 'MEGAHIT_singleEnd.contigs.fasta')



invisible(file.remove(list.files('tmp', full.names = TRUE))) 

# Blast all unique readsa against the vector sequences.
parallelBlat(R1.unique, config$CPUs, 'genomeDB.blat.path', label = 'unique_R1')
parallelBlat(R2.unique, config$CPUs, 'genomeDB.blat.path', label = 'unique_R2')

system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'unique_R1.psl', full.names = TRUE), collapse = ' '), ' > unique_R1.genome.psl'))
system(paste0('cat ', paste0(list.files(path = 'tmp', pattern = 'unique_R2.psl', full.names = TRUE), collapse = ' '), ' > unique_R2.genome.psl'))

stopCluster(cluster)

R1.blat <- parseBLAToutput('unique_R1.genome.psl')
R2.blat <- parseBLAToutput('unique_R2.genome.psl')

R1.blat.flankingDownStream <- subset(R1.blat, tName == config$siteChromosome & tStart >= (config$sitePosition - 500) & tStart <= (config$sitePosition - 25))
R2.blat.flankingDownStream <- subset(R2.blat, tName == config$siteChromosome & tStart >= (config$sitePosition - 500) & tStart <= (config$sitePosition - 25))
R2.blat.flankingUpStream   <- subset(R2.blat, tName == config$siteChromosome & tEnd >= (config$sitePosition + 25) & tEnd <= (config$sitePosition + 500))
R1.blat.flankingUpStream   <- subset(R1.blat, tName == config$siteChromosome & tEnd >= (config$sitePosition + 25) & tEnd <= (config$sitePosition + 500))

unique.downStream.bridgingReads <- c(R1.unique[names(R1.unique) %in% R1.blat.flankingDownStream$qName],
                                     R2.unique[names(R2.unique) %in% R2.blat.flankingDownStream$qName])

unique.upStream.bridgingReads <- c(R1.unique[names(R1.unique) %in% R1.blat.flankingUpStream$qName],
                                   R2.unique[names(R2.unique) %in% R2.blat.flankingUpStream$qName])

invisible(file.remove(list.files('tmp', full.names = TRUE))) 

writeFasta(unique.downStream.bridgingReads, file = 'unique.downStream.bridgingReads.fasta')
writeFasta(unique.upStream.bridgingReads,   file = 'unique.upStream.bridgingReads.fasta')



if(dir.exists('MEGAHIT_singleEnd_downstream_output')) unlink('MEGAHIT_singleEnd_downstream_output', recursive = TRUE)
system(paste0(config$megahit.path, ' -r unique.downStream.bridgingReads.fasta -o MEGAHIT_singleEnd_downstream_output ',
              '--min-count 1 --k-list ', paste0(seq(from = 15, to = 127, by = 6), collapse = ',')))



if(dir.exists('MEGAHIT_singleEnd_upstream_output')) unlink('MEGAHIT_singleEnd_upstream_output', recursive = TRUE)
system(paste0(config$megahit.path, ' -r unique.upStream.bridgingReads.fasta -o MEGAHIT_singleEnd_upstream_output ',
              '--min-count 1 --k-list ', paste0(seq(from = 15, to = 127, by = 6), collapse = ',')))


system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' ', config$R1.path, ' ', config$R2.path, ' > R1R2.vector.sam'))
system('samtools view -S -b R1R2.vector.sam > R1R2.vector.bam')     
system('samtools sort -f R1R2.vector.bam R1R2.vector.sorted.bam')                          
system('samtools index R1R2.vector.sorted.bam') 
file.remove(c('R1R2.vector.sam', 'R1R2.vector.bam'))


if(length(unique.downStream.bridgingReads) > 0){
  system(paste0(config$bwa.path, ' mem -M ', config$genomeDB.bwa.path, ' unique.downStream.bridgingReads.fasta > unique.downStream.bridgingReads.sam'))
  system('samtools view -S -b unique.downStream.bridgingReads.sam > unique.downStream.bridgingReads.bam')     
  system('samtools sort -f unique.downStream.bridgingReads.bam unique.downStream.bridgingReads.sorted.bam')                          
  system('samtools index unique.downStream.bridgingReads.sorted.bam')                    
  file.remove(c('unique.downStream.bridgingReads.sam', 'unique.downStream.bridgingReads.bam'))

  system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' unique.downStream.bridgingReads.fasta > unique.downStream.bridgingReads.vector.sam'))
  system('samtools view -S -b unique.downStream.bridgingReads.vector.sam > unique.downStream.bridgingReads.vector.bam')     
  system('samtools sort -f unique.downStream.bridgingReads.vector.bam unique.downStream.bridgingReads.vector.sorted.bam')                          
  system('samtools index unique.downStream.bridgingReads.vector.sorted.bam') 
  file.remove(c('unique.downStream.bridgingReads.vector.sam', 'unique.downStream.bridgingReads.vector.bam'))
}


if(length(unique.upStream.bridgingReads) > 0){
  system(paste0(config$bwa.path, ' mem -M ', config$genomeDB.bwa.path, ' unique.upStream.bridgingReads.fasta > unique.upStream.bridgingReads.sam'))
  system('samtools view -S -b unique.upStream.bridgingReads.sam > unique.upStream.bridgingReads.bam')     
  system('samtools sort -f unique.upStream.bridgingReads.bam unique.upStream.bridgingReads.sorted.bam')                          
  system('samtools index unique.upStream.bridgingReads.sorted.bam')                    
  file.remove(c('unique.upStream.bridgingReads.sam', 'unique.upStream.bridgingReads.bam'))
  
  system(paste0(config$bwa.path, ' mem -M ', config$vectorDB.bwa.path, ' unique.upStream.bridgingReads.fasta > unique.upStream.bridgingReads.vector.sam'))
  system('samtools view -S -b unique.upStream.bridgingReads.vector.sam > unique.upStream.bridgingReads.vector.bam')     
  system('samtools sort -f unique.upStream.bridgingReads.vector.bam unique.upStream.bridgingReads.vector.sorted.bam')                          
  system('samtools index unique.upStream.bridgingReads.vector.sorted.bam') 
  file.remove(c('unique.upStream.bridgingReads.vector.sam', 'unique.upStream.bridgingReads.vector.bam'))
}

unlink('tmp', recursive = TRUE)
dir.create('alignments')
system('mv *.bam alignments/')
system('mv *.bai alignments/')
invisible(file.remove(list.files(pattern = '*.psl')))

write(c('new',
      'genome canFam3',
      paste0('load  ', getwd(), '/alignments/unique.downStream.bridgingReads.sorted.bam'),
      paste0('load  ', getwd(), '/alignments/unique.upStream.bridgingReads.sorted.bam'),
      paste0('snapshotDirectory ', getwd()),
      paste0('goto ', config$siteChromosome, ':', config$sitePosition-1000, '-',  config$sitePosition+1000),         
      'sort position',
      'collapse',
      'snapshot genome.bridgingReads.png',
      'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/unique.downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/unique.upStream.bridgingReads.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto HeavyChainVector:1-4042',         
        'sort position',
        'collapse',
        'snapshot vector.heavy.bridgingReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/unique.downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/unique.upStream.bridgingReads.vector.sorted.bam'),
        paste0('snapshotDirectory ', getwd()),
        'goto LightChainVector:1-3904',         
        'sort position',
        'collapse',
        'snapshot vector.light.bridgingReads.png',
        'exit'), file = 'batchFile', append = FALSE)

system('xvfb-run --auto-servernum  /home/everett/software/IGV_Linux_2.7.2/igv.sh -b batchFile')


write(c('new',
        'genome /home/everett/releases/Canine_hemophilia_AAV/plasmids/IGV_genomes/VectorSeqs_ITR_to_ITR.genome',
        paste0('load  ', getwd(), '/alignments/unique.downStream.bridgingReads.vector.sorted.bam'),
        paste0('load  ', getwd(), '/alignments/unique.upStream.bridgingReads.vector.sorted.bam'),
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







