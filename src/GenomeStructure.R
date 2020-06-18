# get centromeres position from
# Genome Decoration Page, NCBI.
# Ideogram data for Homo sapience (850 bphs, Assembly GRCh38.p12) 2018
# https://www.ncbi.nlm.nih.gov/genome/tools/gdp
genome_anotation <- read.table("data/ideogram_9606_GCF_000001305.15_850_V1",
                               header=T,
                               fill=T,
                               comment.char = "%",
                               row.names = NULL)
cen_start <- intersect(which(genome_anotation$X.chromosome==1),which(genome_anotation$arm=="p"))
cen_start <- intersect(cen_start,which(genome_anotation$stain=="acen"))
cen_start <- genome_anotation$bp_start[cen_start]
cen_end   <- intersect(which(genome_anotation$X.chromosome==1),which(genome_anotation$arm=="q"))
cen_end   <- intersect(cen_end,which(genome_anotation$stain=="acen"))
cen_end   <- genome_anotation$bp_stop[cen_end]
centromeres <- c(cen_start,cen_end) 
for (chr in 2:22){
  cen_start <- intersect(which(genome_anotation$X.chromosome==chr),which(genome_anotation$arm=="p"))
  cen_start <- intersect(cen_start,which(genome_anotation$stain=="acen"))
  cen_start <- genome_anotation$bp_start[cen_start]
  cen_end   <- intersect(which(genome_anotation$X.chromosome==chr),which(genome_anotation$arm=="q"))
  cen_end   <- intersect(cen_end,which(genome_anotation$stain=="acen"))
  cen_end   <- genome_anotation$bp_stop[cen_end]
  centromeres <- rbind(centromeres, c(cen_start,cen_end) )
}
colnames(centromeres) <- c("start","end")
rownames(centromeres) <- paste0("chr",1:22)

(centromeres)


# read data for hu7man genome from stdpopsim
HumanGenome <- read.table("data/stdpopsimHumanGenome.txt",header=T)
(HumanGenome)

# read recombination map from HapMap files
positions        <- 1
rates <- numeric()
cumulativeLength <- 0.0
chr<-1
chromosomes <- 1
for (chr in 1:22){
  centromeres[chr,] <- centromeres[chr,] + cumulativeLength
  
  rec_map_file <- paste0("data/RecombinationMap2010/genetic_map_GRCh37_chr",chr,".txt")
  r_map <- read.table(rec_map_file,header=T)
  positions <- c(positions, cumulativeLength + r_map$Position.bp.)
  cumulativeLength <- cumulativeLength + r_map$Position.bp.[nrow(r_map)]
  chromosomes <- c(chromosomes, cumulativeLength)
  if (chr!=22) {
    positions   <- c(positions, cumulativeLength + 1)
    chromosomes <- c(chromosomes, cumulativeLength + 1)
  }

  rates     <- c(rates,  HumanGenome$Rec_rate[chr], r_map$Rate.cM.Mb.[1:(nrow(r_map)-1)] * 1e-8)
  if (chr!=22) rates <- c(rates, 0.5) else  rates <- c(rates, 0.0)
  
  #ends  <- r_map$Position.bp.-1
  #rates <- c(0,r_map$Rate.cM.Mb.[1:(nrow(r_map)-1)] * 1e-8)
  #write.table( cbind(ends,rates),
  #             paste0("results/r_map_slim/r_map_chr",chr,"_slim.txt"),
  #             sep="\t",
  #             col.names=FALSE,
  #             row.names=FALSE)
}

length(positions);head(positions);tail(positions)
length(rates);head(rates);tail(rates)

(chromosomes)
(centromeres)


