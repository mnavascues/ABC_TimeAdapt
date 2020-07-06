# get centromeres position from
# Genome Decoration Page, NCBI.
# Ideogram data for Homo sapience (850 bphs, Assembly GRCh38.p12) 2018
# https://www.ncbi.nlm.nih.gov/genome/tools/gdp
genome_anotation <- read.table("data/ideogram_9606_GCF_000001305.15_850_V1",
                               header=T,
                               fill=T,
                               comment.char = "%",
                               row.names = NULL)

centromeres <- matrix(NA, nrow=22,ncol=2,dimnames=list(paste0("chr",1:22),c("start","end")))
for (chr in 1:22){
  cen_start <- intersect(which(genome_anotation$X.chromosome==chr),which(genome_anotation$arm=="p"))
  cen_start <- intersect(cen_start,which(genome_anotation$stain=="acen"))
  cen_start <- genome_anotation$bp_start[cen_start]
  cen_end   <- intersect(which(genome_anotation$X.chromosome==chr),which(genome_anotation$arm=="q"))
  cen_end   <- intersect(cen_end,which(genome_anotation$stain=="acen"))
  cen_end   <- genome_anotation$bp_stop[cen_end]
  centromeres[chr,] <-  c(cen_start,cen_end)
}
(centromeres)


# read data for hu7man genome from stdpopsim
HumanGenome <- read.table("data/stdpopsimHumanGenome.txt",header=T)
(HumanGenome)

# read recombination map from HapMap files
positions        <- 1
rates <- numeric()
cumulativeLength <- 0.0
chr<-1

block_size <- 100

chromosomes <- matrix(NA, nrow=22,ncol=2,dimnames=list(paste0("chr",1:22),c("start","end")))
chromosomes[1,1] <- 1
for (chr in 1:22){
  centromeres[chr,] <- centromeres[chr,] + cumulativeLength
  
  rec_map_file <- paste0("data/RecombinationMap2010/genetic_map_GRCh37_chr",chr,".txt")
  r_map <- read.table(rec_map_file,header=T)
  
  positions <- c(positions, cumulativeLength + r_map$Position.bp.)
  cumulativeLength <- cumulativeLength + r_map$Position.bp.[nrow(r_map)]
  chromosomes[chr,2] <- cumulativeLength
  if (chr!=22) {
    positions   <- c(positions, cumulativeLength + 1)
    chromosomes[chr+1,1] <- cumulativeLength + 1
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

# set positions to start at 0 for SLiM and msprime

positions   <- positions-1
chromosomes <- chromosomes-1
centromeres <- centromeres-1

positions_A <- sort(c(chromosomes))
rates_A <- numeric()
for (chr in 1:22){
  rates_A <- c(rates_A,HumanGenome$Rec_rate[chr],0.5)
}
recombination_map_A <- data.frame(positions=positions_A,
                                 rates=rates_A)

dim(recombination_map_A)
write.table(recombination_map_A,file="data/recombination_map_msprime_A.txt",col.names=F,row.names=F)
recombination_map_A <- cbind(positions_A[-1],rates_A[-length(rates_A)])
write.table(recombination_map_A,file="data/recombination_map_slim_A.txt",col.names=F,row.names=F)


write.table(chromosomes,file="data/chromosomes.txt",col.names=F,row.names=F)
write.table(centromeres,file="data/centromeres.txt",col.names=F,row.names=F)



reference   <- round(log10(rates[1]))
positions_B <- positions[1]
rates_B     <- 10^round(log10(rates[1]))
for (i in 2:length(positions)){
  print(i)
  if (round(log10(rates[i]))!=reference){
    reference <- rates[i]
    positions_B <- c(positions_B,positions[i])
    if (rates[i]==0.5){
      rates_B <- c(rates_B,0.5)
    }else{
      rates_B <- c(rates_B,10^round(log10(rates[i])))
    }
  }
}

  recombination_map_B <- data.frame(positions=positions_B,
                                  rates=rates_B)

dim(recombination_map_B)

dim(recombination_map_B)
write.table(recombination_map_B,file="data/recombination_map_msprime_B.txt",col.names=F,row.names=F)
recombination_map_B <- cbind(positions_B[-1],rates_B[-length(rates_B)])
write.table(recombination_map_B,file="data/recombination_map_slim_B.txt",col.names=F,row.names=F)

