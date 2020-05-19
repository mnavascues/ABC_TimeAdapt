# Read HapMap recombination map file and save it in a table to be read by SLiM
dir.create("results/r_map_slim")

for (chr in 1:22){
  r_map_file <- paste0("data/RecombinationMap2010/genetic_map_GRCh37_chr",chr,".txt")
  r_map <- read.table(rec_map_file,header=T)
  ends  <- r_map$Position.bp.-1
  rates <- c(0,r_map$Rate.cM.Mb.[1:(nrow(r_map)-1)] * 1e-8)
  write.table( cbind(ends,rates),
               paste0("results/r_map_slim/r_map_chr",chr,"_slim.txt"),
               sep="\t",
               col.names=FALSE,
               row.names=FALSE)
}

