# simplyfy recombination map

positions <- c(14431347, 14432618, 14433624, 14433659, 14433758,
              14434713, 14434960, 14435070, 14435171, 14435207,
              14439734, 14441016, 14441342, 14449374, 14452292,
              14456927, 14459493, 14479437, 14480059, 14482325)
rates <- c(7.750400, 7.779539, 7.783596, 7.787343, 7.782681,
           7.697552, 7.697816, 0.468435, 0.470604, 0.474626,
           3.967907, 3.964676, 3.964676, 3.889957, 3.888082,
           3.896374, 3.904113, 0.006432, 0.008905, 0.010738)
factor <- 1e-8

round_log10_rates <- round(log10(rates*factor))
rate_change <- which(diff(round(log10(rates*factor)))!=0)+1

new_positions <- positions[rate_change]
new_rates <- 10**round_log10_rates[rate_change]

simplify_rec_map <- function(positions,rates,factor=1){
  round_log10_rates <- round(log10(rates*factor))
  rate_change <- c(1,which(diff(round_log10_rates)!=0)+1)
  new_positions <- positions[rate_change]
  new_rates <- 10**round_log10_rates[rate_change]
  return(data.frame(positions=new_positions,rates=new_rates))
}

hapmap = read.table("data/hapmap/genetic_map_chr22_b36.txt", header=T)

sim_hapmap = simplify_rec_map(hapmap$position,hapmap$COMBINED_rate.cM.Mb.,1e-8)

head(sim_hapmap)
head(hapmap)

tail(sim_hapmap)
tail(hapmap)

hapmap = read.table("data/hapmap/genetic_map_chr21_b36.txt", header=T)
tail(hapmap)
hapmap = read.table("data/hapmap/genetic_map_chr22_b36.txt", header=T)
head(hapmap)


# read recombination map from image
require(digitize)
cal = ReadAndCal("data/bee_recombination/chromosome10.png")
rec=DigitData(col="green")
data = Calibrate(rec, cal, 0.1,5,20,200)# x1, x2, y1, y2
saveRDS(data,file="data/bee_recombination/chromosome10.rds")


data = readRDS(file = paste0("data/bee_recombination/chromosome1.rds"))
genome <- cbind(Chromosome=rep(1,length(data$x)),
                Position=round(data$x*1000000),
                Recombination_rate=data$x*0.01/1000000)
for (i in 2:16){
  data = readRDS(file = paste0("data/bee_recombination/chromosome",i,".rds"))
  print(data$x[length(data$x)]*1000000)
  genome = rbind(genome, cbind(Chromosome=rep(i,length(data$x)),
                               Position=round(data$x*1000000),
                               Recombination_rate=data$x*0.01/1000000))
}

write.table(genome,file="tests/bee_genome.txt",row.names = F)
