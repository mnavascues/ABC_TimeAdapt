set.seed(1234567890)

# write header
write("i Ne",file="results/Ne.txt")


# READ DATA

rec_map_file <- paste0("data/RecombinationMap2010/genetic_map_GRCh37_chr22.txt")
recombination_map <- read.table(rec_map_file,header=T)
L <- recombination_map$Position.bp.[nrow(recombination_map)]

sample_data <- read.table("data/SampleData.txt",header=T)
total_sample_size <- sum(sample_data$size)
missing <- rbeta(total_sample_size,0.4,1) # to be substituted with real proportion of missing data



na=nrow(sample_data)-1
i=1;N=200;n0=10;mu=1.25e-08
generation_time <- runif(1,20,30)
all_different_generations<-F
while (!all_different_generations){
  ts<-NULL
  for (samp in seq_len(na)+1){
    ts <- c(ts,round(runif(1,sample_data$age_y[samp],sample_data$age_o[samp])/generation_time))
  }
  if (ts[na]>ts[na-1]) all_different_generations <- !any(duplicated(ts))
}
tf        <- max(ts)+1
pyslim_ts <- sort(ts) 
slim_ts   <- sort(tf-ts)+1

# SLiM (DEMOGRAPHY) 
seed.slim   <- round(runif(1,0,2^32-1))
system2(command="slim",
        args=c("-d", paste0("N=", N[i]),
               "-d", paste0("tf=", tf),
               "-d", paste0("ts=\"c("), paste(slim_ts,collapse=","), paste0(")\""),
               "-d", paste0("i=", i),
               "-d", paste0("na=", na),
               "-d", paste0("L=", L),
               "-s", seed.slim,
               "src/model.demography.slim > /tmp/slimout.txt"))
  
  # pyslim (RECAPITATION+MUTATION)
seed.pyslim <- round(runif(1,0,2^32-1))
system2(command="python3",
        args=c("src/add.mutations.py",
               i, N[i], mu[i], seed.pyslim, n0, na, ts,
               " > /tmp/pyslimout.txt"))
  
 
  
  
# awk (MISSIG DATA + PSEUDOHAPLOID aDNA)
seed.awk    <- round(runif(1,0,2^32-1))
source("src/awk_postprocessing_template.R")
system2(command="awk",
        args=awk_args)
# sed (UNPHASE)
system2(command="sed",
        args=c("'s/|/\\//g'","results/genotypes2.vcf > results/genotypes3.vcf"))














  
 
# }
sumstats <- read.table("results/sumstats.txt",header=T)
sumstats <- sumstats[order(sumstats$i),]
Ne <- read.table("results/Ne.txt",header=T)
Ne <- Ne[order(Ne$i),]

# Joins parameters, latent variables and statistics in reference table
referenceTable <- as.data.frame(cbind(N,mu,Ne["Ne"],sumstats[c("pi","S")]))
#saveRDS(referenceTable,"results/referenceTable.rds")
#referenceTable <- readRDS("results/referenceTable.rds")

# learning from reference table & get posterior
library(abcrf)
modelN     <- regAbcrf(log10(N)~pi+S, referenceTable, ntree=2000, paral=T, ncores=ncores)
posteriorN <- predict(modelN, HareSumStats, referenceTable, rf.weights=TRUE, paral=T, ncores=ncores)

modelNe     <- regAbcrf(log10(Ne)~pi+S, referenceTable, ntree=2000, paral=T, ncores=ncores)
posteriorNe <- predict(modelNe, HareSumStats, referenceTable, rf.weights=TRUE, paral=T, ncores=ncores)

library(weights)
library(graphics)
library(gplots)
# plot out-of-bag predictions for verification
pdf(file="results/OOBNe.pdf",width=4,height=4)
par(mar=c(4.5,4.5,0.5,0.5)+0.1)
hist2d(data.frame(log10(referenceTable$Ne),modelNe$model.rf$predictions),
       nbins=40,
       same.scale=T,
       col=grey.colors(40,start=1,end=0),
       FUN=function(x) log(length(x)),
       xlim=c(1,3),ylim=c(1,3),
       xlab=expression(log[10]*italic(N)[e]),
       ylab=expression(log[10]*italic(hat(N))[e]))
abline(0,1,col=rgb(213, 94,  0,150,maxColorValue=255),lwd=2)
dev.off()

# plot prior and posterior together
pdf(file="results/posteriorNe.pdf",width=4,height=3)
par(mar=c(4,4,0.2,0.2)+0.1)
hist(x      = log10(referenceTable$Ne),
     breaks = seq(0.9,3.1,0.06),
     col    = rgb(0, 114, 178, 100, maxColorValue=255),
     freq   = FALSE,
     ylim   = c(0,4),
     xlim   = c(1,3),
     main   = "",
     xlab   = expression(log[10]*italic(N)[e]),
     ylab   = "probability density")
wtd.hist(x      = log10(referenceTable$Ne),
         breaks = seq(0.9,3.1,0.06),
         col    = rgb(213, 94, 0, 150, maxColorValue=255),
         freq   = FALSE,
         add    = TRUE,
         weight = posteriorNe$weights)
points(x=posteriorNe$expectation,y=4,pch=16,
       col=rgb(213, 94, 0, 255, maxColorValue=255))
arrows(x0=posteriorNe$expectation, y0=4,
       x1=posteriorNe$quantiles[1],y1=4,
       length=0.1,lwd=2,
       col=rgb(213, 94, 0, 255, maxColorValue=255))
arrows(x0=posteriorNe$expectation, y0=4,
       x1=posteriorNe$quantiles[2],y1=4,
       length=0.1,lwd=2,
       col=rgb(213, 94, 0, 255, maxColorValue=255))
box()
dev.off()

