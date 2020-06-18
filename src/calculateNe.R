args <- commandArgs()
N        <- args[1]
N
parentID <- args[2]
parentID
prob <- tabulate(parentID)/N
pIBD <- sum(prob*prob*0.25)
Ne   <- 1/pIBD
