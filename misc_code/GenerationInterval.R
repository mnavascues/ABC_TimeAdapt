# Calculates generation interval
# following Fenner (2005) doi:10.1002/ajpa.20188
# data from Howell (1979) & Binford (2001) as reported in Fenner (2005)

Mf <- 0.6/100 # percentage of the female hunter-gatherer reproductive population that dies annually (Dobe Ju/’hoansi)
Mm <- Mf      # percentage of the male hunter-gatherer reproductive population that dies annually (Dobe Ju/’hoansi)
Ff  <- 21.4  # Mean maternal age at first birth (Dobe Ju/’hoansi)
Fl  <- 34.4  # Mean maternal age at last birth (Dobe Ju/’hoansi)
Dhg <- mean(c(8.0,6.0,2.5)) # Hunter-gatherer male/female age differential at first marriage (!kung)

If <- numeric() 
for (i in seq(Ff,Fl,1)){
  If <- c(If, i*(1 - Mf*(i-Ff)))
}
If <- sum(If)/(Fl-Ff+1)

Im <- numeric()
for (i in seq(Ff+Dhg,Fl+Dhg,1)){
  Im <- c(Im, i*(1 - Mm*(i-(Ff+Dhg))))
}
Im <- sum(Im)/(Fl-Ff+1)

Ih <- mean(c(Im,If))

print(paste("If =",If,"; Im =",Im,"; Ih =",Ih))

rscaledbeta <- function(n,shape1,shape2,lower_limit,upper_limit){
  return((upper_limit-lower_limit)*rbeta(n,shape1,shape2) + lower_limit)
}

# Generation interval estimated from aDNA and diverge to Neanderthal
# from Moorjani et al (2016) is estimated to be:
lower_limit <- 26
upper_limit <- 30

Moorjani_estimate <- 28.1  # from Moorjani et al (2016)

target_mode <- mean(c(Moorjani_estimate,Ih))
target_mode <- (target_mode-lower_limit)/(upper_limit-lower_limit)

shape1 <- 2 # value somehow subjective; gives sd values larger to error estimate in Moorjani et al (2016)
shape2 <- (shape1-1)/target_mode - shape1 + 2

(upper_limit-lower_limit)*shape1/(shape1+shape2) + lower_limit           # expected_mean 
(upper_limit-lower_limit)*(shape1-1/3)/(shape1+shape2-2/3) + lower_limit # expected_median
(upper_limit-lower_limit)*(shape1-1)/(shape1+shape2-2) + lower_limit     # expected mode

test <- rscaledbeta(100000,shape1,shape2,lower_limit,upper_limit)
hist(test)
abline(v=Ih,col="red")
abline(v=Moorjani_estimate,col="blue")
sd(test)
var(test)
