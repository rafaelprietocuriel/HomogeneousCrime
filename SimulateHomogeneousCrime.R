
#### code created to simulate uniform crime
#### N individuals suffer a crime rate following a Poisson distribution
#### the rate of the uniform distribution is homogeneous, c
#### with 1,000 simulations of a crime rate in (0, 1) it gives observed
#### concentration metrics for homogeneous crime rates

Rate      <- c()
Top1      <- c()
Top5      <- c()
Top10     <- c()
Top20     <- c()
EntropyR  <- c()
GiniR     <- c()
GiniCorrected <- c()
AverageCrimes <- c()
RECC      <- c()
PoisGamma <- c()

sims <- 1000
for (j in 1:sims){
  rate <- runif(1)
  crimes <- rpois(n = 1000000, lambda = rate)
  sortedCrimes <- sort(crimes, decreasing = TRUE)
  top.1 <- sum(sortedCrimes[1:10000]) / sum(crimes)
  top.5 <- sum(sortedCrimes[1:50000]) / sum(crimes)
  top.10 <- sum(sortedCrimes[1:100000]) / sum(crimes)
  top.20 <- sum(sortedCrimes[1:200000]) / sum(crimes)
  Rate <- c(Rate, rate)
  Top1 <- c(Top1, top.1)
  Top5 <- c(Top5, top.5)
  Top10 <- c(Top10, top.10)
  Top20 <- c(Top20, top.20)
  GiniR <- c(GiniR, Gini(crimes))
  EntropyR <- c(EntropyR, entropy(crimes))
  AverageCrimes <- c(AverageCrimes, sum(crimes)/sum(crimes>0))
  GiniCorrected <- c(GiniCorrected, Gini(sortedCrimes[1:sum(crimes)]))
  CamanModel <- mixalg(crimes, family = "poisson", startk = 3, acc = 10^(-7))
  CamanRates <- rep(CamanModel@t, times = CamanModel@p * 1000000)
  gg <- Gini(CamanRates)
  RECC <- c(RECC, gg)
  Ntotal <- length(crimes)
  pars <- fitdistr( crimes , "negative binomial")$estimate
  simulated_gam <- rgamma( Ntotal , shape=pars[ 1 ] , rate = 1)
  sorted_gam <- sort( simulated_gam , decreasing=T)
  normalized_gam <- sorted_gam/sum( sorted_gam )
  gini =(1 / Ntotal ) * (2 * sum(cumsum(normalized_gam)) - Ntotal - 1)
  PoisGamma <- c(PoisGamma, gini)
}

Results <- data.frame(Rate = Rate,
                      Top1 = Top1,
                      Top5 = Top5,
                      Top10 = Top10,
                      Top20 = Top20,
                      GiniR = GiniR, 
                      EntropyR = EntropyR, 
                      AverageCrimes = AverageCrimes, 
                      GiniCorrected = GiniCorrected,
                      RECC = ~RECC,
                      PoisGamma = PoisGamma)
save(Results, file = "SimulateHomogeneousCrime.RData")
