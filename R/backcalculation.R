# back-calculation methods

# TODO refactor

rm(list=ls())
# libraries
library(dplyr)
library(ggplot2)
library(here)

# define back-calculation formulas

# define biological intercept model (Campana, pulled from Vigliola & Meekan 2009 )
BI <- function(pars, dat){
  # pars should be L0 and R0
  L <- c()
  for(i in 1:length(dat$Lcpt)){
    L[i] <- dat$Lcpt[i] + ((dat$R_i[i] - dat$Rcpt[i]) * ((dat$Lcpt[i] - pars[1])/(dat$Rcpt[i] - pars[2])))
  }
  return(L)
}

# define Fraser-Lee back-calculation model (pulled from Vigliola & Meekan 2009)
FL <- function(pars, dat){ # pars here is just b0
  L_pred <- c()
  for(i in 1:length(dat$idx1)){
    L_pred[i] <- pars + (dat$Lcpt[i]-pars)*(dat$R_i[i]/dat$Rcpt[i])
  }
  return(L_pred)
}


load(here("data-raw", "incdat_2021.Rdata"))

# numbers for indexing
nrec <- length(idat$id)

# parameters - use sum of squares to figure out what Mike used
df <- idat |> 
  select(Lcpt = TL_mm, Rcpt = M_radius, R_i = Radius, L)

# define sum of squares function using biological intercept model
# TODO fix this
ss <- function(pars, dat){
  L_pred <- c()
  for (i in 1:nrec){
    L_pred[i] <- dat$Lcpt[i] + ((dat$R_i[i] - dat$Rcpt[i]) * ((dat$Lcpt[i] - pars[1])/(dat$Rcpt[i] - pars[2])))
  }
  return(sum((L_pred-dat$L)^2))
}

# minimize sum of squared differences
pars=c(L0=0, R0=0) # starting values
fit <- nlminb(pars, ss, dat=df)
fit$par
fit$objective

# 03-22-2021 results
# L0 = 21.6734241
# R0 = 0.1370712
L0p <- fit$par[1]
R0p <- fit$par[2]



# now compare BI model to Mike's back-calculated data
L4 <- c()
L4 <- BI(c(L0p, R0p), df) # back-calculate lengths using R0 and L0 from Mike's data

plot(idat$L, L4)
plot(idat$Age, L4)
resid <- idat$L-L4
hist(resid) # this looks pretty good
summary(resid)

# now add all 3 BCM to the incdat file
rm(list=ls())
load('incdat_2021.Rdata')

# name the columns correctly
df <- data.frame(idat$TL_mm, idat$M_radius, idat$Radius, idat$L, idat$idx1)
colnames(df) <- c('Lcpt', 'Rcpt', 'R_i', 'L', 'idx1')
pars = 30 #b0, from MNDNR handbook
idat$L_FL <- FL(pars, df)
plot(idat$Age, idat$L_FL)
save(idat, file='incdat_2021-BCM1.Rdata')