
library(rootSolve)

#setwd("...")

source("functions_2raters.R")

########################### Simulate an example 2x2 rating data ##########################
N <- 200 ## sample size
p1 <- 0.1; p2 <- 0.3 ## marginal probability of positive response
TRAC <- 0.5 ## true TRAC estimand value

f12 <- calc_f12_from_Fa(Fa=TRAC) ## transform to, true f12 (log cross product ratio) estimand value
## calc_Fa_from_f12(f12=f12) #doublecheck

p11 <- calc_p11(p1=p1, p2=p2, f12=f12) ## true underlying probabilities of the 2x2 table
p10 <- p1-p11
p01 <- p2-p11
p00 <- 1-p11-p10-p01
pa <- p11 + p00
kappa <- (pa-p1*p2-(1-p1)*(1-p2))/(1-p1*p2-(1-p1)*(1-p2)) ## true kappa estimand
phi <- (p11-p1*p2)/sqrt(p1*(1-p1)*p2*(1-p2)) ## true Phi estimand

## Simulate a single 2x2 table
set.seed(1234)
y <- as.numeric(rmultinom(1, N, c(p11, p10, p01, p00)))
y

##########################################################################################


########################### Estimation of f12 and TRAC ###################################

## prespecify a bias-correction pseudo count uniformly added to each 2x2 data cell
CC <- 0.5

## calculate the point and variance estimate of f12 and F(a) using the SA method
f12_est <- calc_f12est(y=y, CC=CC)
f12_var_SA <- calc_f12var_SA(y=y, CC=CC)
c(f12_est, f12_var_SA)

Fa_est <- calc_Fa_from_f12(f12=f12_est)
a <- 0.5 ## a-specification in TRAC
Fa_var_SA <- a^2/4*(1-Fa_est^2)^2*f12_var_SA
c(Fa_est, Fa_var_SA)

## 95% CI construction for f12 using the SA method 
c(f12_est - qnorm(0.975)*sqrt(f12_var_SA),
  f12_est + qnorm(0.975)*sqrt(f12_var_SA))
## 95% CI construction for F(a) using the SA method 
c(Fa_est - qnorm(0.975)*sqrt(Fa_var_SA),
  Fa_est + qnorm(0.975)*sqrt(Fa_var_SA))


## 95% CI construction for f12 using the LRT method 
test_LRTI <- try(calc_CI_LRTI(y=y+CC), silent=T)
if(class(test_LRTI)!="try-error"){
  test_LRTI
}
## 95% CI construction for F(a) using the LRT method
test_LRTI_Fa <- try(calc_CI_LRTI_Fa(y=y+CC), silent=T)
if(class(test_LRTI_Fa)!="try-error"){
  test_LRTI_Fa
}


## 95% CI construction for f12 using the GOF method 
test_GOF_profile <- try(calc_CI_GOF_profile(y=y+CC), silent=T)
if(class(test_GOF_profile)!="try-error"){
  test_GOF_profile
}
## 95% CI construction for F(a) using the GOF method
test_GOF_profile_Fa <- try(calc_CI_GOF_profile_Fa(y=y+CC), silent=T)
if(class(test_GOF_profile_Fa)!="try-error"){
  test_GOF_profile_Fa
}

########################################################################################


######################## Estimation of several existing IRA methods ####################

## phi
phi_CC <- calc_phi(y=y+CC)
phi_CC[1]
phi_CC[1]-qnorm(0.975)*phi_CC[2]
phi_CC[1]+qnorm(0.975)*phi_CC[2]

## kappa
kappa_CC <- calc_kappa(y=y+CC)
kappa_CC[1]
kappa_CC[1]-qnorm(0.975)*kappa_CC[2]
kappa_CC[1]+qnorm(0.975)*kappa_CC[2]

## S
S_CC <- calc_S(y=y+CC)
S_CC[1]
S_CC[1]-qnorm(0.975)*S_CC[2]
S_CC[1]+qnorm(0.975)*S_CC[2]

## AC1
AC1_CC <- calc_AC1(y=y+CC)
AC1_CC[1]
AC1_CC[1]-qnorm(0.975)*AC1_CC[2]
AC1_CC[1]+qnorm(0.975)*AC1_CC[2]

#######################################################################################