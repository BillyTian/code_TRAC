#setwd("...)

source("functions_Rraters.R")

############################ Example code for case R=4 ############################
R <- 4 ## Number of raters
C <- 2 ## Binary rating

## Set true probabilities of 2^R cells
p_eq=0.3
p <- array(0, dim=rep(C, R))
p[1,1,1,1] <- p_eq*2*0.01
p[2,2,2,2] <- p_eq*2*0.99
p[1,1,1,2] <- p[1,1,2,1] <- p[1,2,1,1] <- p[2,1,1,1] <- (1-2*p_eq)/22*2
p[2,2,1,1] <- p[2,1,2,1] <- p[2,1,1,2] <- p[1,2,2,1] <- p[1,2,1,2] <- p[1,1,2,2] <- (1-2*p_eq)/22
p[2,2,2,1] <- p[2,2,1,2] <- p[2,1,2,2] <- p[1,2,2,2] <- (1-2*p_eq)/22*2

## Unify the order of specifying probability parameters
param_grid <- expand.grid(
  rater1 = 1:C,
  rater2 = 1:C,
  rater3 = 1:C,
  rater4 = 1:C
)

## Arrange the p_vector using p set above, based on the ORDER in param_grid
p_vector <- NULL
for (i in 1:nrow(param_grid)){
  p_vector[i] <- p[param_grid[i,1], param_grid[i,2], param_grid[i,3], param_grid[i,4]] 
}


## Calculate the true f12_R and F(a)_R given the true probability pattern
f12_R_true <- calc_f12_R(p_vector = p_vector)[1]
f12_R_true
Fa_R_true <- calc_Fa_from_f12(f12=f12_R_true)
Fa_R_true

## Calculate the true estimand values of the multi-rater kappa (Conger's kappa) and Gwet's AC1 for comparison
conger_true <- calc_conger_fromp_R4(p)
conger_true
gwet_true <- calc_gwet_fromp_R4(p)
gwet_true

## Simulate a 2^R-multinomial dataset 
set.seed(1234)
y <- as.numeric(rmultinom(1, N, p_vector))

CC <- 0.5/4 ## recommended bias-correction term (mimicking 0.5-continuity correction in marginal 2x2 tables)
p_CC4 <- (y+0.5/4)/sum(y+0.5/4) ## p_hat, needed in estimation procedures for F(a)_R

## Calculate f12_R estimate and its corresponding N*var
calculate_4 <- calc_f12_R(p_vector = p_CC4)
est_f12_R <- calculate_4[1]
se_f12_R <- sqrt(1/N*calculate_4[2])
a <- 0.5
est_Fa_R <- calc_Fa_from_f12(f12=est_f12_R, a=a)
se_Fa_R <- sqrt((exp(est_f12_R)/(exp(a*est_f12_R)+1)^4))*se_f12_R

## Point estimate and 95% CI for f12_R and F(a)_R
c(est_f12_R, est_f12_R-qnorm(0.975)*se_f12_R, est_f12_R+qnorm(0.975)*se_f12_R)
c(est_Fa_R, est_Fa_R-qnorm(0.975)*se_Fa_R, est_Fa_R+qnorm(0.975)*se_Fa_R)


### For comparison, using existing software to calculate kappa and AC1 estimates
library(irrCAC)
dat <- param_grid[rep(1:nrow(param_grid), times=y),]
res_conger <- irrCAC::conger.kappa.raw(dat)
res_gwet <- irrCAC::gwet.ac1.raw(dat)

est_conger <- res_conger$est$coeff.val
se_conger <- res_conger$est$coeff.se

est_gwet <- res_gwet$est$coeff.val
se_gwet <- res_gwet$est$coeff.se

c(est_conger, 
  est_conger - qnorm(0.975)*se_conger, 
  est_conger + qnorm(0.975)*se_conger)

c(est_gwet, 
  est_gwet - qnorm(0.975)*se_gwet, 
  est_gwet + qnorm(0.975)*se_gwet)

