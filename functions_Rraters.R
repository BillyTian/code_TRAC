#######################################
#### Functions for the R>2 methods ####
#######################################

## Function for calculating TRAC F(a) from f12
calc_Fa_from_f12 <- function(f12, a=0.5){
  (exp(a*f12)-1)/(exp(a*f12)+1)
}

## Function for calculating f12_R and N*var(f12_R)
calc_f12_R <- function(p_vector, R=4){
  p <- p_vector
  C <- 2
  X <- expand.grid(replicate(R, 1:C, simplify = FALSE)) #size=X x R, order consistent with param_grid
  pairs <- combn(R,2)
  numPairs <- ncol(pairs)
  
  grad <- rep(0, C^R)
  est_cumu <- 0
  for (rp in 1:choose(R,2)) {
    #rp=1
    #pairs[[rp]]
    i <- pairs[1,rp]
    j <- pairs[2,rp]
    
    for(k in 1:C) {
      #k=1
      cond11 <- (X[, i] == k & X[, j] == k)
      cond10 <- (X[, i] == k & X[, j] != k)
      cond01 <- (X[, i] != k & X[, j] == k)
      cond00 <- (X[, i] != k & X[, j] != k)
      
      p11 <- sum(p[cond11])
      p10 <- sum(p[cond10])
      p01 <- sum(p[cond01])
      p00 <- sum(p[cond00])
      
      contr <- rep(0, C^R) ## condxx are positions (indicator functions) where to add the 4 things
      contr[cond11] <- contr[cond11] + 1 / p11
      contr[cond00] <- contr[cond00] + 1 / p00
      contr[cond10] <- contr[cond10] - 1 / p10
      contr[cond01] <- contr[cond01] - 1 / p01
      
      grad <- grad + contr
      
      est_cumu <- est_cumu + log(p11)+log(p00)-log(p10)-log(p01)
      #print(k)
    }
    #print(rp)
  }
  
  grad <- grad/(numPairs * C)
  F_var <- sum(p * grad^2)
  
  F_est <- est_cumu/(numPairs * C)
  
  return(c(F_est, F_var))
}

calc_conger_fromp_R4 <- function(p){
  p11xx <- sum(p[1,1,,])
  p00xx <- sum(p[2,2,,])
  p1x1x <- sum(p[1,,1,])
  p0x0x <- sum(p[2,,2,])
  p1xx1 <- sum(p[1,,,1])
  p0xx0 <- sum(p[2,,,2])
  px11x <- sum(p[,1,1,])
  px00x <- sum(p[,2,2,])
  px1x1 <- sum(p[,1,,1])
  px0x0 <- sum(p[,2,,2])
  pxx11 <- sum(p[,,1,1])
  pxx00 <- sum(p[,,2,2])
  
  p1xxx <- sum(p[1,,,])
  px1xx <- sum(p[,1,,])
  pxx1x <- sum(p[,,1,])
  pxxx1 <- sum(p[,,,1])
  
  pa <- ((p11xx+p00xx) + (p1x1x+p0x0x) + (p1xx1+p0xx0) + (px11x+px00x) + (px1x1+px0x0) + (pxx11+pxx00))/6
  pi1 <- p1xxx
  pi2 <- px1xx
  pi3 <- pxx1x
  pi4 <- pxxx1
  pe2 <- (pi1*pi2+(1-pi1)*(1-pi2) + pi1*pi3+(1-pi1)*(1-pi3) + pi1*pi4+(1-pi1)*(1-pi4) + 
            pi2*pi3+(1-pi2)*(1-pi3) + pi2*pi4+(1-pi2)*(1-pi4) + pi3*pi4+(1-pi3)*(1-pi4))/6
  
  return((pa-pe2)/(1-pe2))
}

calc_gwet_fromp_R4 <- function(p){
  p11xx <- sum(p[1,1,,])
  p00xx <- sum(p[2,2,,])
  p1x1x <- sum(p[1,,1,])
  p0x0x <- sum(p[2,,2,])
  p1xx1 <- sum(p[1,,,1])
  p0xx0 <- sum(p[2,,,2])
  px11x <- sum(p[,1,1,])
  px00x <- sum(p[,2,2,])
  px1x1 <- sum(p[,1,,1])
  px0x0 <- sum(p[,2,,2])
  pxx11 <- sum(p[,,1,1])
  pxx00 <- sum(p[,,2,2])
  
  p1xxx <- sum(p[1,,,])
  px1xx <- sum(p[,1,,])
  pxx1x <- sum(p[,,1,])
  pxxx1 <- sum(p[,,,1])
  
  pa <- ((p11xx+p00xx) + (p1x1x+p0x0x) + (p1xx1+p0xx0) + (px11x+px00x) + (px1x1+px0x0) + (pxx11+pxx00))/6
  pi <- (p1xxx + px1xx + pxx1x + pxxx1)/4
  pe <- 2*pi*(1-pi)
  
  return((pa-pe)/(1-pe))
}
