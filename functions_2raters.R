#######################################
#### Functions for the R=2 methods ####
#######################################

## Function for calculating p11 from (p1, p2, f12)
calc_p11 <- function(p1, p2, f12){
  if (f12==0){
    p11 <- p1*p2
  } else {
    p11 <- (-1-(exp(f12)-1)*(p1+p2)+sqrt((1+(exp(f12)-1)*(p1+p2))^2+4*(1-exp(f12))*exp(f12)*p1*p2))/(2*(1-exp(f12)))
  }
  return(p11)
}

## Function for calculating TRAC F(a) from f12
calc_Fa_from_f12 <- function(f12, a=0.5){
  (exp(a*f12)-1)/(exp(a*f12)+1)
}

## Function for calculating f12 from TRAC F(a)
calc_f12_from_Fa <- function(Fa, a=0.5){
  1/a*log( (1+Fa)/(1-Fa) )
}

## Function for calculating f12's point estimate, given data after continuity correction
calc_f12est <- function(y, CC=0.5){
  y_CC <- y + CC
  y11 <- y_CC[1]
  y10 <- y_CC[2]
  y01 <- y_CC[3]
  y00 <- y_CC[4]
  ## f12
  f12 <- log(y11*y00/y10/y01)
  return(f12)
}

## Function for calculating f12's variance using the simple asymptotics (SA) method
calc_f12var_SA <- function(y, CC=0.5){
  y_CC <- y + CC
  y11 <- y_CC[1]
  y10 <- y_CC[2]
  y01 <- y_CC[3]
  y00 <- y_CC[4]
  ## f12 var
  f12var <- 1/y11 + 1/y00 + 1/y10 + 1/y01
  return(f12var)
}

## Function (and helper functions) for implementing the likelihood-ratio test-based (LRT) method for CI construction
log_L <- function(params, y) {
  f1 <- params[1]
  f2 <- params[2]
  f12 <- params[3]
  
  -sum(y) * log(1 + exp(f1) + exp(f2) + exp(f1 + f2 + f12)) + 
    f1 * (y[1] + y[2]) + 
    f2 * (y[1] + y[3]) + 
    f12 * y[1]
}
LRTI_function <- function(f12, y){ ## LRT invert, profile f1, f2
  obj_fun <- function(params) {
    f1 <- params[1]
    f2 <- params[2]
    log_L(c(f1, f2, f12), y)
  }
  # Optimize f1 and f2 while keeping f12 fixed
  result <- optim(par = c(0, 0), fn = obj_fun, method = "BFGS", control=list("fnscale"=-1))
  profile_loglik <- result$value
  
  # full-lik optim
  result_full <- optim(par = c(0,0,1), fn = log_L, y=y, method = "BFGS", control=list("fnscale"=-1))
  full_loglik <- result_full$value
  
  out <- 2 * (full_loglik - profile_loglik) - qchisq(0.95, df = 1)
  return(out)
}
calc_CI_LRTI <- function(y){
  LRTI_vectorized <- Vectorize(LRTI_function, "f12")
  roots_all <- uniroot.all(function(x) LRTI_vectorized(f12=x, y=y), c(-5, 10))
  lower <- roots_all[1]
  upper <- roots_all[2]
  return(c(lower, upper))
}

log_L_Fa <- function(params, y, a=0.5) {
  f1 <- params[1]
  f2 <- params[2]
  Fa <- params[3]
  
  f12 <- calc_f12_from_Fa(Fa=Fa, a=a)
  
  -sum(y) * log(1 + exp(f1) + exp(f2) + exp(f1 + f2 + f12)) + 
    f1 * (y[1] + y[2]) + 
    f2 * (y[1] + y[3]) + 
    f12 * y[1]
}

LRTI_function_Fa <- function(Fa, y, a=0.5){ ### LRT invert, profile f1, f2
  obj_fun <- function(params) {
    f1 <- params[1]
    f2 <- params[2]
    log_L_Fa(c(f1, f2, Fa), y, a)
  }
  # Optimize f1 and f2 while keeping f12 fixed
  result <- optim(par = c(0, 0), fn = obj_fun, method = "BFGS", control=list("fnscale"=-1))
  profile_loglik <- result$value
  
  # full-lik optim
  result_full <- optim(par = c(0,0,0.5), fn = log_L_Fa, y=y, method = "BFGS", control=list("fnscale"=-1))
  full_loglik <- result_full$value
  
  out <- 2 * (full_loglik - profile_loglik) - qchisq(0.95, df = 1)
  return(out)
}

calc_CI_LRTI_Fa <- function(y){
  LRTI_vectorized_Fa <- Vectorize(LRTI_function_Fa, "Fa")
  roots_all <- uniroot.all(function(x) LRTI_vectorized_Fa(Fa=x, y=y), c(-0.999999, 0.999999))
  
  lower <- roots_all[1]
  upper <- roots_all[2]
  return(c(lower, upper))
}


## Function (and helper functions) for implementing the goodness of fit test-based (GOF) method for CI construction
GOF_profile_function <- function(f12, y){ #helper function
  obj_fun <- function(params) {
    f1 <- params[1]
    f2 <- params[2]
    log_L(c(f1, f2, f12), y)
  }
  # Optimize f1 and f2 while keeping f12 fixed
  result <- optim(par = c(0, 0), fn = obj_fun, method = "BFGS", control=list("fnscale"=-1))
  f1 <- result$par[1]
  f2 <- result$par[2]
  p00 <- 1/(1+exp(f1)+exp(f2)+exp(f1+f2+f12))
  p10 <- exp(f1)*p00
  p01 <- exp(f2)*p00
  p11 <- exp(f1+f2+f12)*p00
  
  stat <- (y[1]-N*p11)^2/(N*p11)+(y[2]-N*p10)^2/(N*p10)+(y[3]-N*p01)^2/(N*p01)+(y[4]-N*p00)^2/(N*p00)
  
  out <- stat - qchisq(0.95, df = 1)
  return(out)
}

calc_CI_GOF_profile <- function(y){
  GOF_profile_vectorized <- Vectorize(GOF_profile_function, "f12")
  roots_all <- uniroot.all(function(x) GOF_profile_vectorized(f12=x, y=y), 
                           c(-5, 10))
  roots_all_rev <- rev(roots_all) 
  
  # use the increasing order to search the first + to - sign change
  for (i in 1:length(roots_all)){
    #i=3
    test_L <- GOF_profile_function(f12=roots_all[i] - 0.01, y=y) 
    test_R <- GOF_profile_function(f12=roots_all[i] + 0.01, y=y) 
    if (test_L>0 & test_R<0){
      lower <- roots_all[i]
      break
    }
  }
  for (i in 1:length(roots_all_rev)){
    #i=3
    test_L <- GOF_profile_function(f12=roots_all_rev[i] - 0.01, y=y) 
    test_R <- GOF_profile_function(f12=roots_all_rev[i] + 0.01, y=y) 
    if (test_L<0 & test_R>0){
      upper <- roots_all_rev[i]
      break
    }
  }
  return(c(lower, upper))
}

### GOF method for Fa
GOF_profile_function_Fa <- function(Fa, y, a=0.5){ #helper function
  obj_fun <- function(params) {
    f1 <- params[1]
    f2 <- params[2]
    log_L_Fa(c(f1, f2, Fa), y, a)
  }
  # Optimize f1 and f2 while keeping f12 fixed
  result <- optim(par = c(0, 0), fn = obj_fun, method = "BFGS", control=list("fnscale"=-1))
  f1 <- result$par[1]
  f2 <- result$par[2]
  f12 <- calc_f12_from_Fa(Fa=Fa, a=a)
  p00 <- 1/(1+exp(f1)+exp(f2)+exp(f1+f2+f12))
  p10 <- exp(f1)*p00
  p01 <- exp(f2)*p00
  p11 <- exp(f1+f2+f12)*p00
  
  stat <- (y[1]-N*p11)^2/(N*p11)+(y[2]-N*p10)^2/(N*p10)+(y[3]-N*p01)^2/(N*p01)+(y[4]-N*p00)^2/(N*p00)
  
  out <- stat - qchisq(0.95, df = 1)
  return(out)
}

calc_CI_GOF_profile_Fa <- function(y){
  GOF_profile_vectorized_Fa <- Vectorize(GOF_profile_function_Fa, "Fa")
  roots_all <- uniroot.all(function(x) GOF_profile_vectorized_Fa(Fa=x, y=y), 
                           c(-0.999999, 0.999999))
  roots_all_rev <- rev(roots_all) 
  
  # use the increasing order to search the first + to - sign change
  for (i in 1:length(roots_all)){
    #i=3
    test_L <- GOF_profile_function_Fa(Fa=roots_all[i] - 0.01, y=y) 
    test_R <- GOF_profile_function_Fa(Fa=roots_all[i] + 0.01, y=y) 
    if (test_L>0 & test_R<0){
      lower <- roots_all[i]
      break
    }
  }
  for (i in 1:length(roots_all_rev)){
    #i=3
    test_L <- GOF_profile_function_Fa(Fa=roots_all_rev[i] - 0.01, y=y) 
    test_R <- GOF_profile_function_Fa(Fa=roots_all_rev[i] + 0.01, y=y) 
    if (test_L<0 & test_R>0){
      upper <- roots_all_rev[i]
      break
    }
  }
  return(c(lower, upper))
}


############################# Functions for implementing existing IRA methods ################################

calc_phi <- function(y){
  N=sum(y)
  p11 <- y[1]/N
  p10 <- y[2]/N
  p01 <- y[3]/N
  p00 <- y[4]/N
  
  p0x <- p00+p01
  p1x <- p11+p10
  px0 <- p00+p10
  px1 <- p11+p01
  
  phi <- (p11 * p00 - p10 * p01)/sqrt(p1x * p0x * px1 * px0)
  delta_vec <- c(p00/sqrt(p1x * p0x * px1 * px0) - (p1x+px1)/(2*p1x*px1)*phi,
                 -p01/sqrt(p1x * p0x * px1 * px0) - (p1x+px0)/(2*p1x*px0)*phi,
                 -p10/sqrt(p1x * p0x * px1 * px0) - (px1+p0x)/(2*px1*p0x)*phi,
                 p11/sqrt(p1x * p0x * px1 * px0) - (p0x+px0)/(2*p0x*px0)*phi)
  
  Sigma <- matrix(c(p11*(1-p11), -p11*p10, -p11*p01, -p11*p00,
                    -p11*p10, p10*(1-p10), -p10*p01, -p10*p00,
                    -p11*p01, -p10*p01, p01*(1-p01), -p01*p00,
                    -p11*p00, -p10*p00, -p01*p00, p00*(1-p00)),nrow=4, ncol=4)
  
  est_phi <- phi
  sigma2_phi <- t(delta_vec)%*%Sigma%*%delta_vec
  return( c(est_phi,
            sqrt(1/N*sigma2_phi)) )
}

calc_S <- function(y){
  pa <- (y[1]+y[4])/sum(y)
  S <- 2*pa - 1
  var.S <- 4/sum(y) * pa*(1-pa)
  return(c(S, sqrt(var.S)))
}

calc_kappa <- function(y){
  tb <- matrix(y, nrow=2, byrow = T)
  weights = diag(1, nrow=2)
  pa <- sum(weights * tb/sum(tb))
  pk. <- (tb %*% rep(1, 2))/sum(tb)
  p.l <- t((t(rep(1, 2)) %*% tb)/sum(tb))
  pe <- sum(weights * (pk. %*% t(p.l)))
  kappa <- (pa - pe)/(1 - pe)
  pkl <- tb/sum(tb)
  pb.k <- weights %*% p.l
  pbl. <- t(weights) %*% pk.
  sum1 <- 0
  for (k in 1:2) {
    for (l in 1:2) {
      sum1 <- sum1 + pkl[k, l] * (weights[k, l] - (1 - kappa) * (pb.k[k] + pbl.[l]))^2
    }
  }
  var.kappa <- (1/(sum(tb) * (1 - pe)^2)) * (sum1 - (pa - 2 * (1 - kappa) * pe)^2)
  return(c(kappa, sqrt(var.kappa)))
}

calc_AC1 <- function(y){
  tb <- matrix(y, nrow=2, byrow = T)
  weights = diag(1, nrow=2)
  pa <- sum(weights * tb/sum(tb))
  pk. <- (tb %*% rep(1, 2))/sum(tb)
  p.l <- t((t(rep(1, 2)) %*% tb)/sum(tb))
  pi.k <- (pk. + p.l)/2
  tw <- sum(weights)
  pe <- tw * sum(pi.k * (1 - pi.k))/(2 * (2 - 1))
  gwet.ac1 <- (pa - pe)/(1 - pe)
  pkl <- tb/sum(tb)
  sum1 <- 0
  for (k in 1:2) {
    for (l in 1:2) {
      sum1 <- sum1 + pkl[k, l] * (weights[k, l] - 2 * (1 - gwet.ac1) * tw * (1 - (pi.k[k] + pi.k[l])/2)/(2 * (2 - 1)))^2
    }
  }
  var.gwet <- (1/(sum(tb) * (1 - pe)^2)) * (sum1 - (pa - 2 * (1 - gwet.ac1) * pe)^2)
  return(c(gwet.ac1, sqrt(var.gwet)))
}

