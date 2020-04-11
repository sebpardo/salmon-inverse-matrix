###
# Functions used in the sample-importance-resampling analysis
# Sebastian Pardo 01/09/2018 
###

# This is the function that defines the structure of the matrix! Very important to define properly
matrix_5 <- function(P1, Prg, S1, S1p, S2, Sk, F1, F2, Fk) {
  A <- matrix(0, nrow = 5, ncol = 5, byrow = T)
  # Assigning vital rates to matrix 
  A[2, 1] <- P1                   # egg-to-smolt survival
  A[3, 2] <- Prg * S1      # transition probability from smolt to grilse (prop grilse * grilse survival * return * survival)
  A[4, 2] <- (1 - Prg) * S1p * S2   # transition probability from smolt to 2SW (prop 2SW * grilse survival * return * survival)
                              # this proportion should IDEALLY include the S1 estimate from the previous year ()
  A[5, 3] <- S2 
  A[5, 4] <- S2 
  A[5, 5] <- S2  # Aug 7: adding probability of kelts returning to breed (before we assumed they all died)
  A[1, 3] <- F1 * 0.5
  A[1, 4] <- F2 * 0.5
  A[1, 5] <- Fk * 0.5
  # demogR::eigen.analysis(A)
  A
}

# As above but including 3SW fish, which is needed for non-NL populations
matrix_6 <- function(P1, Prg, Pr2, S1, S1p, S1p2, S2, S2p, S3, Sk, F1, F2, F3, Fk) {
  A <- matrix(0, nrow = 6, ncol = 6, byrow = T)
  # Assigning vital rates to matrix 
  A[2, 1] <- P1                   # egg-to-smolt survival
  A[3, 2] <- Prg * S1      # transition probability from smolt to grilse (prop grilse * grilse survival * return * survival)
  A[4, 2] <- (1 - Prg) * Pr2 * S1p * S2   # transition probability from smolt to 2SW (prop 2SW * grilse survival * return * survival)
  # this proportion should IDEALLY include the S1 estimate from the previous year (S1p)
  A[5, 2] <- (1 - Prg) * (1 - Pr2) * S1p2 * S2p * S3  # transition probability from smolt to 2SW (prop 2SW * grilse survival * return * survival)
  # this proportion should IDEALLY include the S1 estimate from the previous year (S1p)
  A[6, 3] <- S2 
  A[6, 4] <- S2 
  A[6, 5] <- S2  # Aug 7: adding probability of kelts returning to breed (before we assumed they all died)
  A[1, 3] <- F1 * 0.5
  A[1, 4] <- F2 * 0.5
  A[1, 5] <- Fk * 0.5
  A[1, 6] <- Fk * 0.5
    # demogR::eigen.analysis(A)
  A
}


# As above but matrix is simplified to only include smolts and adults (1SW, 2SW, RS)
matrix_simple <- function(Prg, Prgp, S1, S1p, S2) {
  A <- matrix(0, nrow = 4, ncol = 4, byrow = T)
  # Assigning vital rates to matrix 
  A[2, 1] <- Prg * S1      # transition probability from smolt to grilse (prop grilse * grilse survival * return * survival)
  A[3, 1] <- (1 - Prgp) * S1p * S2   # transition probability from smolt to 2SW. This proportion has to include the S1 and Prg estimates from the previous year
  A[4, 2] <- S2 
  A[4, 3] <- S2 
  A[4, 4] <- S2  # Aug 7: adding probability of kelts returning to breed (before we assumed they all died)
  # demogR::eigen.analysis(A)
  A
}


# As above but matrix is simplified to only include smolts and adults (1SW, 2SW, RS)
matrix_simple <- function(Prg, Prgp, S1, S1p, S2) {
  A <- matrix(0, nrow = 4, ncol = 4, byrow = T)
  # Assigning vital rates to matrix 
  A[2, 1] <- Prg * S1      # transition probability from smolt to grilse (prop grilse * grilse survival * return * survival)
  A[3, 1] <- (1 - Prgp) * S1p * S2   # transition probability from smolt to 2SW. This proportion has to include the S1 and Prg estimates from the previous year
  A[4, 2] <- S2 
  A[4, 3] <- S2 
  A[4, 4] <- S2  # Aug 7: adding probability of kelts returning to breed (before we assumed they all died)
  # demogR::eigen.analysis(A)
  A
}

# As above but matrix is simplified to only include smolts and adults (1SW, 2SW, RS)
matrix_3 <- function(Prg, Prgp, S1, S1p, S2) {
  A <- matrix(0, nrow = 3, ncol = 3, byrow = T)
  # Assigning vital rates to matrix 
  A[2, 1] <- Prg * S1      # transition probability from smolt to grilse (prop grilse * grilse survival * return * survival)
  A[3, 1] <- (1 - Prgp) * S1p * S2   # transition probability from smolt to 2SW. This proportion has to include the S1 and Prg estimates from the previous year
  # A[4, 2] <- S2 
  # A[4, 3] <- S2 
  # A[4, 4] <- S2  # Aug 7: adding probability of kelts returning to breed (before we assumed they all died)
  # demogR::eigen.analysis(A)
  A
}




# This function transforms a population vector (as in those used in/resulting from matrix multiplication) 
# to small and large returns, which is what is needed for the cost function
vec2returns <- function(vec, grilse, SW2, SW3 = NA, kelt = NA) {
  smalls <- vec[grilse]
  larges <- sum(vec[na.omit(c(SW2, SW3, kelt))])
  c(smalls, larges)
}


generate_matrix_dist <- function(iter, P1, Prg, S1, S1p, S2, 
                                 Sk, F1, F2, Fk, dists, randomize = FALSE) {
  lengths <- Map(function(x) length(eval(parse(text = x))), dists) %>% unlist %>% as.numeric()
  if (!all(lengths %in% iter)) stop("Number or iterations and length of distribution vectors differ.")
  # if (length(S1) != iter | length(S2) != iter) stop("Number or iterations and length of distribution vectors differ.")
  # random order of vector
  # NOTE: DON'T USE RANDOMIZE!!! IT BREAKS THE SAMPLING OF THE EXTERNAL VECTORS!!!!
    if (randomize) {
    for (i in dists) {
      assign(i, sample(eval(parse(text = i))))
    }
  } 
  #matarray <- array(NA, dim = c(5, 5, iter))
  # Another option would be to use the replicate() function
  matlist <- Map(matrix_5, P1 = P1[], Prg = Prg[], S1 = S1[], S1p = S1p[],
                 S2 = S2[], 
                # Sr = Sr[], # Sk = Sk[],
                 F1 = F1[], F2 = F2[], Fk = Fk[])
  #return(abind::abind(matarray, along = 3)) # return as array
  matlist #return as list
}

generate_matrix_dist_6 <- function(iter, P1, Prg, Pr2, S1, S1p, S1p2,
                                   S2, S2p, S3,
                                 Sk, F1, F2, F3, Fk, dists) {
  lengths <- Map(function(x) length(eval(parse(text = x))), dists) %>% unlist %>% as.numeric()
  if (!all(lengths %in% iter)) stop("Number or iterations and length of distribution vectors differ.")
  #matarray <- array(NA, dim = c(5, 5, iter))
  # Another option would be to use the replicate() function
  matlist <- Map(matrix_6, P1 = P1[], Prg = Prg[], Pr2 = Pr2[],
                 S1 = S1[], S1p = S1p[], S1p2 = S1p2[],
                 S2 = S2[], S2p = S2p[],
                 S3 = S3[],
                 # Sr = Sr[], # Sk = Sk[],
                 F1 = F1[], F2 = F2[], F3 = F3[], Fk = Fk[])
  #return(abind::abind(matarray, along = 3)) # return as array
  matlist #return as list
}

generate_matrix_dist_simple <- function(iter, Prg, Prgp, Pr2, S1, S1p, S2, 
                                   dists) {
  lengths <- Map(function(x) length(eval(parse(text = x))), dists) %>% unlist %>% as.numeric()
  if (!all(lengths %in% iter)) stop("Number or iterations and length of distribution vectors differ.")
  #matarray <- array(NA, dim = c(5, 5, iter))
  # Another option would be to use the replicate() function
  matlist <- Map(matrix_simple, Prg = Prg[], Prgp = Prgp[], 
                 S1 = S1[], S1p = S1p[],
                 S2 = S2[]
                 # Sr = Sr[], # Sk = Sk[],
                 )
  #return(abind::abind(matarray, along = 3)) # return as array
  matlist #return as list
}


generate_matrix_dist_3 <- function(iter, Prg, Prgp, Pr2, S1, S1p, S2, 
                                        dists) {
  lengths <- Map(function(x) length(eval(parse(text = x))), dists) %>% unlist %>% as.numeric()
  if (!all(lengths %in% iter)) stop("Number or iterations and length of distribution vectors differ.")
  #matarray <- array(NA, dim = c(5, 5, iter))
  # Another option would be to use the replicate() function
  matlist <- Map(matrix_3, Prg = Prg[], Prgp = Prgp[], 
                 S1 = S1[], S1p = S1p[],
                 S2 = S2[]
                 # Sr = Sr[], # Sk = Sk[],
  )
  #return(abind::abind(matarray, along = 3)) # return as array
  matlist #return as list
}


cost_function <- function(returns_obs, returns_est, constant = 1) {
  # # weighting both small and large equally rather than adding the differences:
  # cost <- exp(sum(sqrt((returns_obs[1] - returns_est[1])^2))*-1/(sum(returns_obs[1]) * constant)) *
  #   exp(sum(sqrt((returns_obs[2] - returns_est[2])^2))*-1/(sum(returns_obs[2]) * constant))
  print(returns_obs)
  print(sd(returns_obs))
  print(sum(returns_obs))
  cost <- exp(sum(sqrt((returns_obs - returns_est)^2))*-1/(sum(returns_obs) * constant)) 
  # adapted from Smart et al. (the * -1 in the numerator)
  # print(sum(sqrt((returns_obs - returns_est)^2))/(sum(returns_obs) * constant))
  # print(sum(sqrt((returns_obs - returns_est)^2))*-1/(sum(returns_obs) * constant))
  # print(exp(sum(sqrt((returns_obs - returns_est)^2))*-1/(sum(returns_obs) * constant)))
  #cost <- sum(sqrt((returns_obs - returns_est)^2))/sum(returns_obs)
  # print(returns_obs)
  # print(returns_est)
  # print(Mod(returns_obs - returns_est))
  # print(cost)
  #print(paste("obs:", as.numeric(returns_obs), "est:", as.numeric(returns_est), "cost:", as.numeric(cost)))
  cost #ifelse(cost < 1, cost, 1)
} 

# After revisions, changing sum for sd of denominator
cost_function2 <- function(returns_obs, returns_est, returns_sd, constant = 1) {
  # # weighting both small and large equally rather than adding the differences:
  # cost <- exp(sum(sqrt((returns_obs[1] - returns_est[1])^2))*-1/(sum(returns_obs[1]) * constant)) *
  #   exp(sum(sqrt((returns_obs[2] - returns_est[2])^2))*-1/(sum(returns_obs[2]) * constant))
  cost <-  exp(-1 * sum(((returns_obs - returns_est)^2)/(returns_sd * constant))) 
  # adapted from Smart et al. (the * -1 in the numerator)
  # print(sum(sqrt((returns_obs - returns_est)^2))/(sum(returns_obs) * constant))
  # print(sum(sqrt((returns_obs - returns_est)^2))*-1/(sum(returns_obs) * constant))
  # print(exp(sum(sqrt((returns_obs - returns_est)^2))*-1/(sum(returns_obs) * constant)))
  #cost <- sum(sqrt((returns_obs - returns_est)^2))/sum(returns_obs)
  # print(returns_obs)
  # print(returns_est)
  # print(Mod(returns_obs - returns_est))
  # print(cost)
  #print(paste("obs:", as.numeric(returns_obs), "est:", as.numeric(returns_est), "cost:", as.numeric(cost)))
  cost #ifelse(cost < 1, cost, 1)
} 

testvec <- c(100000, 50000, 2000, 300, 100)


# NEED TO CHANGE THIS FUNCTION TO COMPARE LIFE STAGES (1SW, 2SW, RS) rather than clump into small and large salmon (!!!) 

apply_cost <- function(mats, curvec, pos.smalls = 2, pos.large = 3, pos.rs = NA, returnsobs, returns_sd, constant = 1) {
  # print("Matrix multiplying items like this:")
  # print(mats[[1]])
  # print("with items like this:")
  # print(curvec)
  nextvecs <- lapply(mats, function(x) x %*% curvec)
  #print(newvecs)# does the matrix multiplication for each matrix in the list
  # browser()
  estreturns <- lapply(nextvecs, function(x) c(x[pos.smalls, ], sum(x[pos.large, ]))) #, x[pos.rs, ])) # adding separate life stage for RS
  # str(head(estreturns)))
  # print(class(returnsobs))
  # print(class(returns_sd))
  costs <- lapply(estreturns, function(x) cost_function2(returns_obs = returnsobs, returns_est = x, returns_sd = returns_sd, constant = constant)) #0.5))
  list(vectors = nextvecs, costfun = costs) # negative log
}


testvec <- c(1000000, 10000, 2000, 300, 100)

# vec is the current year population vector for which is multiplied with each matrix to obtain the n + 1 vectors
apply_cost_smalls <- function(mats, curvec, pos.smalls = 2, returnsobs, constant = 1) {
  newvecs <- lapply(mats, function(x) x %*% curvec)
  #print(newvecs) # does the matrix multiplication for each matrix in the list
  #browser()
  estreturns <- lapply(newvecs, function(x) c(x[pos.smalls, ])) #, sum(x[pos.large, ]))) %>% # REMOVED LARGE SALMON FROM COST FUNCTION
  #print(estreturns)
  costs <- lapply(estreturns, function(x) cost_function(returns_obs = returnsobs, returns_est = x, constant = constant))
  list(vectors = newvecs, costfun = costs) # negative log
}

#apply_cost_smalls(mats[1:3], curvec = testvec, returnsobs = c(2350))

# lapply(newvecs, function(x) {
#   print(x[3, ])
#   print(x[c(4,5), ])
#   print(sum(x[c(4,5), ]))
#   c(x[3, ], sum(x[c(4,5), ]))
#   }) 

# This function performs the sample-importance-resampling
sir_matrix <- function(mats, iter_samples, curvec, returnsobs, weights.exp = 1, returns_sd, constant = 1) {
  # calculate n + 1 population vectors and cost value for each
  #print(curvec)
  cost <- apply_cost(mats, curvec, returnsobs = returnsobs, constant = constant, returns_sd = returns_sd) 
  #print(cost)
  weights <- unlist(cost$costfun) # group cost values
  #browser()
  #print(weights)
  resample_pos <- base::sample(iter, iter_samples, replace = TRUE, prob = weights^weights.exp)
  matarray <- abind::abind(mats, along = 3)
  list(matrices = matarray[, , resample_pos], 
       vectors = do.call(cbind, cost$vectors[resample_pos]), 
       pos = resample_pos, 
       cost = weights)
}



# This function runs the SIR analysis after 
# NOTE: inter and iter_samples need to be already specified 
run_sir <- function(simuldat, S1, S2, Prg, constant) {
  
  # Parameters to be estimated and sampled from distributions
  distparams <- c("S1", "S2", "Prg")
  
  popvecs <- simuldat %>%
    #rename(small = SW1_obs) %>%
    mutate(smolt = lag(smolt_obs)) %>%
    filter(year > 3) %>% # & year < 18) %>%
    select(year, smolt, SW1_obs, SW2_obs)
  
  # years for simulated data
  siryears <- head(popvecs$year, -2) # remove the last year
  
  siroutput <- list()
  posteriors <- list()
  
  returnssd <- popvecs %>% 
    select(SW1_obs, SW2_obs) %>% 
    summarise_all(sd, na.rm = TRUE) %>% as.numeric
  
  for (i in seq_along(siryears)) {
    #for (i in 1) {
    print(paste("Start of i =", i))
    if (i == 1)  { 
      S1p <- simuldat[simuldat$year == siryears[i] - 1, "S1_true"]  # S1.empirical[S1.empirical$year == 1990, "S1"] # based on empirical S1 estimate on the first year 
      Prgp <- simuldat[simuldat$year == siryears[i] - 1, "Pr_true"]
    } else {
      S1p <- median(posteriors[[i - 1]]$S1) # Estimating from previous year's SIR output 
      Prgp <- median(posteriors[[i - 1]]$Prg)  # Estimating from previous year's SIR output
    } 
    # print(paste("Previous year's S1 =", S1p))
    # print(paste("Previous year's Prg =", Prgp))
    
    mats <- generate_matrix_dist_3(iter, Prg=Prg, Prgp=Prgp, S1=S1, S1p=S1p, S2=S2, dists = distparams)
    curyear <- siryears[i]
    paste0("Current year: ", curyear) %>% print
    paste0("Popn vector for current year +1 (",curyear + 1,"): ", 
           paste(as.character(popvecs[popvecs$year == (curyear + 1), -1]), collapse = ", ")) %>% print
    paste0("Popn vector for current year +2 (",curyear + 2,"): ", 
           paste(as.character(popvecs[popvecs$year == (curyear + 2), -1]), collapse = ", ")) %>% print
    
    # vector of next year's returns (for 1SW)
    
    
    returnsp1 <- vec2returns(as.numeric(popvecs[popvecs$year == (curyear + 1), -1]), grilse = 2, SW2 = 3, kelt = NA)
    
    print(returnsp1)
    # vector of returns 2 years forward (for 2SW)
    returnsp2 <- vec2returns(as.numeric(popvecs[popvecs$year == (curyear + 2), -1]), grilse = 2, SW2 = 3, kelt = NA)
    obsreturns <- c(returnsp1[1], returnsp2[2])
    print(obsreturns)
    
    paste("Observed returns from", curyear + 1, "small and", curyear + 2, 
          "large salmon:", paste(obsreturns, collapse = ", ")) %>% print
    
    curvec <- as.numeric(popvecs[popvecs$year == curyear,-1])
    # print(popvecs[popvecs$year == curyear,])
    # print(curvec)
    # print(length(curvec))
    # 
    # print(mats[[1]])
    # print(as.numeric(popvecs[popvecs$year == curyear, -1]))
    # print(obsreturns)
    
    siroutput[[i]] <- sir_matrix(mats = mats, iter_samples = iter_samples, 
                                 curvec = curvec, 
                                 returnsobs = obsreturns, #obsreturns,
                                 weights.exp = 1, 
                                 returns_sd = returnssd,
                                 constant = constant)#0.12) # 1 - 0.5) 
    posit <- siroutput[[i]]$pos
    costs <- siroutput[[i]]$cost
    posteriors[[i]] <- data.frame(year = curyear, S1 = S1[posit], S2 = S2[posit], Prg = Prg[posit], 
                                  weights = costs[posit])
  }
  
  list(siroutput = siroutput, posteriors = posteriors, simuldat = simuldat)
}


# Show unique number of posterior draws and maximum number of individual draws drawn
sir_convergence <- function(out) {
  lapply(out$siroutput, function(x) plyr::count(x$pos) %>% 
           arrange(desc(freq)) %>% 
           summarize(max = max(freq), n = length(unique(x)))) %>% 
    bind_rows()
}



pl.beta <- function(a, b, asp = if(isLim) 1, ylim = if(isLim) c(0,1.1)) {
  if(isLim <- a == 0 || b == 0 || a == Inf || b == Inf) {
    eps <- 1e-10
    x <- c(0, eps, (1:7)/16, 1/2+c(-eps,0,eps), (9:15)/16, 1-eps, 1)
  } else {
    x <- seq(0, 1, length = 1025)
  }
  fx <- cbind(dbeta(x, a,b), pbeta(x, a,b), qbeta(x, a,b))
  f <- fx; f[fx == Inf] <- 1e100
  matplot(x, f, ylab="", type="l", ylim=ylim, asp=asp,
          main = sprintf("[dpq]beta(x, a=%g, b=%g)", a,b))
  abline(0,1,     col="gray", lty=3)
  abline(h = 0:1, col="gray", lty=3)
  legend("top", paste0(c("d","p","q"), "beta(x, a,b)"),
         col=1:3, lty=1:3, bty = "n")
  invisible(cbind(x, fx))
}

pl.gamma <- function(a, b, asp = if(isLim) 1, ylim = if(isLim) c(0,1.1)) {
  if(isLim <- a == 0 || b == 0 || a == Inf || b == Inf) {
    eps <- 1e-10
    x <- c(0, eps, (1:7)/16, 1/2+c(-eps,0,eps), (9:15)/16, 1-eps, 1)
  } else {
    x <- seq(0, 1, length = 1025)
  }
  fx <- cbind(dgamma(x, a,b), pgamma(x, a,b), qgamma(x, a,b))
  f <- fx; f[fx == Inf] <- 1e100
  matplot(x, f, ylab="", type="l",
          main = sprintf("[dpq]gamma(x, a=%g, b=%g)", a,b))
  abline(0,1,     col="gray", lty=3)
  abline(h = 0:1, col="gray", lty=3)
  legend("top", paste0(c("d","p","q"), "gamma(x, a,b)"),
         col=1:3, lty=1:3, bty = "n")
  invisible(cbind(x, fx))
}


