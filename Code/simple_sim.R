library(rdlocrand)
library(sensitivitymv)
library(senstrat)
library(MatchIt)
library(rddtools)

## set parameters
beta = 0.0
lambda = c(0,0,0,0)
gamma.val = 1
eta = c(0,0,0,0,5)
n = 100 # a total sample size
n.rep = 500 # the number of independent replicates 
n.ef = 5 # the number of proposed evidence factors
v = 4 # the minimum number of valid evidence factors
cutoffs = c(1,2,3,4,5)

p.val = matrix(NA, n.rep, n.ef)
combined.p = rep(NA, n.rep)

for(r in 1:n.rep){
  set.seed(r)
  
  W = runif(n, 0, 6)
  C = apply(rmultinom(n, 1, c(0.2, 0.2, 0.2, 0.2, 0.2)), 2, function(x) which(x==1))
  D = (W >= C)
  Z = U = matrix(NA, n, n.ef)

  for(k in 1:(n.ef-1)){
    Z[,k] = as.integer(C <= cutoffs[n.ef-k])
  }
  Z[,n.ef] = as.integer(W >= C)
  #D = as.integer(W >= C)
  D = ifelse(W >= C, rbinom(n, 1, 0.9-0.1*C), rep(0,n)) # consider a fuzzy RD
  for(k in 1:(n.ef)){
    U[,k] = as.integer(W >= cutoffs[k] & W < cutoffs[k] + 0.3)
  }
  
  Y0 = Z[,1:(n.ef-1)] %*% lambda + U %*% eta + 1*W + rnorm(n, 0, 1) 
  Y1 = beta + Y0
  obs.Y = Y0*(D==0) + Y1*(D==1)
  res.Y = lm(obs.Y ~ W)$residuals # use the transformed outcome
  dat = data.frame(obs.Y = obs.Y, res.Y = res.Y, Z = Z, C = C, W = W, D = D)
  
  ## for each proposed IV: 
  for(k in 1:n.ef){
    
    
    if(k == 1) match.out = tryCatch(matchit(Z.1 ~ W, method = "nearest", 
                                            exact = c("Z.2"), data = dat), error=function(err) NA) 
    
    if(k>1 & k < n.ef-1) match.out = tryCatch(matchit(eval(parse(text=paste0("Z.", k))) ~ W, method = "nearest", 
                                                      exact = c(paste0("Z.", k-1), paste0("Z.", k+1)), 
                                                      data = dat), error=function(err) NA) 
    if(k == n.ef-1) match.out = tryCatch(matchit(eval(parse(text=paste0("Z.", n.ef-1))) ~ W, method = "nearest", 
                                                 exact = c(paste0("Z.", n.ef-2)), 
                                                 data = dat), error=function(err) NA) 
    
    if(k==n.ef) match.out = tryCatch(matchit(eval(parse(text=paste0("Z.", n.ef))) ~ W, method = "nearest", std.caliper = FALSE, caliper = 0.1, exact = c("C"), data = dat), error=function(err) NA); 
    
    
    outcomes = treatments = stratum = c()
    index = 0
    if(class(match.out) == "matchit"){
      nj = length(table(match.out$subclass))
      for(j in 1:nj){
        index = index + 1
        if(k==n.ef){
          treated = dat$res.Y[which(match.out$subclass == j & eval(parse(text=paste0("dat$Z.", k))) == 1)]
          control = dat$res.Y[which(match.out$subclass == j & eval(parse(text=paste0("dat$Z.", k))) == 0)]
        }else{
          treated = dat$obs.Y[which(match.out$subclass == j & eval(parse(text=paste0("dat$Z.", k))) == 1)]
          control = dat$obs.Y[which(match.out$subclass == j & eval(parse(text=paste0("dat$Z.", k))) == 0)]
        }
        outcomes = c(outcomes, c(treated, control)) 
        stratum = c(stratum, rep(index, sum(match.out$subclass == j, na.rm = TRUE)))
        treatments = c(treatments, rep(1,length(treated)), rep(0, length(control))) 
      }
      test.ranks = hodgeslehmann(y = outcomes, z = treatments, st = stratum, align="hl") # aligned wilcoxon rank test (stratified)
      test.results = senstrat(sc = test.ranks, z = treatments, st = stratum, gamma = gamma.val)
      p.val[r,k] = test.results$Result["P-value"]
    }
  }
  
  p.order = p.val[r,][order(p.val[r,])][(n.ef-v+1):n.ef]
  combined.p[r] = tryCatch(truncatedP(p.order, 1), error=function(err) NA)
}

## rejection rates 
colMeans(p.val <= 0.05, na.rm = TRUE) # p-values from each evidence factor analysis
mean(combined.p <= 0.05, na.rm = TRUE) # a combined p-value
