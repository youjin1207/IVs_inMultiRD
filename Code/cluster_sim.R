library(rdlocrand)
library(sensitivitymv)
library(senstrat)
library(MatchIt)
library(rddtools)

## set parameters
beta = 0.0 
lambda = c(0,0,0,0,0,0,0,0,0)
gamma.val = 1
eta = 0
nc = 10 # the number of clusters
nn = 1000 # a size of each cluster
n = nc*nn # a total sample size
n.rep = 100 # the number of independent replicates 
n.ef = 10 # the number of proposed evidence factors
v = 10

p.val = matrix(NA, n.rep, n.ef)
combined.p = rep(NA, n.rep)

for(r in 1:n.rep){
  set.seed(r)
  
  tmp.C = runif(nc, 0, 1) 
  C = rep(tmp.C, each = nn)
  #W = rnorm(n, -0.5*C, 1)
  W = runif(n, 0, 1)
  Z = matrix(NA, n, n.ef)
 
  cutoffs = tmp.C[order(tmp.C)]
  for(k in 1:(n.ef-1)){
    Z[,k] = as.integer(C <= cutoffs[n.ef-k])
  }
  Z[,n.ef] = as.integer(W >= C)
  D = as.integer(W >= C)
  #D = ifelse(W >= C, rbinom(n, 1, 0.9-0.6*C), 0) # consider a fuzzy a
  U = (W >= tmp.C[1] & W < tmp.C[1] + 0.01)
  
  Y0 = Z[,1:(n.ef-1)] %*% lambda + eta*U + 1*W + rnorm(n, 0, 1) 
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
        treated = dat$res.Y[which(match.out$subclass == j & eval(parse(text=paste0("dat$Z.", k))) == 1)]
        control = dat$res.Y[which(match.out$subclass == j & eval(parse(text=paste0("dat$Z.", k))) == 0)]
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
colMeans(p.val <= 0.05) # p-values from each evidence factor analysis
mean(combined.p <= 0.05) # a combined p-value
