## with the analysis.dat
gamma.val = 2.5
n.ef = 13
p.val = rep(NA, n.ef) # save the p-values
caliper.val = 1.0

## for each proposed IV: 
for(k in 1:n.ef){
  
  if(k == 1) match.out = tryCatch(matchit(X1 ~ W, method = "nearest", 
                                          exact = c("X2"), data = analysis.dat, options = (warn = -1)), error=function(err) NA) 
  
  if(k>1 & k < n.ef-1) match.out = tryCatch(matchit(eval(parse(text=paste0("X", k))) ~ W, method = "nearest", 
                                                    exact = c(paste0("X", k-1), paste0("X", k+1)), 
                                                    data = analysis.dat, options = (warn = -1)), error=function(err) NA) 
  if(k == n.ef-1) match.out = tryCatch(matchit(eval(parse(text=paste0("X", n.ef-1))) ~ W, method = "nearest", 
                                               exact = c(paste0("X", n.ef-2)), 
                                               data = analysis.dat, options = (warn = -1)), error=function(err) NA) 
  
  if(k==n.ef) match.out = tryCatch(matchit(eval(parse(text=paste0("X", n.ef))) ~ W, method = "nearest", std.caliper = FALSE, caliper = caliper.val, 
                                           exact = c("C"), data = analysis.dat, options = (warn = -1)), error=function(err) NA); 
  
  
  outcomes = treatments = stratum = c()
  index = 0
  if(class(match.out) == "matchit"){
    nj = length(table(match.out$subclass))
    for(j in 1:nj){
      index = index + 1
      treated = analysis.dat$res.Y[which(match.out$subclass == j & eval(parse(text=paste0("analysis.dat$X", k))) == 1)]
      control = analysis.dat$res.Y[which(match.out$subclass == j & eval(parse(text=paste0("analysis.dat$X", k))) == 0)]
      outcomes = c(outcomes, c(treated, control)) 
      stratum = c(stratum, rep(index, sum(match.out$subclass == j, na.rm = TRUE)))
      treatments = c(treatments, rep(1,length(treated)), rep(0, length(control))) 
    }
    test.ranks = hodgeslehmann(y = outcomes, z = treatments, st = stratum, align="hl") # aligned wilcoxon rank test (stratified)
    test.results = senstrat(sc = test.ranks, z = treatments, st = stratum, gamma = gamma.val)
    p.val[k] = test.results$Result["P-value"]
  }
}


mean(combine.p)
## Fisher's method
v = 9 # the minimum number of valid evidence factors
p.order = p.val[order(p.val)][(n.ef-v+1):n.ef]
truncatedP(p.order, trunc = 1)

pdf("Figure/realdata.pdf", width = 10, height = 7)
par(mfrow = c(1,1), cex.lab = 1.7, cex.axis = 1.5, 
    mar=c(5,5,3,2), oma = c(0,0,0,0), tcl = 0.5,  xpd = FALSE, cex.main = 2)
dot.color = ifelse(analysis.dat$D == 1, "red", "blue")
plot(analysis.dat$W, analysis.dat$Y, pch = 19, col = dot.color,
     cex = 0.5, ylab = "Baccalaureate exam scores (Y)", xlab = "Transition scores (W)")
legend("bottomright", c("Treated", "Controls"), pch = 19, cex = 1.5, bty = "n",
       col = c("red", "blue"))
abline(v = analysis.dat$C, lwd = 1.5, lty = 2)
dev.off()


## Figure with analysisdat ##
summary(lm(res.Y ~ C, data = analysis.dat))
dot.color = ifelse(analysis.dat$D == 1, "red", "blue")

pdf("Figure/real_standardized.pdf", width = 10, height = 7)
par(mfrow = c(1,1), cex.lab = 1.7, cex.axis = 1.5, 
    mar=c(5,5,3,2), oma = c(0,0,0,0), tcl = 0.5,  xpd = FALSE, cex.main = 2)
plot(analysis.dat$W - analysis.dat$C, analysis.dat$res.Y, pch = 19, col = dot.color,
     cex = 0.5, ylab = "Transformed outcomes (Y)", xlab = "Scores - cutoff (W-C)")
legend("bottomright", c("Treated", "Controls"), pch = 19, cex = 1.5, bty = "n",
       col = c("red", "blue"))
abline(v = 0, lwd = 2)
dev.off()


