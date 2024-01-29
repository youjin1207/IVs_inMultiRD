library(foreign)
library(tidyverse)

##### load original data
data1 = read.dta("Data/data-AER-1.dta") # 2001-2003 cohort administrative data (multiples by school-cutoff)

##### decode town identification and transition score
data1$ua = unlist(strsplit(data1$uazY,split = "z"))[seq(1,7243810,2)] # extract Town ID
data1$grade = data1$dzag + data1$zga # extract transition score 

##### Towns with 3 schools in 2001
data1_2001 = data1[data1$Y == 1,]
data1_2001_nusua3 = data1_2001 %>%
  filter(nusua == 3 & Y == 1) %>%
  distinct(sid2,.keep_all = T) 


##### reformat the data into RD design analysis
towns_cutoff = data1_2001 %>%
  filter(nusua == 3) %>%
  distinct(z,ua,.keep_all = T) %>%
  select(z,zga,ua) %>%
  group_by(ua) %>%
  summarise(c = max(zga),
            z = which.max(zga)) # extract the maximum cutoff in each town

colnames(data1_2001_nusua3)[c(1,2,7,8,9,11,18,20,21)] = c("C_index","W_normalized","Y","D","Year","C","Students","Towns","W")
data1_2001_nusua3 = data1_2001_nusua3 %>%
  select(c("C_index","W_normalized","Y","D","C","Students","Towns","W"))

for(i in 1:nrow(towns_cutoff)){
  data1_2001_nusua3$C[data1_2001_nusua3$Towns == towns_cutoff$ua[i]] = towns_cutoff$c[i]
  data1_2001_nusua3$D[data1_2001_nusua3$Towns == towns_cutoff$ua[i]] = (data1_2001_nusua3$D[data1_2001_nusua3$Towns == towns_cutoff$ua[i]] == towns_cutoff$z[i] + 1)
}

data1_2001_nusua3 = data1_2001_nusua3 %>%
  select(-c(C_index,W_normalized)) %>%
  distinct(Students,.keep_all = T)

##### Data Dictionary
# Y: Baccalaureate exam score
# D: Indicator of entering the school above cutoff (i.e. the better school)
# C: Cutoff
# Students: Student ID
# Towns: Towns ID
# W: Transition score 

original.dat = data1_2001_nusua3
## exclude the towns with two-sided compliance
tmp = as.numeric(names(table(original.dat$Towns)))
for(t in 1:21){
  subdat = original.dat[original.dat$Towns ==  tmp[t],]
  nrow(subdat)
  print(c(t, sum(subdat$D[subdat$W < subdat$C])))
}

## exclude towns 11, 12, 16, 20
subdat1 = original.dat[original.dat$Towns %in% tmp[-c(11,12,16,20)],]

sum(subdat1$D[subdat1$W < subdat1$C])
mean(subdat1$D[subdat1$W >= subdat1$C]) ## still fuzzy RD
length(unique(subdat1$C)) ## 17 cutoffs

## delete missing observations (temporarily)
mean(is.na(subdat1$Y))
subdat2 = na.omit(subdat1) # ignore missing outcomes

cutoffs = unique(subdat2$C); cutoffs = cutoffs[order(cutoffs)]
IVs = matrix(NA, nrow(subdat2), 17)
for(k in 1:16){
  for(i in 1:nrow(subdat2)){
    IVs[i,k] = as.integer(subdat2$C[i] <= cutoffs[17-k])
  }
}
IVs[,17] = as.integer(subdat2$W >= subdat2$C)
iv.corr = iv.pval = Y.pval = c()
for(k in 1:17){
  iv.corr[k] = cor(IVs[,k], subdat2$D)
  iv.pval[k] = summary(glm(subdat2$D ~ IVs[,k], family = binomial()))$coefficients[2,4]
  Y.pval[k] = summary(lm(subdat2$Y ~ IVs[,k]))$coefficients[2,4]
}



summary(glm(subdat2$D ~ IVs, family = binomial()))

## delete IV 3, 5, 11, 14
subdat3 = subdat2[!(subdat2$C %in% cutoffs[17-c(3,5,11,14)]), ]
cutoffs = unique(subdat3$C); cutoffs = cutoffs[order(cutoffs)]

IVs = matrix(NA, nrow(subdat3), 13)
for(k in 1:12){
  for(i in 1:nrow(subdat3)){
    IVs[i,k] = as.integer(subdat3$C[i] <= cutoffs[13-k])
  }
}
IVs[,13] = as.integer(subdat3$W >= subdat3$C)
iv.corr = iv.pval = Y.pval = c()
for(k in 1:13){
  iv.corr[k] = cor(IVs[,k], subdat3$D)
  iv.pval[k] = summary(glm(subdat3$D ~ IVs[,k], family = binomial()))$coefficients[2,4]
  Y.pval[k] = summary(lm(subdat3$Y ~ IVs[,k]))$coefficients[2,4]
}


summary(glm(subdat3$D ~ IVs, family = binomial()))
analysis.dat = data.frame(subdat3, IVs)
outcome.model = lm(Y ~ poly(W,3), data = analysis.dat)
analysis.dat$res.Y = outcome.model$residuals

save(analysis.dat, file = "Data/analysis.dat")
