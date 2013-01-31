remove(list=ls())
load('../data/sandraDat.rdat')
source('func_rust.r')

tmp <- rust.clst.fit(gData=dat$g[dat$ind,], tData=dat$t, lambda=0.05, n.states=c(2,3), fit.as='log2Dat', rSmpls=c(5,10), fit.pll=TRUE)
str(tmp)

save(tmp, file='test.rdat')


## Testing single gene fit
source('func_rust.r')
rust.fit.kStt(gData=g.tmp[1,], tData=dat$t, lambda = 0.01, n.states = 4, fit.as='log2Dat', fix.w=TRUE, w=w)


## Mclust for clustering
library(mclust)
tmp <- Mclust(data=dat$g[dat$ind[1:1000],], G=1:30, modelNames=c('EII', 'VII', 'VVI', 'VEI'))
plot(tmp)
