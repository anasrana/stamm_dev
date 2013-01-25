remove(list=ls())
load('../data/sandraDat.rdat')
source('func_rust.r')
tmp <- rust.clst.fit(gData=dat$g[dat$ind,], tData=dat$t, lambda=0.1, n.states=3, fit.as='log2Dat', rSmpls=c(20))

save(tmp, file='test.rdat')
