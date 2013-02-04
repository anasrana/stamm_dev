remove(list=ls())
load('../data/sandraDat.rdat')
source('func_rust.r')

tmp <- rust.clst.fit(gData=dat$g[dat$ind,], tData=dat$t, lambda=0.05, n.states=3, fit.as='log2Dat', rSmpls=5, fit.pll=F)
str(tmp)

save(tmp, file='test.rdat')
