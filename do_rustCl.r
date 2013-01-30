remove(list=ls())
load('../data/sandraDat.rdat')
source('func_rust.r')
tmp <- rust.clst.fit(gData=dat$g[dat$ind,], tData=dat$t, lambda=0.1, n.states=4, fit.as='log2Dat', rSmpls=c(20))
str(tmp)

save(tmp, file='test.rdat')

betas <- NULL
rss <- NULL
fit.gnes <- tmp
for(il in 1:length(fit.gnes)){
  betas <- rbind(betas, fit.gnes[[il]]$beta)
  rss <- c(rss, fit.gnes[[il]]$obj)
}

## Testing single gene fit
source('func_rust.r')
rust.fit.kStt(gData=g.tmp[1,], tData=dat$t, lambda = 0.01, n.states = 3, fit.as='log2Dat', fix.w=TRUE, w=w)
