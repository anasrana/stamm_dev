remove(list=ls())
load('../data/sandraDat.rdat')
source('func_rust.r')
tmp <- rust.clst.fit(gData=dat$g[dat$ind,], tData=dat$t, lambda=0.1, n.states=4, fit.as='log2Dat', rSmpls=c(20))
str(tmp)

## testing normalisation
raw.dat <- read.table("/Users/anas/Dropbox/complexity/phd/data/Sandra2012/gene_exp.txt", header=TRUE) 
g.raw <- as.matrix(raw.dat[2:21])
g.raw <- g.raw[,c(2,12,13:20,1,3:7, 8:11)]
rownames(g.raw) <- raw.dat$external_gene_id
colnames(g.raw) <- unlist(lapply(data.frame(strsplit(colnames(g.raw), '_'))[2,], as.character))

library(edgeR)
library(ggplot2)
nF <- NULL
ctof <- seq(0,0.45,0.05)
for(i in ctof){
  tmp <- calcNormFactors(g.raw, method='TMM', refColumn=1, logratioTrim=i)
  nF <- rbind(nF, tmp)
}

matplot(ctof, nF, t='l')



sa <- stack(as.data.frame(nF))
sa$cutof <- ctof
qplot(cutof, values, data=sa, group=ind, colour=ind, geom='line')

## Testing single gene fit
source('func_rust.r')
rust.fit.kStt(gData=g.tmp[1,], tData=dat$t, lambda = 0.01, n.states = 4, fit.as='log2Dat', fix.w=TRUE, w=w)


## Mclust for clustering
library(mclust)
tmp <- Mclust(data=dat$g[dat$ind[1:1000],], G=1:30, modelNames=c('EII', 'VII', 'VVI', 'VEI'))
plot(tmp)


## Check simulations
library(ggplot2)
source('func_simRust.r')

dt <- 0.01
n.stt <- 4
n.gn <- 10
betas <- rnbinom(n.stt*n.gn, 1, 0.3 )


tmp.sim <- rust.simSnglCl(nCells = 1000,nGenes = n.gn, tau = c(3.5,8,14.5), nStates = n.stt,
                          betaVec=betas, dt = dt, endTime = 30, avNoise=0.01)
str(tmp.sim)

sim1 <- tmp.sim
save(sim1,file = '../data/sim1.rdat')


sa <- stack(as.data.frame(t(tmp.sim$gsim)))
sa$t <- seq(1*dt, 3000*dt, dt)

qplot(t, values, data=sa, geom='line', colour=ind, group=ind)
