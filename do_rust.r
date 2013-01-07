setwd("~/complexity/phd/projects/rust/code")
## ----------[loading libraries and function files]--------------------
source("func_rust.r")
source("func_plotting.r")

## ----------[ loading data ]--------------------
load('resultsSim/parametersSim3sttGD.rdat')
load('resultsSim/plotSet2.rdat')

## ----------[ running simulation ]--------------------
dt <- 0.01; tau <- c(3.5, 5)
nCells <- 1000; nStates <- 3; avNoise <- 0.3
sttNoise <- rep(0.0,nStates)
t <- seq(1,30, 3); tPoints <- t/dt
tAll <- seq(dt,30, dt)

## Forward simulation and writing transition matrix in proper form
w <- matrix(rep(0,3*3),3,3)
w[2,1] <- 1/tau[1]; w[3,2] <- 1/tau[2]

simltn <- singelCellSML(nCells = nCells, nStates = 3, dt=dt,sttNoise=sttNoise, avNoise=avNoise, betaVec=betaVec)

save(dt, tau, nCells, nStates, sttNoise, avNoise, t, w, betaVec , file = 'parameters.rdat')

remove(datasim)
datasim <- addNoise(simltn$gsim, 0.1, 1)

save(t, w, dt, nCells, betaVec, tPoints, simltn, avNoise, sttNoise, tau, nStates, 
     file = 'resultsSim/parametersSim3sttGD.rdat')

mdl <- rustksttFwd(w, matrix(betaVec, 9, nStates), t)

res <- fitnSttModel(gData=simltn$dataSim[,t/dt], tData=t, tuning=0.1)

remove(estpar, wFit, betaFit, mdl, mdlFit)
estpar <- res$par
wFit <- matrix(0, nStates, nStates)
diag(wFit[-1,]) <- estpar[1:(nStates - 1)]
lnx <- length(estpar)
betaFit <- matrix( estpar[-c(1:(nStates - 1))], (lnx - nStates + 1)/ nStates, nStates)
mdl <- rustksttFwd(w, simltn$beta, seq(dt,30,dt))
mdlFit <- rustksttFwd(wFit, betaFit,seq(dt,30,dt))
res

save(t, w, dt, datasim, nCells, simltn, avNoise, sttNoise, tau, wFit, betaFit, mdlFit, file = 'resultsSim/plotSet2.rdat')

plotRustTraj(t, dt, beta=simltn$beta, betaFit=betaFit)

betaFileName <- paste("plotsSim/beta_ScaledNs", avNoise, "_set2.pdf", sep="")

plotRustTraj(t, dt, datasim[,t/dt], simltn$gsim, dataFit=mdlFit$S, mdlLegend=0)

trajctFileName <- paste("plotsSim/traject_ntScaledNs", avNoise, "_set2.pdf", sep="")
if(file.exists(trajctFileName)){
  warning(paste("file: <", trajctFileName, "> already exists! Output not saved"))
} else {
  dev.copy(pdf,trajctFileName)
  dev.off()
}

lambdas <- seq(0,0.4, 0.01)
bic <- NULL
aic <- NULL
obj <- NULL
rss <- NULL
for(i in 1:6){
  tmp <- rust.fit.nStt(dat$g[cl.nms,], dat$t, lambda = 0.1, n.states = i, fit.as='log2Dat')
  bic <- rbind(bic, tmp$bic)
  aic <- rbind(aic, tmp$aic)
  obj <- rbind(obj, tmp$fit$objective)
  rss <- rbind(rss, tmp$rss)
}

plot(1:6, aic, t='b', ylim=c(min(aic,bic), max(aic,bic)), ylab='objective function', xlab='states')
lines(1:6, bic, t='b', col='red')

tmp <- rust.fit.nStt(dat$g[cl.nms,], dat$t, lambda = 0.1, n.states = 4, fit.as='log2Dat')

1/tmp$w
1/fit4$w
fit4$beta[cl.nms,]
tmp$beta

## k-means clustering and fitting
dat.g <- dat$g[dat$ind,]

tmp.obj <- 10^20
for(j in 1:100){
  tmp <- kmeans(dat.g, 6, iter.max=100)
  tmp.obj <- min(tmp.obj, tmp$tot.withinss)
  if(tmp.obj == tmp$tot.withinss){
    cl.kmns <- tmp
  }
}
rep.gns <- NULL
rep.gns <- sapply(as.list(as.data.frame(t(tmp$centers))),
                  function(x){
                    rep.gns <- cbind(rep.gns, which.min(apply(abs(dat.g-x), 1,sum)))
                  } )
par(mfrow=c(2,1))
matplot(t(dat$t), t(dat.g[rep.gns,]), t='l')

## Plot centroids
matplot(t(dat$t), t(cl.kmns$centers), t='l')

cl.bic <- NULL
cl.rss <- NULL
for(i in 1:4){
  cl.fit <- rust.fit.nStt(gData=dat.g[rep.gns,], tData=dat$t[-3], lambda=0.1, n.states=4, fit.as='log2Dat')
  cl.bic <- cbind(cl.bic, cl.fit$bic)
  cl.rss <- cbind(cl.rss, cl.fit$rss)
}
