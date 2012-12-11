setwd("~/complexity/phd/code/")
## ----------[loading libraries and function files]--------------------
source("func_rust.r")
source("func_plotting.r")

## ----------[ loading data ]--------------------
load("~/complexity/phd/Rdata/wranaData.Rdata")
load("~/complexity/phd/Rdata/wrana4sttFit.Rdata")
load('resultsSim/parametersSim3sttGD.rdat')
load('resultsSim/plotSet1.rdat')

## ----------[ running simulation ]--------------------
dt <- 0.01
tau <- c(3.5, 5, 14.5)
nCells <- 1000
nStates <- 3
sttNoise <- rep(0.0,nStates)
avNoise <- 0.3
t <- seq(1,30, 3)
tPoints <- t/dt

## Forward simulation and writing transition matrix in proper form
w <- matrix(rep(0,3*3),3,3)
w[2,1] <- 1/tau[1]; w[3,2] <- 1/tau[2]

simltn <- singelCellSML(nCells = nCells, nStates = 3, dt=dt,sttNoise=sttNoise, avNoise=avNoise)

remove(datasim)
datasim <- addNoise(simltn$gsim, 0.2, 1)

save(t, w, dt, nCells, betaVec, tPoints, simltn, avNoise, sttNoise, tau, nStates, 
     file = 'resultsSim/parametersSim3sttGD.rdat')

mdl <- rustksttFwd(w, simltn$beta, seq(dt,30,dt))

res <- fitnSttModel(gData=simltn$dataSim[,tPoints], tData=t, tuning=0.1)

res <- fitnSttModel(gData=datasim[,tPoints], tData=t, tuning=0.1)

remove(estpar, wFit, betaFit, mdl, mdlFit)
estpar <- res$par
wFit <- matrix(0, nStates, nStates)
diag(wFit[-1,]) <- estpar[1:(nStates - 1)]
lnx <- length(estpar)
betaFit <- matrix( estpar[-c(1:(nStates - 1))], (lnx - nStates + 1)/ nStates, nStates)
mdl <- rustksttFwd(w, simltn$beta, seq(dt,30,dt))
mdlFit <- rustksttFwd(wFit, betaFit, seq(dt,30,dt))
res

save(t, w, dt, datasim, nCells, simltn, avNoise, sttNoise, tau, wFit, betaFit, mdlFit, file = 'resultsSim/plotSet2.rdat')

plotRustTraj(t, dt, beta=simltn$beta, betaFit=betaFit)

betaFileName <- paste("plotsSim/beta_ScaledNs", avNoise, "_set2.pdf", sep="")

plotRustTraj(t, dt, datasim[,tPoints], simltn$gsim, dataFit=mdlFit$S, mdlLegend=0)

trajctFileName <- paste("plotsSim/traject_ntScaledNs", avNoise, "_set2.pdf", sep="")
if(file.exists(trajctFileName)){
  warning(paste("file: <", trajctFileName, "> already exists! Output not saved"))
} else {
  dev.copy(pdf,trajctFileName)
  dev.off()
}


plotRustTraj(t, dt, simltn$dataSim[,tPoints], simltn$gsim, dataFit=mdlFit$S, mdlLegend=0) 
plotRustTraj(t, dt, simltn$dataSim[,tPoints], dataMdl=mdl$S, dataFit=mdlFit$S, mdlLegend=1)
