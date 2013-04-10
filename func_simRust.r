library(stats)
## library(pracma)                       #library for lsqnonlin
## library(Matrix)                       #matrix specific functions
library(matrixStats)                  #use rowSds from library
#library(pls)                            #use stdize
library(BB)                             #use for spg
library(expm)


## ----------[ Simulations on single cells]--------------------
##' Forward simulation of model using single cells
##'
##' @title RustSim
##' @param nCells number of cells used in simulation (default 200)
##' @param nGenes number of genes simulated (default 9)
##' @param tau average jumpTime per state
##' @param nStates number of states simulated (default 2)
##' @param betaVec predefined beta values for simulation (default randomly generated)
##' @param dt timestep used throughout simulation
##' @param endTime final timepoint (in "real" time)
##' @param avNoise standard deviation of Gaussian noise added to the average trajectory (default 0.001)
##' @param sttNoise vector of standard deviations for Gaussian noise per state (default 0)
RustSim <- function(nCells = 200,nGenes = 9, tau = c(3.5,5,14.5), nStates = 2,
                    betaVec = NULL, dt = 0.01, endTime = 30, avNoise = 0.001,
                    sttNoise = rep(0, nStates)){

  ## Checking if some of the arguments passed are the right format and adjusting if possible
  if (!is.vector(betaVec) & !is.null(betaVec)) {
    stop('beta must be passed as vector')
  }
  if (length(sttNoise) != nStates & is.vector(sttNoise)) {
    sttNoise <- rep(sttNoise, nStates)
  }
  ## get jump time from the state occupation time $/tau$
  jumpTime <- NULL
  for (i in 1:(nStates - 1)) {
    jumpTime <- rbind(jumpTime, rexp(nCells, 1 / tau[i]))
  }
  ## Checks if beta values are passed in the argument (as vector)
  ## Assigns randomvalues to beta if not passed as argument
  if (is.null(betaVec)) {
    betaVec <- rnorm(nGenes * nStates, sd = 5)
    betaVec[betaVec<0]  <- 10^(-40)
  }
  ## Reshape betavector as matrix/
  betaVals <- matrix(betaVec, nGenes, nStates)
  ## Initialise results matrix with zeros
  gSim <- matrix(rep(0, nGenes * endTime / dt), nGenes, endTime / dt)
  nPoints <- endTime / dt
  for (iCells in 1:nCells) {
    gSimC <- NULL
    for (jStates in 1:(nStates - 1)) {
      nSentries <- ceiling(jumpTime[jStates, iCells]/dt)
      val_tmp <- betaVals[, rep.int(jStates, nSentries)]
      gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp), sd = sttNoise[jStates]) )
      rm(val_tmp)
    }
    nSentries <- (nPoints - ncol(gSimC))
    if (nSentries > 0) {
      val_fin <- betaVals[, rep.int(nStates, nSentries)]
      gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin),sd = sttNoise[nStates]))
      rm(val_fin)
    }
    gSim <- gSim + gSimC[, 1:nPoints] / nCells
  }
  ## Add gaussian noise to the data
  datasim <- log2(gSim) + matrix(rnorm(length(gSim), sd=avNoise), dim(gSim))
  dataSim <- 2^(datasim)
  ## Return values are full simulated data all time points, beta and if t given gSim
  ## with t-pts return all parameters used for this simulation
  return(list(gsim = gSim, beta = betaVals, dataSim = dataSim, n.cells = nCells, n.gns = nGenes,
              tau = tau, dt = dt, n.stt = nStates, ns.av = avNoise))
}

##'  Function that add noise to normal data
##'
##' @title AddNoise
##' @param sim
##' @param ns.sd
##' @param ns.type
##' @param gData simulated "expression" data
##' @return datasim
##' @author anas ahmad rana
AddNoise <- function(sim = NULL, ns.sd, ns.type, gData = NULL){
  if (is.null(gData) & exists('gsim', where = sim) )
    gData <- sim$gsim



  if (ns.type == 1) {
    datsim <- log2(gData) + abs(rowMeans(log(gData))) *
      matrix(rnorm(length(gData), sd = ns.sd), dim(gData))
    datasim <- 2^(datsim)
  } else if (ns.type==2)  {
    datsim <- log2(gData) + matrix(rnorm(length(gData), sd = ns.sd), dim(gData))
    datasim <- 2^(datsim)
  } else if (ns.type == 3)  {
    datasim <- gData + matrix(rnorm(length(gData), sd = ns.sd), dim(gData))
  }
  return(datasim)
}
