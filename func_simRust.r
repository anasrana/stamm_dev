library(stats)
library(BB)                             #use for spg
library(expm)


## ----------[ Simulations on single cells]--------------------
##' Forward simulation of model using single cells
##'
##' @title RustSim
##' @param n.cells number of cells used in simulation (default 200)
##' @param n.genes number of genes simulated (default 9)
##' @param tau average jump.time per state
##' @param n.states number of states simulated (default 2)
##' @param beta.vec predefined beta values for simulation (default randomly generated)
##' @param dt timestep used throughout simulation
##' @param end.time final timepoint (in "real" time)
##' @param av.noise standard deviation of Gaussian noise added to population average (default 0.01)
##' @param stt.noise vector of standard deviations for Gaussian noise per state (default 0)
RustSim <- function(n.cells = 200, n.genes = 9, tau = c(3.5, 5, 14.5), n.states = 2,
                    beta.vec = NULL, dt = 0.01, end.time = 30, av.noise = 0.01,
                    stt.noise = rep(0, n.states), jump.dist = 'exp', sd.par = NULL) {

  ## Checking if some of the arguments passed are the right format and adjusting if possible
  if (!is.vector(beta.vec) & !is.null(beta.vec)) {
    stop('beta must be passed as vector')
  }
  if (length(stt.noise) != n.states & is.vector(stt.noise)) {
    stt.noise <- rep(stt.noise, n.states)
  }
  ## get jump time from the state occupation time use function JumpTime
  if (jump.dist != 'exp' & is.null(sd.par))
    sd.par  <- end.time / n.states
  jump.time <- JumpTime(tau, n.states, n.cells, jump.dist = jump.dist, sd.par = sd.par)
  ## Checks if beta values are passed in the argument (as vector)
  ## Assigns randomvalues to beta if not passed as argument
  if (is.null(beta.vec)) {
    beta.vec <- rnorm(n.genes * n.states, sd = 5)
    beta.vec[beta.vec < 0]  <- 10^(-40)
  }
  ## Reshape betavector as matrix
  betaVals <- matrix(beta.vec, n.genes, n.states)
  ## Initialise results matrix with zeros
  gSim <- matrix(rep(0, n.genes * end.time / dt), n.genes, end.time / dt)
  nPoints <- end.time / dt
  for (iCells in 1:n.cells) {
    gSimC <- NULL
    for (jStates in 1:(n.states - 1)) {
      nSentries <- ceiling(jump.time[jStates, iCells] / dt)
      val_tmp <- betaVals[, rep.int(jStates, nSentries)]
      gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp), sd = stt.noise[jStates]) )
      rm(val_tmp)
    }
    nSentries <- (nPoints - ncol(gSimC))
    if (nSentries > 0) {
      val_fin <- betaVals[, rep.int(n.states, nSentries)]
      gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin), sd = stt.noise[n.states]))
      rm(val_fin)
    }
    gSim <- gSim + gSimC[, 1:nPoints] / n.cells
  }
  ## Add gaussian noise to the data
  datasim <- log2(gSim) + matrix(rnorm(length(gSim), sd = av.noise), dim(gSim))
  dataSim <- 2^(datasim)
  ## Return values are full simulated data all time points, beta and if t given gSim
  ## with t-pts return all parameters used for this simulation
  return(list(gsim = gSim, beta = betaVals, dataSim = dataSim, n.cells = n.cells, n.gns = n.genes,
              tau = tau, dt = dt, n.stt = n.states, ns.av = av.noise))
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
AddNoise <- function(sim = NULL, ns.sd, ns.type, gData = NULL) {
  if (is.null(gData) & exists('gsim', where = sim) )
    gData <- sim$gsim

  if (ns.type == 1) {
    datsim <- log2(gData) + abs(rowMeans(log(gData))) *
      matrix(rnorm(length(gData), sd = ns.sd), dim(gData))
    datasim <- 2^(datsim)
  } else if (ns.type == 2)  {
    datsim <- log2(gData) + matrix(rnorm(length(gData), sd = ns.sd), dim(gData))
    datasim <- 2^(datsim)
  } else if (ns.type == 3)  {
    datasim <- gData + matrix(rnorm(length(gData), sd = ns.sd), dim(gData))
  }
  return(datasim)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title JumpTime
##' @param tau average transition times
##' @param n.states number of states simulated
##' @param n.cells number of cells simulated
##' @param jump.dist Distribution of jump time (default = 'exp')
##' @param sd.par second parameter to use for non exponential distributions
##' @return jump.time
##' @author anas ahmad rana
JumpTime  <- function(tau, n.states, n.cells, jump.dist = 'exp', sd.par = NULL) {
  jump.time <- NULL
  if (jump.dist == 'exp') {
    for (i in 1:(n.states - 1)) {
      jump.time <- rbind(jump.time, rexp(n.cells, 1 / tau[i]))
    }
  } else if (jump.dist == 'gaus') {
    for (i in 1:(n.states -1)) {
      jump.time  <- rbind(jump.time, abs(rnorm(n = n.cells, mean = 1 / tau[i], sd = sd.par)))
    }
  }
  return(jump.time)
}
