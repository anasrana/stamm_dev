require(ggplot2)
require(grid)

##' Function to show beta values for sim next to beta values for fit
##'
##'
##'
##' @title rust.comp.B
##' @param b.sim beta matrix that was used in the simulation
##' @param b.fit beta matrix from fitting
##' @return
##' @author anas ahmad rana
plot.comp.beta.rust <- function(b.sim, b.fit){

  p <- nrow(b.sim)
  k <- ncol(b.sim)
  b.val <- data.frame(beta=c(as.vector(b.sim), as.vector(b.fit)), tp=factor(rep(c('sim', 'fit'), each=p*k)), stt=rep(rep(1:k, each=p),2), gn=factor(rep(c( (rownames(b.sim)),  (rownames(b.fit)), k))))

b.g <- ggplot(b.val) +
    geom_bar(aes(x=stt, y=beta, alpha=tp, fill=stt),
             position='dodge', stat='identity') +
          facet_wrap(~gn, scales='free_y') +
      scale_alpha_discrete(range=c(0.5,1)) +
        xlab('States')+
          ylab('beta value') +
                theme_bw() +
                labs(fill='States', alpha='') +
                  theme(legend.key.size=unit(0.3, 'cm'),
                        legend.text = element_text(size=10, face='bold'),
                        axis.title.x = element_text(face='bold', size=20),
                        axis.title.y = element_text(face='bold', size=20),
                        strip.text.x = element_text(size=12),
                        strip.background = element_rect(colour = NA),
                        axis.text.x = element_text(size=10),
                        axis.text.y = element_text(size=10))
  return(b.g)
}


##' Function to plot trajectories of data and fit on top of each other
##'
##'
##' @title rust.comp.traj
##' @param g.dat data used in fitting
##' @param t time points of g.dat
##' @param b.fit fitted beta values
##' @param w fitted w matrix
##' @return
##' @author anas ahmad rana
plot.comp.traj.rust <- function(g.dat, t, b.fit, w){


  p <- nrow(b.fit)
  t.fit <- seq(0,max(t),0.01)
  g.fit <- rust.kStt(w, b.fit, t=t.fit)
  fit.dat <- data.frame(g=as.vector(g.fit$y), gn=factor(rep( (rownames(b.fit)), length(t.fit))),
                        t=rep(t.fit, each=p))
  sim.dat <- data.frame(g=as.vector(g.dat), gn = factor(rep( (rownames(b.fit)), length(t))),
                        t=rep(t, each=p))

  t.g <- ggplot(sim.dat, aes(x=t, y=g)) +
    geom_point(size=1.5, col='darkgreen') +
      geom_line(size=0.2, col='darkgreen') +
        geom_line(data=fit.dat, aes(x=t, y=g), col='blue') +
        facet_wrap(~gn, scales='free_y') +
          theme_bw() +
            xlab('time') +
              ylab('gene expression') +
  theme(legend.key.size=unit(0.3, 'cm'),
  legend.text = element_text(size=10, face='bold'),
  axis.title.x = element_text(face='bold', size=20),
  axis.title.y = element_text(face='bold', size=20),
                  strip.text.x = element_text(size=12),
                  strip.background = element_rect(colour = NA),
  axis.text.x = element_text(size=10),
  axis.text.y = element_text(size=10))

  return(t.g)
}

plot.traj.rust <- function(g.dat, t){

  p <- nrow(g.dat)
  sim.dat <- data.frame(g=as.vector(g.dat), gn = factor(rep( (rownames(g.dat)), length(t))),
                        t=rep(t, each=p))

  t.g <- ggplot(sim.dat, aes(x=t, y=g)) +
    geom_point(size=1.5, col='darkgreen') +
      geom_line(size=0.2, col='darkgreen') +
        facet_wrap(~gn, scales='free_y') +
          theme_bw() +
            xlab('time') +
              ylab('gene expression') +
  theme(legend.key.size=unit(0.3, 'cm'),
  legend.text = element_text(size=10, face='bold'),
  axis.title.x = element_text(face='bold', size=20),
  axis.title.y = element_text(face='bold', size=20),
                  strip.text.x = element_text(size=12),
                  strip.background = element_rect(colour = NA),
  axis.text.x = element_text(size=10),
  axis.text.y = element_text(size=10))

  return(t.g)
}

plot.beta.rust <- function(b.sim){


  p <- nrow(b.sim)
  k <- ncol(b.sim)
  if(is.null(rownames(b.sim)))
        rownames(b.sim) <- 1:p
  b.val <- data.frame(beta=as.vector(b.sim), stt=factor(rep(1:k, each=p)), gn=factor(rep( (rownames(b.sim)), k)))

b.g <- ggplot(b.val) +
    geom_bar(aes(x=stt, y=beta, fill=(stt)),
             position='dodge', stat='identity') +
          facet_wrap(~gn, scales='free_y') +
      scale_alpha_discrete(range=c(0.5,1)) +
        xlab('States')+
          ylab('beta value') +
                theme_bw() +
                labs(fill='States', alpha='') +
                  theme(legend.key.size=unit(0.3, 'cm'),
                        legend.text = element_text(size=10, face='bold'),
                        axis.title.x = element_text(face='bold', size=20),
                        axis.title.y = element_text(face='bold', size=20),
                        strip.text.x = element_text(size=12),
                        strip.background = element_rect(colour = NA),
                        axis.text.x = element_text(size=10),
                        axis.text.y = element_text(size=10))
  return(b.g)
}


vplayout <- function(x, y)viewport(layout.pos.row = x, layout.pos.col =y)

plot.cv.lambda.rust <- function(dat.mat, lambda, x.lab='', y.lab='', n.run=1, l.sz=1.2){

  if(n.run==1){
    dat.l <- data.frame(y.val=apply(dat.mat, 1, sum), lambda=lambda)
    plot.p <- ggplot(dat.l, aes(x=lambda, y=y.val)) +
      geom_point(size=2) +
        geom_line(size=1.2) +
          xlab(x.lab) +
            ylab(y.lab) +
              theme_bw() +
                theme(legend.position = 'none',
                      axis.title.x = element_text(face='bold', size=20),
                      axis.title.y = element_text(face='bold', size=20),
                      axis.text.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      strip.background = element_rect(colour = NA),
                      plot.title = element_text(face='bold'))
  } else if(n.run>1){
    dat.l <- data.frame(y.val=as.vector(dat.mat), lambda=rep(lambda, ncol(dat.mat)),
                        nrun=as.factor(rep(1:ncol(dat.mat), each=nrow(dat.mat))))
    plot.p <- ggplot(dat.l, aes(x=lambda, y=y.val, col=nrun, group=nrun)) +
      geom_point(size=2) +
        geom_line(size=l.sz) +
          xlab(x.lab) +
            ylab(y.lab) +
              theme_bw() +
                theme(axis.title.x = element_text(face='bold', size=20),
                      axis.title.y = element_text(face='bold', size=20),
                      axis.text.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      strip.background = element_rect(colour = NA),
                      plot.title = element_text(face='bold'))
  }

  return(plot.p)
}

plot.cv.facet.lambda.rust <- function(dat.mat, lambda, x.lab='predicted t-point', y.lab='RSS', t.dat,
                                      title.g='RSS of t-pt knockout, facet lambda'){

  rss.tk <- data.frame(rss=as.vector(rss.mat), lambda=rep(lambda, length(t.dat)-1),
                       t.ko = rep(2:(length(t.dat) ), each=nrow(dat.mat)))
  plot.p <- ggplot(rss.tk, aes(y=rss)) +
  geom_point(aes(x=t.ko)) +
  geom_line(aes(x=t.ko)) +
  facet_wrap(~lambda) +
  xlab(x.lab) +
  ylab(y.lab) +
  scale_x_continuous(limits=c(2,15), breaks=seq(2,16,2) ) +
  ggtitle(title.g) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x = element_text(face='bold', size=20),
        axis.title.y = element_text(face='bold', size=20),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=10),
        strip.background = element_rect(colour = NA),
        plot.title = element_text(face='bold'))
  return(plot.p)
}

plot.cv.facet.t.rust <- function(dat.mat, lambda, x.lab='lambda', y.lab='RSS', t.dat,
                                      title.g='RSS of t-pt knockout, facet predicted t'){
  rss.tk <- data.frame(rss=as.vector(rss.mat), lambda=rep(lambda, length(t.dat)-1),
                       t.ko = rep(2:(length(t.dat) ), each=nrow(dat.mat)))

  plot.p <- ggplot(rss.tk, aes(y=rss)) +
    geom_point(aes(x=lambda)) +
      geom_line(aes(x=lambda)) +
        facet_wrap(~t.ko) +
          xlab(x.lab) +
            ylab(y.lab) +
              scale_x_continuous(limits=c(0,0.3), breaks=seq(0,0.29,0.05) ) +
                ggtitle(title.g) +
                  theme_bw() +
                    theme(legend.position = 'none',
                          axis.title.x = element_text(face='bold', size=20),
                          axis.title.y = element_text(face='bold', size=20),
                          axis.text.x = element_text(size=7),
                          axis.text.y = element_text(size=14),
                          strip.text.x = element_text(size=10),
                          strip.background = element_rect(colour = NA),
                          plot.title = element_text(face='bold'))
  return(plot.p)
}

plot.beta.scatter.rust <- function(beta.sc, beta.al, title.g='Scatter plot comparing beta values',
                                   x.lab, b.scl = 'log', lmbd.vec, n.stt=4, n.gn=12){
  if(b.scl=='log'){  #All the beta values below are shifted by one,
                     #the assumption is that they contain 0 values
    beta.dm <- data.frame(beta0=rep(beta.sc +1 , ncol(beta.al)), beta=as.vector(beta.al +1),
                          lambda = rep(lmbd.vec, each=nrow(beta.al)),
                          stt=as.factor(rep(rep(1:n.stt, each =n.gn),ncol(beta.al))),
                          gn=as.factor(rep(1:n.gn, n.stt*ncol(beta.al))) )

    beta.scl <- c(1, 10^seq(2, ceiling(max(log10(beta.al))), 2) +1)
    beta.lbl <- c(0, 10^seq(2, ceiling(max(log10(beta.al))), 2) )

    ggplot(beta.dm, aes(x=beta0, y=beta)) +
      geom_point(aes(colour=gn, size=stt)) +
        scale_colour_brewer(palette="Paired") +
          facet_wrap(~lambda) +
            coord_trans(x='log', y='log') +
              scale_x_continuous(breaks= beta.scl, label=beta.lbl ) +
                scale_y_continuous(breaks= beta.scl, label=beta.lbl ) +
                  scale_size_discrete(range=c(1,2.5)) +
                    xlab(x.lab) +
                      ylab('beta') +
                        ggtitle(title.g) +
                          theme_bw() +
                            theme(axis.title.x = element_text(face='bold', size=20),
                                  axis.title.y = element_text(face='bold', size=20),
                                  axis.text.x = element_text(size=7),
                                  axis.text.y = element_text(size=7),
                                  strip.text.x = element_text(size=10),
                                  strip.background = element_rect(colour = NA),
                                  plot.title = element_text(face='bold'))
  }
}

## ********************************************************************************
##   functions that perform very simple calculations but are not
##   directly related to fitting or simulating
##   ********************************************************************************


rust.cv.rss <- function(fit.file, n.stt=4, n.gn=12, t.ko=1){
  load(fit.file)
  load(paste('~/complexity/phd/projects/rust/results/5mar13/sim_varScl_4.rdat',sep='') )

  beta.f <- matrix(NA, n.stt*n.gn, length(fit))
  for(i in 1:length(fit)){
    beta.f[,i] <- as.vector(fit[[i]][[t.ko]]$beta)
  }

  w.f <- matrix(NA, n.stt*n.gn, length(fit))
  for(i in 1:length(fit)){
    w.f[,i] <- as.vector(fit[[i]][[t.ko]]$w)
  }


  w.sim <- matrix(0, nrow(fit[[1]][[1]]$w),ncol(fit[[1]][[1]]$w))
  diag(w.sim[-1,]) <- 1/sim.dat$sim$tau

  rss.mat <- matrix(NA, length(fit), length(t.dat)-1)
  beta <- matrix(NA, length(fit), length(t.dat)-1)
  w <- matrix(NA, length(fit), length(t.dat)-1)
  rss.mat <- matrix(NA, length(fit), length(t.dat)-1)
  for( i.l in 1:length(fit) ){
    for( i.t in 2:length(t.dat) ){
      beta[i.l, i.t-1] <- sum((fit[[i.l]][[i.t]]$beta - sim.dat$beta)^2)
      w[i.l, i.t-1] <- sum((fit[[i.l]][[i.t]]$w - w.sim)^2)
      bt <- fit[[i.l]][[i.t]]$beta
      w.fit <- fit[[i.l]][[i.t]]$w
      rep.fit <- rust.kStt(w.fit, bt,  t.dat)
      rss.mat[i.l, i.t-1] <- sum((log2(g.dat[,i.t] + 1) - log2(rep.fit$y[,i.t] +1) )^2 )
    }
  }

  return(list(rss = rss.mat, beta = beta, w = w, sim.dat=sim.dat, fit.dat = fit,
              lambda=lmbd.mat[,1], beta.fit=beta.f, w.fit=w.f))
}
