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

plot.cv.lambda.rust <- function(dat.mat, lambda, x.lab='', y.lab=''){

  dat.l <- data.frame(y.val=apply(dat.mat, 1, sum), lambda=lambda)
  ggplot(dat.l, aes(x=lambda, y=y.val)) +
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



}
