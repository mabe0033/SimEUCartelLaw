#' Simulate the Legal Exemption Game
#' 
#' \code{LEgame} simulates the legal exemption game.
#'
#' \code{LEgame} simulates the deterrent effect of the European cartel law 
#' based on a game-theoretic model for the legal exemption system.
#'
#' @param params named list containing numeric vectors Phi, Rho, Chi, 
#' Ksi, M, G and A with the ranges for the input parameters.
#' @param m numeric scalar containing the number of Monte Carlo
#' replications. Defaults to \code{1e5}.
#' @param corrMat matrix containing the correlation matrix for the
#' simulation. Defaults to a 7x7 identity matrix.
#' @param QMC logical scalar. If \code{TRUE}, an equidistant grid is
#' generated, if \code{FALSE}, uniformly distributed random numbers are
#' simulated.
#' @param seed numeric scalar containing the random seed for the
#' simulation. Defaults to \code{1} in order to make results reproducible.
#' @return A dataframe  containing the realized output of the simulation.
#' @importFrom rgl decorate3d view3d rotationMatrix
#' @importFrom plot3D scatter3D
#' @importFrom plot3Drgl scatter3Drgl
#' @importFrom stats sd
#' @rdname SimEUCartelLawfunction
#' @examples
#' Par <- list(Phi=c(0.1,0.5), Rho=c(0.5,0.9), Ksi=c(0.05,0.3), Chi=c(0.1,0.4),
#'             M=c(0.2,1.2), G=c(0.05,0.2), A=c(0.1,0.3))
#' res <- LEgame(params=Par, m=100000)
#' print(aggResults(res))
#'
#' @export LEgame
LEgame <- function(params, m=1e5, corrMat=diag(7), QMC=FALSE, seed=1) {
  if (!is.null(seed)) set.seed(seed)                                            # set random seed (if available)
  if (QMC) {                                                                    # adjust number of (Q)MC runs for QMC applications
    nsteps <- round(m^(1/7))
    m      <- nsteps^7
  }
  corrMat <- 2*sin(pi*corrMat/6)                                                # prepare correlation matrix for simulation
  diag(corrMat) <- 1
  corrMat <- chol(corrMat)

  res <- .Call(C_Simulate, range(params$Rho), range(params$Phi),                # start simulation
               range(params$Ksi), range(params$Chi), range(params$M), 
               range(params$G), range(params$A), if (QMC) nsteps else m, 
               corrMat, QMC)

  dim(res) <- c(m,11)                                                           # set dimension and dimnames for result
  colnames(res) <- c("p1","p2","c","gg","rrho","rphi","rksi","rchi","rM","rG",
                     "rA")
  invisible(data.frame(res))                                                    # return results
}

#' Matrix containing the correlation structure
#' 
#' \code{corrStruct} contains the correlation structure of the input parameters.
#'
#' \code{corrStruct} contains the correlation structure of the input parameters.
#' The actual correlation matrix used in the simulation is calculated as the
#' corresponding identity maxtrix + r times this matrix.
#'
#' @export corrStruct
corrStruct  <- matrix(data=c(  0,   0,    0,    0,    0,   0,    0,             # correlation structure (apart from diagonal entries)
                               0,   0,   -1,    1,    1,   0,    1,
                               0,  -1,    0,   -1,   -1,   0,   -1,
                               0,   1,   -1,    0,    1,   0,    1,
                               0,   1,   -1,    1,    0,   0,    1,   
                               0,   0,    0,    0,    0,   0,    0,
                               0,   1,   -1,    1,    1,   0,    0
                               ),nrow=7,ncol=7)

#' Investigate the effect of correlated input parameters
#' 
#' \code{CorrStudy} investigates the effect of correlated input parameters
#'
#' \code{CorrStudy} performs repeated simulations via \code{LEgame} with
#' different values for the correlation intensity in order to illustrate
#' the effect of correlation on the deterrent effect of the legal exemption
#' system.
#'
#' @param params named list containing numeric vectors Phi, Rho, Chi, 
#' Ksi, M, G and A with the ranges for the input parameters.
#' @param m numeric scalar containing the number of Monte Carlo
#' replications (for each correlation intensity). Defaults to \code{1e5}.
#' @param rho a numeric vector containing correlation intensities. Defaults to
#' \code{seq(0.1,0.9,by=0.2)}.
#' @param QMC logical scalar. If \code{TRUE}, an equidistant grid is
#' generated, if \code{FALSE}, uniformly distributed random numbers are
#' simulated.
#' @param seed numeric scalar containing the random seed for each
#' simulation. Defaults to \code{1} in order to make results reproducible.
#' @return A matrix containing the results of the repeated simulations.
#' @examples
#' Par <- list(Phi=c(0.1,0.5), Rho=c(0.5,0.9), Ksi=c(0.05,0.3), Chi=c(0.1,0.4),
#'             M=c(0.2,1.2), G=c(0.05,0.2), A=c(0.1,0.3))
#' res <- CorrStudy(params=Par, m=10000)
#' print(res)
#' @export CorrStudy
CorrStudy <- function(params,m=1e5,rho=seq(0.1,0.9,by=0.2),QMC=FALSE,seed=1) {
  out <- matrix(0,nrow=11,ncol=length(rho))                                     # setup maxtrix for results 
  for (i in seq_along(rho)) {
    res <- LEgame(params,m,diag(7)+rho[i]*corrStruct,QMC,seed)                  # simulate with rho[i]
    out[,i] <- c(100*mean(res$gg==1), 100*mean(res$gg==2), 100*mean(res$gg==3), # collect results
                 100*mean(res$p1), 100*sd(res$p1), 100*mean(res$p2),
                 100*sd(res$p2), 100*mean(res$c), 100*sd(res$c),
                 100*mean((1-res$c)*res$rA), 100*sd((1-res$c)*res$rA))
  }               
  colnames(out) <- rho                                                          # set dimnames
  rownames(out) <- c("Fraction cases EQ 1", "Fraction cases EQ 3.1", 
                     "Fraction cases EQ 5", "$p_I$ (mean)", "$p_I$ (sd)", 
                      "$p_{II}$ (mean)", "$p_{II}$ (sd)",
                      "Compliance (mean)", "Compliance (sd)",
                      "Expected illegal gain (mean)", 
                      "Expected illegal gain (sd)")
  out                                                                           # return results
}

#' Investigate the effect of correlated input parameters depending on illegal 
#' gain
#' 
#' \code{CorrStudySplit} investigates the effect of correlated input parameters
#' and its dependence on the illegal gain \code{A}.
#'
#' \code{CorrStudySplit} performs repeated simulations via \code{LEgame} with
#' different values for the correlation intensity and reports results for 
#' compliance and expected illegal gain for various subsets of simulated
#' illegal gains \code{A} in order to further illustrate the effect of 
#' correlation on the deterrent effect of the legal exemption system.
#'
#' @param params named list containing numeric vectors Phi, Rho, Chi, 
#' Ksi, M, G and A with the ranges for the input parameters.
#' @param m numeric scalar containing the number of Monte Carlo
#' replications (for each correlation intensity). Defaults to \code{1e5}.
#' @param rho a numeric vector containing correlation intensities. Defaults to
#' \code{seq(0.1,0.9,by=0.2)}.
#' @param breaks a numeric vector with breaks for the construction of the 
#' intervals for the illegal gain \code{A}. Defaults to 
#' \code{seq(0.1,0.3,by=0.04)}.
#' @param QMC logical scalar. If \code{TRUE}, an equidistant grid is
#' generated, if \code{FALSE}, uniformly distributed random numbers are
#' simulated.
#' @param seed numeric scalar containing the random seed for each
#' simulation. Defaults to \code{1} in order to make results reproducible.
#' @return A matrix containing the results of the repeated simulations.
#' @examples
#' Par <- list(Phi=c(0.1,0.5), Rho=c(0.5,0.9), Ksi=c(0.05,0.3), Chi=c(0.1,0.4),
#'             M=c(0.2,1.2), G=c(0.05,0.2), A=c(0.1,0.3))
#' res <- CorrStudySplit(params=Par, m=10000)
#' print(res)
#' @export CorrStudySplit
CorrStudySplit <- function(params,m=1e5,rho=seq(0.1,0.9,by=0.2),
                           breaks=seq(0.1,0.3,by=0.04),QMC=FALSE,seed=1) {
  nb  <- length(breaks)-1
  out <- matrix(0,nrow=(nb+1)*2,ncol=length(rho))                               # setup maxtrix for results 
  for (i in seq_along(rho)) {
    res <- LEgame(params,m,diag(7)+rho[i]*corrStruct,QMC,seed)                  # simulate with rho[i]
    idx <- cut(res$rA,breaks=breaks)
    for (j in 1:nb) {
      out[j+nb+1,i] <- 100*mean(((1-res$c)*res$rA)[as.numeric(idx)==j])         # collect results
      out[j,i]      <- 100*mean(res$c[as.numeric(idx)==j])
    }
  out[2*(nb+1),i] <- 100*mean(((1-res$c)*res$rA))                               # overall results
  out[nb+1,i]     <- 100*mean(res$c)
  }
  colnames(out) <- rho                                                          # set dimnames
  rownames(out) <- paste(rep(c("Compliance (mean)","Exp. illegal gain (mean)"),
                             each=nb+1),
                         rep(c(levels(idx),"Overall"),times=2))
  out                                                                           # return results
}

#' Visualize results of simulation of legal exemption system
#' 
#' \code{RglPlot} visualizes the results of the simulation of the legal 
#' exemption system using 3D-projections and corresponding 3D-plots.
#'
#' \code{RglPlot} visualizes the results of the simulation of the legal 
#' exemption system using 3D-projections and corresponding 3D-plots using
#' \code{rgl}/\code{GL} to produce real 3D plots which can be rotated or
#' zoomed in/out by the user, even in browser windows via \code{WebGL}.
#'
#' @param res dataframe containing results of simulation using \code{LEgame}.
#' @param xvar character scalar containing variable for the x-axis.
#' Defaults to \code{"rA"}, the simulated illegal gain.
#' @param yvar character scalar containing variable for the y-axis.
#' Defaults to \code{"rM"}, the simulated fine.
#' @param zvar character scalar containing variable for the z-axis.
#' Defaults to \code{"c"}, the complicance level.
#' @param xf numeric scalar containing scaling constant for the x-axis. Defaults
#' to \code{1}.
#' @param yf numeric scalar containing scaling constant for the y-axis. Defaults
#' to \code{1}.
#' @param zf numeric scalar containing scaling constant for the z-axis. Defaults
#' to \code{1}.
#' @param userMatrix matrix containing information about the initial perspective
#' used for the plot. Defaults to \code{rotationMatrix(1.3,-1,0.28, 0.4)}.
#' @param fov numeric scalar containing the field-of-view angle in degrees.
#' Defaults to \code{30}.
#' @param zoom numeric scalar containing the zoom factor. Defaults to 
#' \code{0.95}.
#' @return Nothing useful, function called for its side effects.
#' @examples
#' \donttest{
#' Par <- list(Phi=c(0.1,0.5), Rho=c(0.5,0.9), Ksi=c(0.05,0.3), Chi=c(0.1,0.4),
#'             M=c(0.2,1.2), G=c(0.05,0.2), A=c(0.1,0.3))
#' RglPlot(LEgame(params=Par, m=10000))
#' }
#' @export RglPlot
RglPlot <- function(res,xvar="rA",yvar="rM",zvar="c",xf=1,yf=1,zf=1,
                    userMatrix=rotationMatrix(1.3,-1,0.28, 0.4),fov=30,
                    zoom=0.95) {
  res <- cbind(res,oc=(1-res[,"c"])*res[,"rA"])                                 # calculate expected illegal gain
  NiceNames <- c(rA="illegal gain",rM="fine",c="compliance",p1=expression(p[I]),# define nice labels
                 p2=expression(p[II]),rrho=expression(rho),rphi=expression(phi),
                 rksi=expression(xi),rG="G",rchi=expression(chi),
                 oc="expected illegal gain")
  scatter3Drgl(xf*res[,xvar],yf*res[,yvar],zf*res[,zvar],pch=".",               # create the plot
               xlab=NiceNames[xvar],ylab=NiceNames[yvar],main=NiceNames[zvar],
               zlab="")
  view3d(userMatrix=userMatrix,fov=fov,zoom=zoom)                               # start with good perspective
}

#' Visualize results of simulation of legal exemption system
#' 
#' \code{NoRglPlot} visualizes the results of the simulation of the legal 
#' exemption system using 3D-projections and corresponding 3D-plots.
#'
#' \code{NoRglPlot} visualizes the results of the simulation of the legal 
#' exemption system using 3D-projections and corresponding plots without 
#' using \code{rgl}/\code{GL}.
#'
#' @param res dataframe containing results of simulation using \code{LEgame}.
#' @param xvar character scalar containing variable for the x-axis.
#' Defaults to \code{"rA"}, the simulated illegal gain.
#' @param yvar character scalar containing variable for the y-axis.
#' Defaults to \code{"rM"}, the simulated fine.
#' @param zvar character scalar containing variable for the z-axis.
#' Defaults to \code{"c"}, the complicance level.
#' @param xf numeric scalar containing scaling constant for the x-axis. Defaults
#' to \code{1}.
#' @param yf numeric scalar containing scaling constant for the y-axis. Defaults
#' to \code{1}.
#' @param zf numeric scalar containing scaling constant for the z-axis. Defaults
#' to \code{1}.
#' @param pch numeric or character scalar containing the plot character used 
#' for the individual points. Defaults to \code{16}.
#' @param phi numeric scalar containing the phi angle (colatitude)
#' for the perspective in degrees. Defaults to \code{20}.
#' @param theta numeric scalar containing the theta angle (azimuthal direction)
#' for the perspective in degrees. Defaults to \code{-30}.
#' @param d numeric scalar for the strenth of the perspective effect. 
#' Defaults to \code{2}.
#' @return Nothing useful, function called for its side effects.
#' @examples
#' \donttest{
#' Par <- list(Phi=c(0.1,0.5), Rho=c(0.5,0.9), Ksi=c(0.05,0.3), Chi=c(0.1,0.4),
#'             M=c(0.2,1.2), G=c(0.05,0.2), A=c(0.1,0.3))
#' NoRglPlot(LEgame(params=Par, m=10000))
#' }
#' @export NoRglPlot
NoRglPlot <- function(res,xvar="rA",yvar="rM",zvar="c",xf=1,yf=1,zf=1,
                      pch=16,phi=20,theta=-30,d=2) {
  res <- cbind(res,oc=(1-res[,"c"])*res[,"rA"])                                 # calculate expected illegal gain
  NiceNames <- c(rA="illegal gain",rM="fine",c="compliance",p1=expression(p[I]),# define nice labels
                 p2=expression(p[II]),rrho=expression(rho),rphi=expression(phi),
                 rksi=expression(xi),rG="G",rchi=expression(chi),
                 oc="expected illegal gain")
  scatter3D(xf*res[,xvar],yf*res[,yvar],zf*res[,zvar],pch=pch,phi=phi,          # create the plot
            theta=theta,xlab=NiceNames[xvar],ylab=NiceNames[yvar],
            zlab=NiceNames[zvar],d=d,axes=TRUE,ticktype="detailed")
}

#' Aggregate results of the legal exemption game simulation
#'
#' \code{aggResults} aggregates the results of \code{LEgame}.
#'
#' \code{aggResults} aggregates the results of \code{LEgame} to a matrix
#' containing information about the fractions for the potential equilibria as
#' well as the means and standard deviations of the error probabilities, the
#' compliance level, and the expected illegal gains.
#' 
#' @param res dataframe containing results of simulation using \code{LEgame}.
#' @return A matrix containing the aggregated results.
#' @examples
#' Par <- list(Phi=c(0.1,0.5), Rho=c(0.5,0.9), Ksi=c(0.05,0.3), Chi=c(0.1,0.4),
#'             M=c(0.2,1.2), G=c(0.05,0.2), A=c(0.1,0.3))
#' res <- LEgame(params=Par, m=100000)
#' print(aggResults(res))
#'
#' @export aggResults
aggResults <- function(res) {
  fn <- function(x) 100*c(nrow(x)/nrow(res), mean(x$p1), sd(x$p1),              # function for calculation of columns
                          mean(x$p2), sd(x$p2), mean(x$c),sd(x$c),
                          mean((1-x$c)*x$rA), sd((1-x$c)*x$rA))
  out <- cbind(fn(res[res$gg==1,]), fn(res[res$gg==2,]), fn(res[res$gg==3,]),   # aggregate results
               fn(res))
  colnames(out) <- c("EQ 1", "EQ 3.1", "EQ 5", "Overall")                       # set dimnames
  rownames(out) <- c("Fraction cases", "$p_I$ (mean)", "$p_I$ (sd)", 
                      "$p_{II}$ (mean)", "$p_{II}$ (sd)",
                      "Compliance (mean)", "Compliance (sd)",
                      "Expected illegal gain (mean)", 
                      "Expected illegal gain (sd)")
  out                                                                           # return results
}
