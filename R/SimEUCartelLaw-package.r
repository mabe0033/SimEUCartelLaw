#' Simulation of Legal Exemption System for European Cartel Law
#'
#' SimEUCartelLaw implements simulation methods for the legal exemption system
#' fot the European cartel law.
#'
#' SimEUCartelLaw implements Monte Carlo simulations of a game-theoretic model 
#' for the legal exemption system of the European cartel in order to 
#' estimate the (mean) deterrent effect of this system.
#' The input and output parameters of the simulated cartel opportunities 
#' can be visualized by three-dimensional projections.
#'
#' @rdname SimEUCartelLawpackage
#' @docType package
#' @name SimEUCartelLaw
#' @useDynLib SimEUCartelLaw, .registration = TRUE, .fixes = "C_"
#' @examples
#' Par <- list(Phi=c(0.1,0.5), Rho=c(0.5,0.9), Ksi=c(0.05,0.3), Chi=c(0.1,0.4),
#'             M=c(0.2,1.2), G=c(0.05,0.2), A=c(0.1,0.3))
#' res <- LEgame(params=Par,m=100000)
#' print(aggResults(res))
#' print(CorrStudy(params=Par, m=10000))
#' print(CorrStudySplit(params=Par, m=10000))
#' \donttest{
#' RglPlot(LEgame(params=Par, m=10000))
#' NoRglPlot(LEgame(params=Par, m=10000))
#' }
#'
NULL
