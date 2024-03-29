% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/switchNet.R
\name{switchNet}
\alias{switchNet}
\title{Once tau.hat is available, refit the model to obtain the precision matrices}
\usage{
switchNet(X, tau, lambda = NULL, cstar = c(1, 1))
}
\arguments{
\item{X}{A p by N data matrix}

\item{tau}{Estimate of tau}

\item{lambda}{Optimal tuning parameters for each regime}
}
\description{
Once tau.hat is available, refit the model to obtain the precision matrices
}
