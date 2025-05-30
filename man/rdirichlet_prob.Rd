% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculator_f.R
\name{rdirichlet_prob}
\alias{rdirichlet_prob}
\title{Draw from a dirichlet distribution based on mean transition probabilities and standard errors}
\usage{
rdirichlet_prob(n = 1, alpha, se, seed = NULL)
}
\arguments{
\item{n}{Number of draws (must be >= 1). If n>1, it will return a list of matrices.}

\item{alpha}{A matrix of transition probabilities}

\item{se}{A matrix of standard errors}

\item{seed}{An integer which will be used to set the seed for this draw.}
}
\value{
A transition matrix. If n>1, it will return a list of matrices.
}
\description{
Draw from a dirichlet distribution based on mean transition probabilities and standard errors
}
\examples{
rdirichlet_prob(n=1,alpha= matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7),3,3),
se=matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7)/10,3,3))

rdirichlet_prob(n=2,alpha= matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7),3,3),
se=matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7)/10,3,3))
}
