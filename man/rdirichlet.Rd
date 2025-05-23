% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculator_f.R
\name{rdirichlet}
\alias{rdirichlet}
\title{Draw from a dirichlet distribution based on number of counts in transition. Adapted from brms::rdirichlet}
\usage{
rdirichlet(n = 1, alpha, seed = NULL)
}
\arguments{
\item{n}{Number of draws (must be >= 1). If n>1, it will return a list of matrices.}

\item{alpha}{A matrix of alphas (transition counts)}

\item{seed}{An integer which will be used to set the seed for this draw.}
}
\value{
A transition matrix. If n>1, it will return a list of matrices.
}
\description{
Draw from a dirichlet distribution based on number of counts in transition. Adapted from brms::rdirichlet
}
\examples{
rdirichlet(n=1,alpha= matrix(c(1251, 0, 350, 731),2,2))
rdirichlet(n=2,alpha= matrix(c(1251, 0, 350, 731),2,2))
}
