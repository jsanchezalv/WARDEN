% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculator_f.R
\name{qcond_weibullPH}
\alias{qcond_weibullPH}
\title{Conditional quantile function for WeibullPH (flexsurv)}
\usage{
qcond_weibullPH(rnd = 0.5, shape, scale, lower_bound = 0)
}
\arguments{
\item{rnd}{Vector of quantiles (between 0 and 1)}

\item{shape}{Shape parameter of WeibullPH}

\item{scale}{Scale (rate) parameter of WeibullPH (i.e., as in hazard = scale * t^(shape - 1))}

\item{lower_bound}{Lower bound (current time)}
}
\value{
Estimate(s) from the conditional weibullPH distribution based on given parameters
}
\description{
Conditional quantile function for WeibullPH (flexsurv)
}
\examples{
qcond_weibullPH(rnd = 0.5, shape = 2, scale = 0.01, lower_bound = 5)
}
