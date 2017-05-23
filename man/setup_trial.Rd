\name{setup_trial}
\alias{setup_trial}
\title{Setup clinical trial (CT) for simulation. }
\description{
The function sets up CT's parameters for simulation. An output is a Rcpp object of a class Trial.
}
\usage{
setup_trial(w, nsbj, nsim, cohort, proc, proc_param)
}
\arguments{
w
}
\examples{
trial <- setup_trial(c(2,1,1), 100, 1000, 1, "DBCD", 2)
}
