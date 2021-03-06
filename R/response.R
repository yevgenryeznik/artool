#' Random response generator
#' 
#' Generates random efficacy responses of subjects randomized in a clinical trial
#' 
#' @author Yevgen Ryeznik (\email{yevgen.ryeznik@gmail.com}), Oleksandr Sverdlov
#' 
#' @param distr distribution of responses. The following values are supported:
#'    \itemize{
#'       \item \code{"binary"} -- \emph{Binary (Bernoulli)} distribution, 
#'       \item \code{"uniform"} -- \emph{Uniform} distribution, 
#'       \item \code{"normal"} -- \emph{Normal} distribution,
#'       \item \code{"exponential"} -- \emph{Exponential} distribution, 
#'       \item \code{"weibull"} -- \emph{Weibull} distribution, 
#'       \item \code{"loglogistic"} -- \emph{Log-logistic} distribution, 
#'       \item \code{"lognormal"} -- \emph{Log-normal} distribution.
#'    }
#' @param distr_param parameter(s) of the response distribution (a named list with numeric vector(s) 
#'    of length \eqn{K} each, where \eqn{K} is a total number of treatments studied in a clinical 
#'    trial. It must have the following structure depending on a value of the \code{distr} parameter:
#'    \itemize{
#'       \item \code{distr_param = list(prob)} if \code{distr == "binary"}; in this case response is
#'       a binary variable with values \eqn{{0, 1}}. It has a parameter \code{prob}\eqn{ = p}, and 
#'       \deqn{Pr(response = 1) = p,} 
#'       \item \code{distr_param = list(min, max)} if \code{distr == "uniform"}; the \emph{Unifrom} 
#'       distribution has a p.d.f. \deqn{f(x) = 1/(max-min),}
#'       \item \code{distr_param = list(mean, sd)} \code{"normal"}; it is assumed 
#'       that a p.d.f. of \emph{Normal} distribution with \code{mean}\eqn{ = \mu} and \code{sd}\eqn{ = \sigma}
#'       equals to \deqn{f(x) = 1/(\sigma\sqrt(2\pi))exp(-(x-\mu)^2/(2\sigma^2)),}
#'       \item \code{distr_param = list(rate)} if \code{distr == "exponential"}; it is assumed 
#'       that a p.d.f. of \emph{Exponential} distribution with \code{rate}\eqn{ = \lambda} equals to
#'       \deqn{f(x) = \lambda exp(-\lambda x),}  
#'       \item \code{distr_param = list(shape, scale)} if \code{distr == "weibull"}; it is assumed 
#'       that a p.d.f. of \emph{Weibull} distribution with \code{shape}\eqn{ = k} and \code{scale}\eqn{ = \lambda}
#'       equals to \deqn{f(x) = (k/\lambda)(x/k)^(k-1)exp(-(x/\lambda)^k),}
#'       \item \code{distr_param = list(shape, scale)} if \code{distr == "loglogistic"}; it is assumed 
#'       that a p.d.f. of \emph{Log-logistic} distribution with \code{shape}\eqn{ = \beta} and \code{scale}\eqn{ = \alpha}
#'       equals to \deqn{f(x) = (\beta/\alpha)(x/\alpha)^(\beta-1)/(1+(x/\alpha)^\beta)^2,}
#'       \item \code{distr_param = list(mu, sigma)} \code{"lognormal"}; it is assumed 
#'       that a p.d.f. of \emph{Log-normal} distribution with \code{mu}\eqn{ = \mu} and \code{sigma}\eqn{ = \sigma}
#'       equals to \deqn{f(x) = 1/(x\sigma\sqrt(2\pi))exp(-(lnx-\mu)^2/(2\sigma^2)).}
#'    }
#' @param treatment vector of treatment assignments (contains integers from 1 to \eqn{K}, where 
#'    \eqn{K} is a total number of treatments studied in a clinical trial).
#' 
#' @return a numeric vector of rsponses' values 
#' 
#' @examples
#'    library(artool)
#'    
#'    # generate random treatment assignments (4 treatments, equal probabilities)
#'    treatment <- sample(1:4, 100, replace = TRUE, prob = rep(0.25,4))
#'    
#'    # generate normal responses with the same means = 0 and all sds = 1
#'    response("normal", list(mean = rep(0,4), sd = rep(1, 4)), treatment)
#'    
#'    # generate normal responses with means c(0, 0.2, 0.4, 0.6), and sd's = c(0.1, 0.5, 0.7, 1)
#'    response("normal", list(mean = c(0, 0.2, 0.4, 0.6), sd = c(0.1, 0.5, 0.7, 1)), treatment)
#'    
#'    # generate Weibul responses with shapes = c(0.5, 1, 1.5, 2) and scale = 1
#'    response("weibull", list(shape = c(0.5, 1, 1.5, 2), scale = rep(1, 4)), treatment)
#'     
#' @useDynLib artool         
#' 
#' @export

response <- function(distr, distr_param, treatment) {
  .response(distr, distr_param, treatment)
}