
## HMM `families': generic specifiction of a HMM model with
## re-estimation rules

## note: demission needs to have state index (k) as first argument for
## convenience when calling outer.


##' Generate one of several pre-defined families for fitting a HMM to
##' count data.
##'
##' Two emission distributions are supported: Poisson and Negative
##' Binomial (with an additional size parameter). The number of states
##' are determined by the length of the \code{state.scale} parameter;
##' these essentially define the mean of the emission distribution, in
##' combination with the additional scale parameter \code{mu}. Their
##' interpretation depends on \code{states.free}.
##'
##' When \code{states.free = TRUE}, the mean of the emission
##' distribution for all states are allowed to vary freely. The
##' initial values of these means are given by \code{mu *
##' state.scale}, but subsequently, the parameters \code{mu} and
##' \code{state.scale} have no separate meaning.
##'
##' On the other hand, when \code{states.free = FALSE}, the
##' \code{state.scale} parameter is constrained to be fixed at its
##' initially specified value, and only \code{mu} is estimated. This
##' option is motivated by biological applications where state means
##' correspond to copy number, which are typically constrained to be
##' simple integer ratios.
##' 
##' Furthermore, if \code{mixture = TRUE}, then the mean of the
##' emission distribution is taken to be a convex combination of 2
##' (motivated as the \sQuote{normal} copy number in a biological
##' context) and the mean defined by \code{mu * state.scale}, with
##' proportions \code{1-p} and \code{p} respectively.
##' 
##' @title HMM Family with Poisson and Negative Binomial Emission
##'     Distributions
##' @param model Character string, \code{"poisson"} or \code{"nbinom"}
##' @param prob Initial estimate of the transition matrix. Can be
##'     omitted, in which case the default is to assume that all
##'     transitions are equally likely, which is usually fine as an
##'     initial estimate. If explicitly specified, cells with 0
##'     probability will remain unchanged; this enables selected
##'     transitions to be disallowed.
##' @param mu Scalar numeric giving initial estimate of a scale factor
##'     controlling the means of the emission distribution. See
##'     details.
##' @param states.scale Numeric vector giving (relative) initial
##'     estimates of the mean of the emission distribution for each
##'     state. The actual mean is obtained by multiplying by
##'     \code{mu}. See details.
##' @param size The size parameter for the Negative Binomial model.
##' @param p A mixture parameter controlling proportion of the modeled
##'     population that can potentially have multiple states. See
##'     details.
##' @param states.free Logical flag specifying whether states can vary
##'     freely (otherwise they are relatively fixed at the initial
##'     estimates given by \code{states.scale}, and only \code{mu} can
##'     vary).
##' @param mixture Logical flag specifying whether the model should
##'     include a mixture parameter.
##' @return A list with components
##' \itemize{
##' \item \code{demission}: Function giving emission probabilities.
##' \item \code{pars}: list of parameters of the model (some possibly fixed). 
##' 
##' \item \code{updated.var}: Function computing the M-step, i.e.,
##' giving updated parameter values given current estimates.
##' 
##' \item \code{mean.states} Function giving vector of means of the
##' emission distribution for each state.
##'
##' }
##' 
##' @author Deepayan Sarkar
hmm.family <-
    function(model = c("poisson", "nbinom"),
             prob, mu, states.scale, size, p,
             states.free = FALSE, mixture = FALSE)
{
    model <- match.arg(model)
    len <- length(states.scale)
    if (missing(prob)) prob <- matrix(1 / len, len, len)
    if (model == "poisson")
    {
        if (states.free)
        {
            variable <-
                list(states = states.scale * mu,
                     prob = prob)
            fixed <- list()
            mean.states <- function(pars) pars$states
            demission <- function(k, x, pars, ...)
            {
                dpois(x, lambda = pars$states[k], log = TRUE)
            }
            updated.var <- function(xlist, Akl, Bk, Nkb, pars, ...)
            {
                ## needs to be checked
                ## print(Akl)
                list(states = Nkb / Bk, prob = Akl / rowSums(Akl))
            }
        }
        else
        {
            if (mixture) stop("Mixtures not yet supported with Poisson model")
            variable <- list(mu = mu, prob = prob)
            fixed <- list(states = states.scale)
            mean.states <- function(pars) pars$mu * pars$states
            demission <- function(k, x, pars, ...)
            {
                dpois(x, lambda = pars$mu * pars$states[k],
                      log = TRUE)
            }
            updated.var <- function(xlist, Akl, Bk, Nkb, pars, ...)
            {
                list(mu = sum(unlist(xlist)) / sum(pars$states * Bk),
                     prob = Akl / rowSums(Akl))
            }
        }
    }
    else ## negative binomial
    {
        if (states.free)
        {
            variable <-
                list(states = states.scale * mu,
                     size = size,
                     prob = prob)
            fixed <- list()
            mean.states <- function(pars) pars$states
            demission <- function(k, x, pars, ...)
            {
                dnbinom(x, mu = pars$states[k],
                        size = pars$size,
                        log = TRUE)
            }
            updated.var <-
                function(xlist, Akl, Bk, Nkb, pars,
                         control = list(), ...)
                {
                    states <- Nkb / Bk
                    xx <- unlist(xlist)
                    Q <- function(sigma, pars)
                    {
                        ans <- 
                            length(xx) * (sigma * log(sigma) - lgamma(sigma)) +
                                sum(lgamma(xx + sigma)) +
                                    sum(Nkb * (log(states) - log(sigma + states)) -
                                        sigma * Bk * log(sigma + states))
                        ##print(ans)
                        -ans ## to be minimized
                    }
                    res <-
                        nlminb(pars$size,
                               Q,
                               control = control,
                               lower = 1e-5, upper = Inf,
                               pars = pars)
                    if (res$convergence != 0)
                        warning(paste("nlminb did not converge:",
                                      res$message))
                    list(states = states,
                         size = res$par,
                         prob = Akl / rowSums(Akl))
            }
        }
        else # fixed states: option of normal+disease mixture
        {
            if (mixture) # normal + disease mixture
            {
                ## extra parameter p: proportion of disease cells
                ## (assumed homogeneous). Normal "copy number" is
                ## assumed to be 2 in this case
                variable <- list(mu = mu,
                                 size = size,
                                 p = p,
                                 prob = prob)
                fixed <- list(states = states.scale)
                mean.states <- function(pars) pars$mu * (2 * (1-pars$p) + pars$p * pars$states)
                demission <- function(k, x, pars, ...)
                {
                    dnbinom(x, mu = pars$mu * (2 * (1-pars$p) + pars$p * pars$states[k]),
                            size = pars$size,
                            log = TRUE)
                }
                updated.var <-
                    function(xlist, Akl, Bk, Nkb, pars,
                             control = list(), ...)
                {
                    xx <- unlist(xlist)
                    Q <- function(optpars, pars)
                    {
                        ## print(optpars)
                        mu <- optpars[1]
                        sigma <- optpars[2]
                        p <- optpars[3]
                        nbmean <- mu * (2 * (1-p) + p * pars$states)
                        ans <- 
                            length(xx) * (sigma * log(sigma) - lgamma(sigma)) +
                            sum(lgamma(xx + sigma)) +
                            sum(Nkb * (log(nbmean) -
                                       log(sigma + nbmean)) -
                                sigma * Bk * log(sigma + nbmean))
                        -ans ## to be minimized
                    }
                    res <-
                        nlminb(c(pars$mu, pars$size, pars$p),
                               Q,
                               control = control,
                               lower = 1e-5, upper = c(Inf, Inf, 1),
                               pars = pars)
                    if (res$convergence != 0)
                        warning(paste("nlminb did not converge:",
                                      res$message))
                    list(prob = Akl / rowSums(Akl),
                         mu = res$par[1],
                         size = res$par[2],
                         p = res$par[3])
                }
            }
            else # not mixture, pure population with (relatively) fixed states
            {
                variable <-
                    list(mu = mu,
                         size = size,
                         prob = prob)
                fixed <-
                    list(states = states.scale)
                mean.states <- function(pars) pars$mu * pars$states
                demission <- function(k, x, pars, ...)
                {
                    dnbinom(x, mu = pars$mu * pars$states[k],
                            size = pars$size,
                            log = TRUE)
                }
                updated.var <-
                    function(xlist, Akl, Bk, Nkb, pars,
                             control = list(), ...)
                {
                    xx <- unlist(xlist)
                    Q <- function(optpars, pars)
                    {
                        ## print(optpars)
                        mu <- optpars[1]
                        sigma <- optpars[2]
                        nbmean <- mu * pars$states
                        ans <- 
                            length(xx) * (sigma * log(sigma) - lgamma(sigma)) +
                            sum(lgamma(xx + sigma)) +
                            sum(Nkb * (log(nbmean) -
                                       log(sigma + nbmean)) -
                                sigma * Bk * log(sigma + nbmean))
                        ##print(ans)
                        -ans ## to be minimized
                    }
                    res <-
                        nlminb(c(pars$mu, pars$size),
                               Q,
                               control = control,
                               lower = 1e-5, upper = Inf,
                               pars = pars)
                    if (res$convergence != 0)
                        warning(paste("nlminb did not converge:",
                                      res$message))
                    list(prob = Akl / rowSums(Akl),
                         mu = res$par[1],
                         size = res$par[2])
                }
            }
        }
    }
    list(demission = demission,
         pars = c(variable, fixed),
         updated.var = updated.var,
         mean.states = mean.states)
}






## special case: neg binom conditioned on another dataset

nb.cond.family <-
    function(prob, mu, lambda, sigma)
{
    len <- length(lambda)
    if (missing(prob)) prob <- matrix(1 / len, len, len)
    variable <-
        list(states = lambda,
             sigma = sigma,
             prob = prob)
    fixed <- list(mu = mu)
    mean.states <- function(pars) pars$states ## not really mean
    demission <- function(k, x, y, pars, ...)
    {
        dnbinom(x,
                mu = pars$states[k] * pars$mu,
                size = y + pars$sigma,
                log = TRUE)
    }
    updated.var <-
        function(xlist, ylist, Akl, Bk, Nkb, Mkb, pars,
                 control = list(), ...)
        {
            str(Nkb)
            str(Mkb)
            xx <- unlist(xlist)
            yy <- unlist(ylist)

            Q <- function(sigma, pars)
            {
## FIXME: remove later
stopifnot(length(Nkb) == length(Mkb), length(Nkb) == length(Bk))
                sbpm <- (sigma * Bk + Mkb)
                ans <-
                    ## FIXME: first line reduces to simpler sums of sums
                    sum(lgamma(xx + yy + sigma)) - sum(lgamma(yy + sigma)) -
                        sum( sbpm * log1p(Nkb / sbpm) + Nkb * log1p(sbpm / Nkb) )
                -ans ## to be minimized
            }
            res <-
                nlminb(pars$sigma,
                       Q,
                       control = control,
                       lower = 1e-5, upper = Inf,
                       pars = pars)
            if (res$convergence != 0)
                warning(paste("nlminb did not converge:",
                              res$message))
            list(states = (1 + res$par / pars$mu) * Nkb / (res$par * Bk + Mkb),
                 mu = pars$mu,
                 sigma = res$par,
                 prob = Akl / rowSums(Akl))
        }
    list(demission = demission,
         pars = c(variable, fixed),
         updated.var = updated.var,
         mean.states = mean.states)
}


