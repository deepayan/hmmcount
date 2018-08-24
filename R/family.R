

## HMM `families': generic specifiction of a HMM model with
## re-estimation rules


## note: demission needs to have state index (k) as first argument for
## convenience when calling outer.

hmm.family <-
    function(model = c("poisson", "nbinom"),
             prob, mu, states.scale, size,
             states.free = FALSE)
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
        else
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
                        ##print(optpars)
                        mu <- optpars[1]
                        sigma <- optpars[2]
                        ans <- 
                            length(xx) * (sigma * log(sigma) - lgamma(sigma)) +
                                sum(lgamma(xx + sigma)) +
                                    sum(Nkb * (log(pars$states * mu) -
                                               log(sigma + pars$states * mu)) -
                                        sigma * Bk * log(sigma + pars$states * mu))
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


