

## At one point, we considered re-estimating the states as part of the
## Baum-Welch iteration.  Eventually, we decided not to use it, since
## the converged states are not always well-behaved (and we haven't
## formally investigated theoretical implications either).

## However, a partial re-estimation may make sense.  Here's the idea:
## the states will remain unchanged, but they need not necessarily
## define the emission probabilities fully.  Rather, the emission
## probabilities are defined by some (usually very few) further
## parameters in addition to the states.  These parameters will be
## re-estimated.

## We are interested in two specific examples: (1) X ~ Poisson(state *
## lambda), where lambda is the parameter to be estimated, but the
## states remain fixed (i.e., relative emission expectations are
## fixed). (2) X ~ neg.binom(mu * state, sigma), where mu and sigma
## are to be re-estimated.

## Of these, the Poisson re-estimation has a closed form solution, but
## the negative binomial (probably) doesn't.  We may want to
## experiment with other choices as well.  So, we use a somewhat
## general `family'-like (as in glm's) idea.



### forward.hmm, backward.hmm, baumwelch.hmm, etc may benefit from a C
### implementation.  The current R implementations use vectorization
### fairly efficiently in some cases (not forward.hmm though), so the
### extent of speed benefit is uncertain.  However, the R
### implementation is probably not very readable.


forward.hmm <- function(x, y = NULL, family)
    ## x: observed sequence
    ## family: HMM `family'
{
    prob <- family$pars$prob
    len <- length(family$pars$states)

    ## Used to be:
    ## ans <- outer(seq(length = len), x, family$demission, pars = family$pars)

    g <- expand.grid(k = seq(length = len), x = x)
    if (!is.null(y)) g$y <- y
    ans <-
        family$demission(k = g$k, x = g$x, y = g$y,
                         pars = family$pars)
    dim(ans) <- c(len, length(x))
    for(i in seq(along = x)[-1]) # takes care of length-1 x
    {

        ## adjusting to make sure that numbers being exponentiated are
        ## never too small.  The problem is, they should also not be
        ## too large, since then they could become Inf, which is much
        ## worse than becoming 0 --- one Inf can mess things up, while
        ## all have to be 0 before NaN's get produced (I used to take
        ## M <- min(.), but that caused this sort of problem).
        ## Choosing this optimally would need more work (I'm not sure
        ## how much), for now I'm going to take the easy way out and
        ## just subtract the maximum.  This comment applies to similar
        ## constructs below

        M <- max(ans[, i-1])
        ans[, i] <-
            ans[, i] + M + log(exp(ans[, i-1] - M) %*% prob)
    }
    ans
}




backward.hmm <- function(x, y = NULL, family)
{
    prob <- family$pars$prob
    len <- length(family$pars$states)

    ## we're never going to need e(x_1).  However, in the i-th step,
    ## we'll need both e(x_{i+1}) and b(i+1), so it's not as trivial
    ## to use the same storage for both (changing the i-th column from
    ## e to b in the i-th step).

    ## What we will do is: at the start of the step changing the i-th
    ## column, ans[, i] will contain e(x_{i+1}), at the end it will
    ## contain b(i). Thus the 1st column is initially e(x_2), and ends
    ## up with b(1). We'll never store e(x_1).

    ## To achieve this, we start with the matrix of probabilities
    ## skipping x_1, and append a column of 1's.

    g <- expand.grid(k = seq(length = len), x = c(x[-1], 0))
    if (!is.null(y)) g$y <- c(y[-1], 0)
    ans <-
        family$demission(k = g$k, x = g$x, y = g$y,
                         pars = family$pars)
    dim(ans) <- c(len, length(x))

    ## last column was just to prevent copying later, and is
    ## immediately set to 1 (b(L))
    ans[, length(x)] <- 0
    for (i in rev(seq(length = length(x) - 1)))
    {
        M <- max(ans[, i+1])
        ans[, i] <- M + log(prob %*% (exp( ans[, i] + ans[, i+1] - M ) ))
    }
    ans
}




posterior.hmm <-
    function(f, b, ...)
    ## f, b: forward and backward probabilities
    ## log: logical, whether f and b are in log scale
    ## (result will be in the same scale)
{
    M <- max(f[, ncol(f)])
    log.p <- M + log( sum ( exp( f[, ncol(f)] - M  ))  )
    f + b - log.p
}



baumwelch.hmm <-
    function(xlist, ylist = NULL,
             family,
             flist, blist,
             control = list())

    ## re-estimates transition probability matrix
    ## f, b in log scale. Note that prob is never in the
    ## log scale; not when input, not when output

{
    stopifnot (is.list(xlist), is.list(flist), is.list(blist))

    prob <- family$pars$prob
    nstates <- length(family$pars$states)
    ## print(states)

    ## Need several quantities. These are (all expectations over
    ## possible paths given x and current estimates):

    ## A_kl = expected number of transitions from states k -> l
    ## B_k = expected number of sites with state k
    ## N_k(b) = expected number of emissions b in state k
    ## (Actually, we need only \sum_b b N_k(b) for each k, and that's
    ## what we will henceforth mean by Nkb)
    ## N(b) = expected number of emissions b (not random)


    Akl <- matrix(0, nrow = nrow(prob), ncol = ncol(prob))
    Bk <- numeric(nstates)
    Nkb <- numeric(nstates)
    Mkb <- if (!is.null(ylist)) numeric(nstates) else NULL
    
    Akl[,] <- 0
    for (i in seq(along = xlist))
    {
        ## Note: prob should be unchanged for different i

        x <- xlist[[i]]
        y <- ylist[[i]] ## could be NULL
        f <- flist[[i]]
        b <- blist[[i]]

        lenx <- length(x) ## should equal ncol(f)
        M <- max(f[, lenx])
        log.p <- M + log( sum ( exp( f[, lenx] - M  ))  )  # p = P(x)
        logprob <- log(prob) ## prob remains unchanged (current parameter)

        ## Akl: small loop
        for (k in seq(length = nstates))
            for (l in seq(length = nstates))
            {
                tmp <- f[k, -lenx] + b[l, -1] + 
                    family$demission(k = l, x = x[-1], y = y[-1], pars = family$pars)
                M <- max(tmp)
                Akl[k,l] <-
                    Akl[k,l] +
                        exp(logprob[k, l] - log.p + 
                            M + log(sum(exp(tmp - M))))
## if (!is.finite(Akl[k,l])) browser()
                ## log(sum(f.e.b))
            }
        tmp <- f + b
        M <- apply(tmp, 1, max)
        Bk <- Bk +
            exp( M + log(rowSums(exp(tmp - M))) - log.p)
        ##       <-----  log(sum(f.b)) ------->


        ## Nkb: Want x * f_k * b_k / p
        ## Can't work on log scale because x can be 0

        M <- apply(f + b - log.p, 1, max)
        tmp <- exp(f + b - log.p - M) *  
            matrix(x, nrow = nstates, ncol = lenx, byrow = TRUE)
        Nkb <- Nkb + exp(M) * rowSums(tmp)

        ## similarly Mkb: y * f_k * b_k / p
        if (!is.null(y))
        {
            tmp <- exp(f + b - log.p - M) *  
                matrix(y, nrow = nstates, ncol = lenx, byrow = TRUE)
            Mkb <- Mkb + exp(M) * rowSums(tmp)
        }
    }
    family$updated.var(xlist = xlist,
                       ylist = ylist,
                       Akl = Akl,
                       Bk = Bk,
                       Nkb = Nkb,
                       Mkb = Mkb,
                       pars = family$pars)
}





## the rest of the code is for trying to do useful things with the
## tools defined above given one or more sequences







coverageHmm <-
    function(xlist,
             ...,
             ylist = NULL,
             binlist = NULL, ## locations of bins
             family,
             iterations = 0)
    ## ... : one or more coverage sequences
{
    nstates <- length(family$pars$states)

    if (!is.list(xlist)) xlist <- list(xlist, ...)
    ## remove length-1 sequences, as they mess things up
    xLengths <- sapply(xlist, length)
    xlist <- xlist[xLengths > 1]
    lenx <- length(xlist)
    flist <- vector(mode = "list", length = lenx)
    blist <- vector(mode = "list", length = lenx)
    postlist <- vector(mode = "list", length = lenx) ## necessary?
    names(flist) <- names(blist) <- names(postlist) <- names(xlist)
    for (i in seq(length = lenx))
    {
        flist[[i]] <-
            forward.hmm(x = xlist[[i]], y = ylist[[i]], family = family)
        blist[[i]] <-
            backward.hmm(x = xlist[[i]], y = ylist[[i]], family = family)
        postlist[[i]] <- 
            posterior.hmm(f = flist[[i]], b = blist[[i]])
    }
    ans <-
        list(xlist = xlist,
             ylist = ylist,
             binlist = binlist,
             flist = flist,
             blist = blist,
             postlist = postlist,
             family = family)
    class(ans) <- "coverageHmm"
    if (iterations > 0) update(ans, iterations = iterations)
    else ans
}





write.hmm <-
    function(x, bins = x$binlist, decoded = "mode",
             file = stop("file must be specified (\"\" for console output)"),
             seq.labels = names(x$xlist),
             sep = ",", quote = FALSE, row.names = FALSE,
             combine = FALSE,
             ...)
{
    nobs <- sapply(x$xlist, length)
    if (is.null(seq.labels)) seq.labels <- as.character(seq(along = nobs))
    tmp <- data.frame(chr = I(rep(seq.labels, nobs)))
    if (is.list(bins))
    {
        tmp$start <- round(unlist(lapply(bins, function(x) x[-length(x)])))
        tmp$end <- round(unlist(lapply(bins, function(x) x[-1])))
    }
    tmp$observed <- unlist(x$xlist)
    if (is.character(decoded)) decoded <- decode(x, method = decoded)
    if (is.list(decoded))
    {
        tmp$decoded <- round(unlist(decoded), digits = 4)
        ## tmp$prob <- unlist(lapply(decoded, attr, "prob")) # won't work this way
    }
    ## else what?  Why even have the if?
    if (combine)
    {
        if (is.null(bins)) stop("bins = NULL case not meaningful")
        tmp$observed <- NULL
        keep <- !logical(nrow(tmp))
        for (i in rev(2:nrow(tmp)))
        {
            if (tmp$chr[i] == tmp$chr[i-1] && tmp$decoded[i] == tmp$decoded[i-1])
            {
                keep[i] <- FALSE
                tmp$end[i-1] <- tmp$end[i]
            }
        }
        tmp <- subset(tmp, keep)
    }
    write.table(tmp,
                file = file,
                sep = sep,
                quote = quote,
                row.names = row.names,
                ...)
}




print.coverageHmm <- function(x, ...)
{
    cat(paste("HMM with", length(x$xlist), "sequences:"), fill = TRUE)
    str(x$xlist)
    cat("Model:", fill = TRUE)
    str(x$family)
    invisible(x)
}



summary.coverageHmm <-
    function(object, ...)
{
    ## works only for poisson and negbin families
    ans <-
        list(tp = round(1000 * object$family$pars$prob),
             sp = stationary.distribution(object$family$pars$prob),
             ms = mean.emissions(object),
             size = object$family$pars$size,
             ll = logLik(object))
    class(ans) <- "summary.coverageHmm"
    ans
}


print.summary.coverageHmm <-
    function(x, ...)
{
    print(x$ll)
    cat("\nEstimated transition probability matrix (thousandths):\n")
    print(x$tp)
    cat("\nEstimated means and stationary distribution:\n")
    print(cbind(mean = x$ms, prob = x$sp))
    cat(paste("\nMean of stationary distribution:",
              format(sum(x$ms * x$sp)),
              "\n"))
    if (!is.null(x$size))
        cat(paste("\nEstimated 'size' parameter", format(x$size), "\n"))
    invisible(x)
}



anova.coverageHmm <-
    function(object, object2, ...)
{
    ll1 <- logLik(object)
    ll2 <- logLik(object2)
    x <- abs(as.numeric(2 * (unclass(ll2) - unclass(ll1))))
    df <- abs(attr(ll2, "df") - attr(ll1, "df"))
    cat(paste("P-value:", round(1 - pchisq(x, df = df), 5)), fill = TRUE)
    invisible(pchisq(x, df = df))
}



## get mean emissions in the various states

mean.emissions <- function(x)
{
    ans <- x$family$mean.states(x$family$pars)
    if (any(diff(ans) < 0)) warning("mean states not ordered")
    ans
}


logLik.coverageHmm <-
    function(object, ...)
{
    ans <- 
        sum(sapply(object$flist,
                   function(f) {
                       lenx <- ncol(f)
                       M <- max(f[, lenx])
                       M + log( sum ( exp( f[, lenx] - M  ))  )  # log(P(x))
                   }))
    class(ans) <- "logLik"
    attr(ans, "nobs") <- # Needed for BIC. Number of normalized counts.
        sum(sapply(object$xlist, length))
    attr(ans, "df") <- # rows of $prob add up to 1
        length(unlist(object$family$pars)) - nrow(object$family$pars$prob)
    ans
}



update.coverageHmm <- 
    function(object,
             iterations = 10,
             stop.pc = 0, # stop if %change in likelihood drops below this
             verbose = interactive(),
             control = list())
{
    olik <- logLik(object)
    for (dummy in seq(length = iterations))
    {
        tmp <- # get updated parameters
            baumwelch.hmm(xlist = object$xlist,
                          ylist = object$ylist,
                          family = object$family,
                          flist = object$flist,
                          blist = object$blist,
                          control = control)
        object$family$pars[names(tmp)] <- tmp
        for (i in seq(along = object$xlist))
        {
            object$flist[[i]][] <-
                forward.hmm(x = object$xlist[[i]], y = object$ylist[[i]], family = object$family)
            object$blist[[i]][] <-
                backward.hmm(x = object$xlist[[i]], y = object$ylist[[i]], family = object$family)
        }
        nlik <- logLik(object)
        pc.change <- 100 * (nlik-olik) / abs(olik)
        if (verbose)
            cat(paste("Iteration: ",
                      dummy, "/", iterations,
                      " (", round(nlik, digits = 4), ")\t",
                      "Percent change: ",
                      round(pc.change, digits = 10),
                      sep = ""), fill = TRUE)
        olik <- nlik
        if (pc.change < stop.pc) break;
    }
    if (verbose) cat("\n")
    for (i in seq(along = object$xlist))
    {
        object$postlist[[i]][] <- 
            posterior.hmm(object$flist[[i]],
                          object$blist[[i]])
    }
    object
}



stationary.distribution <-
    function(prob) ## transition matrix
{
    ## getting stationary pi as lim pi0 P^n as n -> Inf

    nstates <- nrow(prob)
    stopifnot(all(dim(prob) == nstates))

    pi0 <- t(rep(1/nstates, nstates))
    pi1 <- pi0 %*% prob
    while (!identical(all.equal(pi0, pi1), TRUE))
    {
        pi0 <- pi1
        pi1 <- pi0 %*% prob
    }
    as.vector(pi1)
}





## this could (and should) be easily re-implemented in C

decode.viterbi <-
    function(x, y = NULL, family)
{
    prob <- family$pars$prob
    len <- length(family$pars$states)
    logprob <- log(prob)
    ## WAS: v <- outer(seq(length = len), x, FUN = family$demission, pars = family$pars)
    g <- expand.grid(k = seq(length = len), x = x)
    if (!is.null(y)) g$y <- y
    v <-
        family$demission(k = g$k, x = g$x, y = g$y,
                         pars = family$pars)
    dim(v) <- c(len, length(x))
    ptr <- v
    ptr[,] <- 0
    for(i in seq(along = x)[-1]) # takes care of length-1 x
    {
        tmp <- v[, i-1] + logprob
        tmax <- apply(tmp, 2, max)
        wmax <- apply(tmp, 2, which.max)
        v[, i] <- v[, i] + tmax
        ptr[, i] <- wmax
    }
    ans <- integer(length(x))
    ans[length(ans)] <- which.max(v[, length(x)])
    for (i in rev(seq(along = x))[-1])
    {
        ans[i] <- ptr[ans[i+1], i+1]
    }
    family$mean.states(family$pars)[ans]
}





decode.mode <-
    function(x, y = NULL, family)
{
    f <- forward.hmm(x = x, y = y, family = family)
    b <- backward.hmm(x = x, y = y, family = family)
    post <- posterior.hmm(f, b)
    ans <- family$mean.states(family$pars)[apply(post, 2, which.max)]
    attr(ans, "prob") <- exp(apply(post, 2, max))
    ans
}




decode.mean <-
    function(x, y = NULL, family)
{
    f <- forward.hmm(x = x, y = y, family = family)
    b <- backward.hmm(x = x, y = y, family = family)
    logpost <- posterior.hmm(f, b)
    post <- exp(logpost)
    states <- family$mean.states(family$pars)
    mean <- colSums(post * states)
    attr(mean, "sd") <- sqrt(colSums(post * states^2) - mean^2)
    attr(mean, "entropy") <- colSums(post * logpost, na.rm = TRUE)
    mean
}





decode <-
    function(hmm,
             xlist = hmm$xlist,
             method = c("viterbi", "mode", "mean"))
{
    stopifnot(class(hmm) == "coverageHmm")
    method <- match.arg(method)
    decode.fun <-
        switch(method,
               viterbi = function(x, ...) decode.viterbi(x, ...),
               mode = function(x, ...) decode.mode(x, ...),
               mean = function(x, ...) decode.mean(x, ...))
    lapply(xlist, decode.fun, 
           family = hmm$family)
}


