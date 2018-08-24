



xyplot.coverageHmm <-
    function(x,
             use.locs = !is.null(x$binlist),
             abscissa = NULL,
             abscissa.unit = if (use.locs) 1e6 else 1,
             decoded = "mode",
             aux = NULL,
             type = 'l',
             layout = c(1, length(nm)),
             strip = FALSE,
             ...)
{
    if (missing(abscissa))
        abscissa <-
            if (use.locs) 
            {
                if (is.null(x$binlist)) stop("Location information not available")
                lapply(x$binlist,
                       function(x) (x[-1] + x[-length(x)]) / (2 * abscissa.unit))
            }
            else 
                lapply(x$xlist, function(x) seq_along(x) / abscissa.unit)
    nm <- names(x$xlist)
    if (is.null(nm)) nm <- as.character(seq(along = x$xlist))
    dat <-
        data.frame(x = unlist(x$xlist, use.names = FALSE),
                   which = factor(rep(nm, sapply(x$xlist, length)),
                                  levels = nm),
                   location = unlist(abscissa))
    ## add decoded states
    if (is.character(decoded)) decoded <- decode(x, method = decoded)
    if (!is.null(decoded)) dat$decoded <- unlist(decoded)

    ## add auxiliary info, if any (e.g. true states in simulation)
    if (!is.null(aux)) dat$aux <- unlist(aux)

    ## generate formula
    form <- c("x", "decoded", "aux")
    form <- paste(form[form %in% names(dat)], collapse = " + ")
    form <- as.formula(paste(form, "~ location | which"))

    xyplot(form, dat,
           type = type,
           layout = layout,
           strip = strip,
           default.scales = list(axs = "i"),
           ...)
}




barchart.coverageHmm <-
    function(x,
             use.locs = !is.null(x$binlist),
             abscissa = NULL,
             abscissa.unit = if (use.locs) 1e6 else 1,
             gaps = NULL,   # location of seq gaps
             prepanel = function(x) list(xlim = range(x, finite = TRUE)),
             probs = do.call(cbind, lapply(x$postlist, exp)),
             layout = c(1, length(nm)),
             strip = FALSE,
             ylim = c(-0.1, 1.1), ylab = "",
             colors,
             ...,
             subscripts)
{
    if (missing(abscissa))
        abscissa <-
            if (use.locs) 
            {
                if (is.null(x$binlist)) stop("Location information not available")
                lapply(x$binlist,
                       function(x) (x[-1] + x[-length(x)]) / (2 * abscissa.unit))
            }
            else 
                lapply(x$xlist, function(x) seq_along(x) / abscissa.unit)
    nm <- names(x$xlist)
    if (is.null(nm)) nm <- as.character(seq(along = x$xlist))
    if (!is.null(gaps)) gaps <- gaps[nm]
    if (missing(colors))
    {
        colors <- as.list(seq(length = nrow(probs)))
        names(colors) <-
            rep(as.character(trellis.par.get("superpose.polygon")$col),
                length = nrow(probs))
    }
    dat <-
        data.frame(which = factor(rep(nm, sapply(x$xlist, length)),
                                  levels = nm),
                   location = unlist(abscissa))
    densityplot(~location | which, dat,
                prepanel = prepanel,
                ylim = ylim, ylab = ylab,
                default.scales = list(axs = "i"),
                layout = layout,
                strip = strip,
                probs = probs,
                gaps = gaps,
                panel = function(x, ..., probs, subscripts, gaps = NULL) {
                    cum.prob <- apply(probs[, subscripts], 2, cumsum)
                    ## add a 0 row to avoid special case-ing first row
                    cum.prob <- rbind(0, cum.prob)
                    wid <- min(diff(x)) ## width[subscripts]
                    for (col in names(colors))
                    {
                        for (row in colors[[col]])
                        {
                            lrect(xleft = c(0, x[-length(x)]),
                                  xright = x,
                                  y = 0.5 * (cum.prob[row+1,] + cum.prob[row,]),
                                  height = cum.prob[row+1,] - cum.prob[row,],
                                  col = col, border = "transparent")
                        }
                    }
                    if (!is.null(gaps))
                    {
                        g <- gaps[[packet.number()]]
                        require(grid, quietly = TRUE)
                        grid.rect(x = 0.5 * (g$start + g$end),
                                  width = g$end - g$start,
                                  default.units = "native",
                                  gp = gpar(fill = "lightgrey", col = "lightgrey"))
                    }
                },
                ...,
                subscripts = TRUE)
}



rootogram.coverageHmm <-
    function(x,
             strip = FALSE,
             ...,
             marginal = TRUE,
             separate = FALSE,
             subscripts)
{
    nm <- names(x$xlist)
    if (is.null(nm)) nm <- as.character(seq(along = x$xlist))
    dat <-
        data.frame(x = unlist(x$xlist),
                   which = factor(rep(nm, sapply(x$xlist, length)), levels = nm))
    family <- x$family
    if (marginal)
    {
        marginal.dfun <- function(x, family, stat.logprob) ## x vector
        {
            colSums(exp(outer(seq(along = stat.logprob), x,
                              family$demission,
                              pars = family$pars) + stat.logprob))
        }
        rootogram(~x | which, dat,
                  dfun = marginal.dfun,
                  family = family,
                  stat.logprob = log(stationary.distribution(family$pars$prob)),
                  strip = strip,
                  ...)
    }
    else if (!separate)
    {

        conditional.dfun <- function(x, family, post.logprob, subscripts)
        {
            states <- seq(length = nrow(post.logprob))
            ans <- numeric(length(x))
            for (i in subscripts)
            {
                ans <- ans +
                    colSums(exp(outer(states, x,
                                      family$demission,
                                      pars = family$pars) + post.logprob[,i]))
            }
            ans / length(subscripts)
        }
        rootogram(~x | which, dat,
                  dfun = conditional.dfun,
                  post.logprob = do.call(cbind, x$postlist),
                  family = family,
                  strip = strip,
                  ...,
                  subscripts = TRUE)
    }
    else NULL ## not written yet
}




## A generic copy-number plot with objects
## obs=data.frame(chrom, x, y)
## fit=data.frame(chrom, start, end, copy)
## log=TRUE/FALSE (if log, take log2 of both)

copyNumberPlot <- function(obs, fit, log = TRUE,
                           col.line = "black", lwd = 1, ...)
{
    if (isTRUE(log)) {
        log <- 2
        fit$copy <- log2(fit$copy)
    }
    xyplot(y ~ x | factor(chrom), data = obs,
           default.scales = list(y = list(log = log, equispaced.log = FALSE)),
           fit = split(fit, fit$chrom),
           panel = function(x, y, ..., fit) {
               panel.xyplot(x, y, ...)
               i <- packet.number()
               with(fit[[i]],
                    panel.segments(x0 = start, y0 = copy,
                                   x1 = end, y1 = copy,
                                   col = col.line, lwd = lwd))
           },
           ...)
}




