
computeIntervals <-
    function(xref, xobs,
             target.obs = 15,
             prop = length(xobs) / length(xref),
             target.ref = NULL,  # overrides target.obs
             counts = TRUE)
{
    if (is.null(target.ref)) target.ref <- target.obs / prop
    nbins <- ceiling(length(xref) * target.ref)
    bins <- quantile(xref, prob = ppoints(nbins, a = 1), names = FALSE)
    attr(bins, "counts") <-
        unname(table(table(cut(x = xobs, breaks = bins, labels = NULL))))
    bins
}

