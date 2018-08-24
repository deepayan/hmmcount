
library(hmmcount)

run1 <- rpois(50, lambda = c(70, 10, 30))
run2 <- rpois(50, lambda = c(70, 10, 30))

lambdav <- c(0.25, 0.5, 1, 2, 4) * 10

mu1 <- rep(sample(lambdav, 50, TRUE, c(1, 3, 10, 3, 1)), run1)
mu2 <- rep(sample(lambdav, 50, TRUE, c(1, 2, 5, 10, 2)), run2)

x1 <- rpois(length(mu1), lambda = mu1)
x2 <- rpois(length(mu2), lambda = mu2)


y1 <- rnbinom(length(mu1), mu = mu1, size = 5)
y2 <- rnbinom(length(mu2), mu = mu2, size = 5)



fm1 <-
    coverageHmm(x1, x2,
                family =
                hmm.family("poisson",
                           mu = 1, 
                           states.scale = c(5, 7.5, 10, 15, 20),
                           states.free = TRUE))

system.time(fm1.up <- update(fm1, 30))
xyplot(fm1.up, decode = "mode")
xyplot(fm1.up, decode = "mean", aux = list(mu1, mu2))
trellis.last.object(col = c("transparent", "grey", "black"),
                    lwd = c(0, 3, 1))
system.time(print(xyplot(fm1.up, decode = "viterbi", aux = list(mu1, mu2))))
fm1.up$family$pars

rootogram(fm1.up, aspect = "xy")
rootogram(fm1.up, aspect = "xy", marginal = FALSE)



fm2 <-
    coverageHmm(list(y1 = y1, y2 = y2),
                family =
                hmm.family("nbinom",
                           mu = 1, size = 30,
                           states.scale = c(5, 7.5, 10, 15, 20),
                           states.free = TRUE))

system.time(fm2.up <- update(fm2, 30))
xyplot(fm2.up, decode = "mode", strip = TRUE)
xyplot(fm2.up, decode = "mean", aux = list(mu1, mu2))
trellis.last.object(col = c("transparent", "grey", "black"),
                    lwd = c(0, 3, 1))
system.time(print(xyplot(fm2.up, decode = "viterbi", aux = list(mu1, mu2))))
fm2.up$family$pars

rootogram(fm2.up, aspect = "xy")


fm3 <-
    coverageHmm(y2,
                family =
                hmm.family("poisson",
                           mu = 1, 
                           states.scale = c(5, 10, 15, 20),
                           states.free = TRUE))

system.time(fm3.up <- update(fm3, 40))
xyplot(fm3.up, decode = "mode")
xyplot(fm3.up, decode = "mean", aux = list(mu2))

rootogram(fm3.up, aspect = "xy")
rootogram(fm3.up, aspect = "xy", marginal = FALSE)






fm4 <-
    coverageHmm(y2,
                family =
                hmm.family("nbinom",
                           mu = 1, size = 30,
                           states.scale = c(5, 10, 15, 20),
                           states.free = TRUE))

system.time(fm4.up <- update(fm4, 40))
xyplot(fm4.up, decode = "mode")
xyplot(fm4.up, decode = "mode", aux = list(mu2))

rootogram(fm4.up, aspect = "xy")
rootogram(fm4.up, aspect = "xy", marginal = FALSE)


fm5 <-
    coverageHmm(list(y1 = y1, y2 = y2),
                family =
                hmm.family("nbinom",
                           mu = 1, size = 30,
                           states.scale = c(5, 10, 15, 20),
                           states.free = TRUE))

system.time(fm5.up <- update(fm5, 25))
xyplot(fm5.up, decode = "mode", aux = list(mu1, mu2), strip = TRUE)
rootogram(fm5.up, aspect = "xy", strip = TRUE)
rootogram(fm5.up, aspect = "xy", marginal = FALSE, strip = TRUE)
