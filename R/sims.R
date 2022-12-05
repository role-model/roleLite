#' @title Unified Neutral Theory of Biodiversity Simulation
#'
#' @description Simulate a local community under the Unified Neutral Theory of
#' Biodiversity.
#'
#' @param J
#' @param Sm
#' @param Jm
#' @param nu
#' @param m
#' @param niter
#' @param niterTimestep
#' @param initLocal
#' @param initMeta
#'
#' @return a matrix with rows for timesteps and columns for each individual in
#' the local community. The cells of the matrix are the species IDs. NOTE: it's
#' recommended that you don't work directly with this output but rather pass it
#' to \code{getAbund}
#'
#' @export


untbSim <- function(J, Sm, Jm, nu, m,
                    niter, niterTimestep = floor(niter / 10),
                    initLocal = NULL, initMeta = NULL) {
    # neutral death function
    dfun <- function(x, ...) rep(1, length(x))

    .simWorkHorse(J, Sm, Jm, nu, m, dfun,
                  niter, niterTimestep, initLocal, initMeta)
}




#' @title Multi-species Competitive Coexistence Theory Simulation
#'
#' @description Simulate a local community under Competitive Coexistence Theory.
#'
#' @param J number of individuals in the local community
#' @param Sm number of species in the metacommunity
#' @param Jm number of individuals in the metacommunity
#' @param nu speciation probability (must be between 0 and 1)
#' @param m immigration probability (must be between 0 and 1)
#' @param alpha competition matrix
#' @param niter number of iterations to run
#' @param niterTimestep number of iterations to save output on
#' @param initLocal initial state of the local community; if \code{NULL} (the
#' default) it will be randomly drawn from metacommunity
#' @param initMeta initial state of the meta community; if \code{NULL} (the
#' default) it will be simulated as a log series parameterized by \code{Sm} and
#' \code{Jm}
#'
#' @return a matrix with rows for timesteps and columns for each individual in
#' the local community. The cells of the matrix are the species IDs. NOTE: it's
#' recommended that you don't work directly with this output but rather pass it
#' to \code{getAbund}
#'
#' @export

compSim <- function(J, Sm, Jm, nu, m, alpha,
                    niter, niterTimestep = floor(niter / 10),
                    initLocal = NULL, initMeta = NULL) {
    # competitive death function
    dfun <- function(x, a) {
        nn <- tabulate(x)
        ialive <- which(nn > 0)
        nn <- nn[ialive]

        # dd <- ((a[ialive, ialive] %*% nn) * nn)[, 1]
        dd <- (a[ialive, ialive] %*% nn)[, 1]

        dout <- numeric(max(x))
        dout[ialive] <- dd

        return(dout[x])
    }

    .simWorkHorse(J, Sm, Jm, nu, m, dfun,
                  niter, niterTimestep, initLocal, initMeta,
                  a = alpha)
}




# note: `death` should be a function that calculates death probs
# `...` are extra params passed to `death`
.simWorkHorse <- function(J, Sm, Jm, nu, m, death,
                          niter, niterTimestep,
                          initLocal, initMeta, ...) {
    # meta comm SAD
    if(is.null(initMeta)) {
        metaSAD <- .lseriesFromSN(Sm, Jm)
    } else {
        metaSAD <- initMeta / sum(initMeta)
        if(missing(Sm)) Sm <- length(metaSAD)
    }

    # set up initial state for local comm
    if(is.null(initLocal)) {
        n <- sample(Sm, size = J, replace = TRUE, prob = metaSAD)
    } else {
        n <- initLocal
    }

    # matrix to hold sim results
    nmat <- matrix(0,
                   nrow = floor(niter / niterTimestep) + 1,
                   ncol = J)

    # record initial state
    j <- 1
    nmat[j, ] <- n

    # keep track of max num of spp
    Smax <- Sm

    # run the sim
    for(i in 2:(niter + 1)) {
        # death
        dprob <- death(n, ...)
        idead <- sample(J, 1, prob = dprob)

        # birth or immigration
        r <- runif(1)
        if(r < m) { # immigration
            inew <- sample(Sm, 1, prob = metaSAD)
        } else { # local birth
            inew <- sample(n, 1)
        }

        n[idead] <- inew

        # speciation?
        s <- runif(1)
        if(s < nu) { # yes speciation
            Smax <- Smax + 1
            n[idead] <- Smax
        }

        # record output
        if(i %% niterTimestep == 0) {
            j <- j + 1
            nmat[j, ] <- n
        }
    }

    attr(nmat, 'niterTimestep') <- niterTimestep

    return(nmat)
}


.lseriesFromSN <- function(S, N) {
    # solve for alpha paramter
    asol <- uniroot(interval = c(.Machine$double.eps^0.25,
                                 .Machine$integer.max),
                    f = function(a) {
                        a * log(1 + N / a) - S
                    })

    # calculate p parameter and beta (as used by pika)
    p <- 1 - exp(-S / asol$root)
    beta <- -log(p)

    # calculate idealized SAD from parameter
    thisSAD <- pika::sad(model = 'lseries', par = beta)
    thisSAD <- pika::sad2Rank(thisSAD, S = S)

    # return relative abundances
    return(thisSAD / sum(thisSAD))
}

