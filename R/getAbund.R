#' @title Get abundances
#'
#' @description Extract a data.frame of species abundances across timesteps from
#' the output of simulation functions
#'
#' @param x the output from either \code{untbSim} or \code{compSim}
#' @param tstep integer giving the specific timestep(s) requested; if \code{NULL}
#' (the default) then all timesteps are used
#'
#' @return a data.frame with three columns: \code{tstep} gives the timestep,
#' \code{spID} gives the species ID, \code{abund} gives the abundance of that
#' species in that timestep. Note: this is a tidy-style data.frame
#'
#'
#' @export

getAbund <- function(x, tstep = NULL) {
    if(is.null(tstep)) tstep <- 1:nrow(x)

    if(length(tstep) > 1) {
        o <- lapply(tstep, function(thist) {
            getAbund(x, tstep = thist)
        })

        return(do.call(rbind, o))
    } else {
        y <- table(x[tstep, ])

        # browser()
        out <- data.frame(tstep = attr(x, 'niterTimestep') * (tstep - 1),
                          spID = as.numeric(names(y)),
                          abund = as.numeric(y))

        return(out)
    }
}
