#' Calculate all joint probabilities from a pair of ancestral state
#' reconstructions
#'
#' @param x, y marginal probabilities; generally the \code{states} element from
#' a \code{rayDISC} object
#' @param include.tips [optional] logical indicating whether or not joint
#' probabilities for terminal nodes should be included in output. Default is
#' \code{TRUE}
#'
#' @references a list with a \code{states} matrix of joint probabilities for
#' ancestral state combinations at internal nodes and if \code{include.tips =
#' TRUE}, a \code{tip.states} matrix of joint probabilities of the state
#' combination at the terminal nodes.
joint_probs <- function(x, y, include.tips = TRUE) {
  # Start by finding out how large the resultant matrix will need to be
  num.combinations <- ncol(x$states) * ncol(y$states)

  states <- matrix(data = NA,
                   nrow = nrow(x$states),
                   ncol = num.combinations)
  colnames(states)[1:ncol(states)] <- NA

  if (include.tips) {
    tip.states <- matrix(data = NA,
                         nrow = nrow(x$tip.states),
                         ncol = num.combinations)
    colnames(tip.states)[1:ncol(tip.states)] <- NA
  }

  state.col <- 1
  for (i.x in 1:ncol(x$states)) {
    for (i.y in 1:ncol(y$states)) {
      states[, state.col] <- x$states[, i.x] * y$states[, i.y]
      colnames(states)[state.col] <- paste0(colnames(x$states)[i.x], colnames(y$states)[i.y])

      if (include.tips) {
        tip.states[, state.col] <- x$tip.states[, i.x] * y$tip.states[, i.y]
        colnames(tip.states)[state.col] <- paste0(colnames(x$tip.states)[i.x], colnames(y$tip.states)[i.y])
      }

      state.col <- state.col + 1
    }
  }

  joint.probs <- list(states = states)

  if (include.tips){
    joint.probs$tip.states <- tip.states
  }

  return(joint.probs)
}
