#' Calculate all joint probabilities from a pair of ancestral state
#' reconstructions
#'
#' @param x,y marginal probabilities; generally the \code{states} element from
#' a \code{rayDISC} object.
#' @param include_tips logical indicating whether or not joint probabilities
#' for terminal nodes should be included in output.
#'
#' @returns a list with a \code{states} matrix of joint probabilities for
#' ancestral state combinations at internal nodes and if \code{include_tips =
#' TRUE}, a \code{tip_states} matrix of joint probabilities of the state
#' combination at the terminal nodes.
#'
#' @export
joint_probs <- function(x, y, include_tips = TRUE) {
  # Start by finding out how large the resultant matrix will need to be
  num_combinations <- ncol(x$states) * ncol(y$states)

  states <- matrix(data = NA,
                   nrow = nrow(x$states),
                   ncol = num_combinations)
  colnames(states)[1:ncol(states)] <- NA

  if (include_tips) {
    tip_states <- matrix(data = NA,
                         nrow = nrow(x$tip_states),
                         ncol = num_combinations)
    colnames(tip_states)[1:ncol(tip_states)] <- NA
  }

  state_col <- 1
  for (i_x in 1:ncol(x$states)) {
    for (i_y in 1:ncol(y$states)) {
      states[, state_col] <- x$states[, i_x] * y$states[, i_y]
      colnames(states)[state_col] <- paste0(colnames(x$states)[i_x], colnames(y$states)[i_y])

      if (include_tips) {
        tip_states[, state_col] <- x$tip_states[, i_x] * y$tip_states[, i_y]
        colnames(tip_states)[state_col] <- paste0(colnames(x$tip_states)[i_x], colnames(y$tip_states)[i_y])
      }

      state_col <- state_col + 1
    }
  }

  joint_probs <- list(states = states)

  if (include_tips){
    joint_probs$tip_states <- tip_states
  }

  return(joint_probs)
}
