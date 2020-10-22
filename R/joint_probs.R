#' Calculate all joint probabilities from a pair of ancestral state
#' reconstructions
#'
#' @param x,y matrices of marginal probabilities for ancestral states;
#' generally the \code{states} element from a \code{rayDISC} object.
#' @param x_tips,y_tips matrices of tip states; generally the \code{tip.states}
#' element of a \code{rayDISC} object. If included, output will
#' include the joint probabilities for terminal nodes
#'
#' @param include_tips logical indicating whether or not joint probabilities
#' for terminal nodes should be included in output. \strong{DEPRECATED}
#'
#' @returns a list with a \code{states} matrix of joint probabilities for
#' ancestral state combinations at internal nodes and if \code{include_tips =
#' TRUE}, a \code{tip_states} matrix of joint probabilities of the state
#' combination at the terminal nodes.
#'
#' @export
joint_probs <- function(x, y, x_tips = NULL, y_tips = NULL, include_tips = TRUE) {
  # Start by finding out how large the resultant matrix will need to be
  num_combinations <- ncol(x) * ncol(y)

  states <- matrix(data = NA,
                   nrow = nrow(x),
                   ncol = num_combinations)
  colnames(states)[1:ncol(states)] <- NA

  include_tips <- FALSE
  if (!is.null(x_tips) & !is.null(y_tips)) {
    tip_states <- matrix(data = NA,
                         nrow = nrow(x_tips),
                         ncol = num_combinations)
    colnames(tip_states)[1:ncol(tip_states)] <- NA
    include_tips <- TRUE
  }

  state_col <- 1
  for (i_x in 1:ncol(x)) {
    for (i_y in 1:ncol(y)) {
      states[, state_col] <- x[, i_x] * y[, i_y]
      colnames(states)[state_col] <- paste0(colnames(x)[i_x], colnames(y)[i_y])

      if (include_tips) {
        tip_states[, state_col] <- x_tips[, i_x] * y_tips[, i_y]
        # Assumes same order of states in states and tip states matrices :-\
        colnames(tip_states)[state_col] <- paste0(colnames(x)[i_x], colnames(y)[i_y])
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
