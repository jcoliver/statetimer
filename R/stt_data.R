#' Data frame of character state through time
#'
#' @param tree a phylogenetic tree object.
#' @param anc_states matrix of ancestral state probabilities, generally the
#' \code{states} element of a \code{corHMM} object.
#' @param tip_states matrix of states for tip taxa, can be the \code{tip.states}
#' element of a \code{corHMM} object.
#' @param state_index integer the index representing the state of interest in
#' the \code{corHMM$tip.states} and \code{corHMM$states} objects. Default is
#' 1.
#' @param contemporary_age numeric cutoff for determining age of nodes to
#' consider as contemporary. Defaults to 1.0e-10.
#'
#' @return data.frame with \code{node_age} and \code{lineage_count} columns for
#' plotting the character state of interest through time.
#'
#' @importFrom ape node.depth.edgelength
#' @import dplyr
#'
#' @export
stt_data <- function(tree, anc_states, tip_states, state_index = 1, contemporary_age = 1.0e-10) {
  # Want a table with:
  # Node    node_age    marg_state_prob
  # So we start by identifying which nodes are terminals and which are internals
  internals <- unique(tree$edge[order(tree$edge[, 1]), 1])
  terminals <- base::setdiff(x = tree$edge[, 2], y = internals)
  terminals <- terminals[order(terminals)]

  # Data frame we'll use for states
  # Calculate node depth, but convert to node age (ape::node.depth.edgelength
  # calculate depth starting from the *root*, so root depth = 0).
  age_state_df <- data.frame(node = c(terminals, internals),
                             node_age = max(ape::node.depth.edgelength(phy = tree)) - ape::node.depth.edgelength(phy = tree))

  # Want to add the marginal probabilities of the state of interest
  age_state_df$marg_state_prob <- c(tip_states[, state_index],
                                    anc_states[, state_index])

  # Reverse order of node age, we will be doing calculations starting with the
  # youngest (i.e. contemporary) nodes first
  age_state_df <- age_state_df[rev(order(age_state_df$node_age)), ]

  # Column to store count of lineages with state of interest
  age_state_df$lineage_count <- NA

  for (i in (nrow(age_state_df)):1) {
    if (i == nrow(age_state_df)) {
      # Looking at last (i.e. youngest) node, which shouldn't have any
      # descendents; initialize the sum of probabilities with marginal
      # probability at that (terminal) node
      age_state_df$lineage_count[i] <- age_state_df$marg_state_prob[i]
    } else {
      # Not the youngest node, so we already have a running tally in the
      # next row (i + 1) of the data frame; update accordingly
      node_i <- age_state_df$node[i]
      sum_marg <- 0
      # Find all descendants of that node (if any)
      descendents_i <- tree$edge[tree$edge[,1] == node_i, 2]
      if (length(descendents_i) > 0) {
        # Sum the marginal probabilities for those descendants
        sum_marg <- sum(age_state_df$marg_state_prob[age_state_df$node %in% descendents_i])
      }

      # Find marginal probability of node i
      node_i_marg <- age_state_df$marg_state_prob[age_state_df$node == node_i]

      # Calculate amount to change our tally of descendants
      change <- -sum_marg + node_i_marg
      age_state_df$lineage_count[i] <- age_state_df$lineage_count[i + 1] + change
    }
  }

  # We need a single entry for the contemporary taxa, so we need to (1) decide
  # how to identify "contemporary", (2) identify those contemporary rows, (3)
  # find the max value from those and use that for the entry, (4) drop all
  # contemporary rows, and (5) add a single entry to represent all contemporary
  # nodes

  # TODO: Combine steps 2 & 3 with pipes
  # (2) identify those contemporary rows
  contemporary_df <- age_state_df %>%
    filter(node_age <= contemporary_age)

  # (3) find the max value from those
  contemporary_point <- max(contemporary_df$lineage_count)

  # (4) drop all contemporary rows
  plot_age_states <- age_state_df %>%
    filter(node_age > contemporary_age)

  # (5) add a single entry to represent all contemporary nodes
  plot_age_states[nrow(plot_age_states) + 1, ] <- NA
  plot_age_states$lineage_count[nrow(plot_age_states)] <- contemporary_point
  plot_age_states$node_age[nrow(plot_age_states)] <- contemporary_age

  return(plot_age_states)
}
