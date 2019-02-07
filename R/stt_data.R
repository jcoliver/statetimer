#' Create a data frame object for character state through time plot
#'
#' @param tree a phylogenetic tree object
#' @param anc.recon list with marginal ancestral state probabilities including
#' \code{states} and \code{tip.states} matrices. Can be a \code{rayDISC} object.
#' @param state.index integer the index representing the state of interest in
#' the \code{rayDISC$tip.states} and \code{rayDISC$states} objects. Default is
#' 1.
#' @param contemporary.age numeric cutoff for determining age of nodes to
#' consider as contemporary. Defaults to 1.0e-10.
#' @return data.frame with \code{node.age} and \code{lineage.count} columns for
#' plotting the character state of interest through time
stt_data <- function(tree, anc.recon, state.index = 1, contemporary.age = 1.0e-10) {
  # Want a table with:
  # Node    Node.age    marg.state.prob
  # So we start by identifying which nodes are terminals and which are internals
  internals <- unique(tree$edge[order(tree$edge[, 1]), 1])
  terminals <- base::setdiff(x = tree$edge[, 2], y = internals)
  terminals <- terminals[order(terminals)]

  # Data frame we'll use for states
  # Calculate node depth, but convert to node age (ape::node.depth.edgelength
  # calculate depth starting from the *root*, so root depth = 0).
  age.state.df <- data.frame(node = c(terminals, internals),
                             node.age = max(node.depth.edgelength(phy = tree)) - node.depth.edgelength(phy = tree))

  # Want to add the marginal probabilities of the state of interest
  age.state.df$marg.state.prob <- c(anc.recon$tip.states[, state.index],
                                    anc.recon$states[, state.index])

  # Reverse order of node age, we will be doing calculations starting with the
  # youngest (i.e. contemporary) nodes first
  age.state.df <- age.state.df[rev(order(age.state.df$node.age)), ]

  # Column to store count of lineages with state of interest
  age.state.df$lineage.count <- NA

  for (i in (nrow(age.state.df)):1) {
    if (i == nrow(age.state.df)) {
      # Looking at last (i.e. youngest) node, which shouldn't have any
      # descendents; initialize the sum of probabilities with marginal
      # probability at that (terminal) node
      age.state.df$lineage.count[i] <- age.state.df$marg.state.prob[i]
    } else {
      # Not the youngest node, so we already have a running tally in the
      # next row (i + 1) of the data frame; update accordingly
      node.i <- age.state.df$node[i]
      sum.marg <- 0
      # Find all descendants of that node (if any)
      descendents.i <- tree$edge[tree$edge[,1] == node.i, 2]
      if (length(descendents.i) > 0) {
        # Sum the marginal probabilities for those descendants
        sum.marg <- sum(age.state.df$marg.state.prob[age.state.df$node %in% descendents.i])
      }

      # Find marginal probability of node i
      node.i.marg <- age.state.df$marg.state.prob[age.state.df$node == node.i]

      # Calculate amount to change our tally of descendants
      change <- -sum.marg + node.i.marg
      age.state.df$lineage.count[i] <- age.state.df$lineage.count[i + 1] + change
    }
  }

  # We need a single entry for the contemporary taxa, so we need to (1) decide
  # how to identify "contemporary", (2) identify those contemporary rows, (3)
  # find the max value from those and use that for the entry, (4) drop all
  # contemporary rows, and (5) add a single entry to represent all contemporary
  # nodes

  # (2) identify those contemporary rows
  contemporary.df <- age.state.df[age.state.df$node.age <= contemporary.age, ]
  # (3) find the max value from those
  contemporary.point <- max(contemporary.df$lineage.count)
  # (4) drop all contemporary rows
  plot.age.states <- age.state.df[age.state.df$node.age > contemporary.age, ]
  # (5) add a single entry to represent all contemporary nodes
  plot.age.states[nrow(plot.age.states) + 1, ] <- NA
  plot.age.states$lineage.count[nrow(plot.age.states)] <- contemporary.point
  plot.age.states$node.age[nrow(plot.age.states)] <- contemporary.age

  return(plot.age.states)
}
