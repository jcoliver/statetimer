% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joint_probs.R
\name{joint_probs}
\alias{joint_probs}
\title{Calculate all joint probabilities from a pair of ancestral state
reconstructions}
\usage{
joint_probs(x, y, x_tips = NULL, y_tips = NULL, include_tips = TRUE)
}
\arguments{
\item{x, y}{matrices of marginal probabilities for ancestral states;
generally the \code{states} element from a \code{rayDISC} object.}

\item{x_tips, y_tips}{matrices of tip states; generally the \code{tip.states}
element of a \code{rayDISC} object. If included, output will
include the joint probabilities for terminal nodes}

\item{include_tips}{logical indicating whether or not joint probabilities
for terminal nodes should be included in output. \strong{DEPRECATED}}
}
\value{
a list with a \code{states} matrix of joint probabilities for
ancestral state combinations at internal nodes and if \code{include_tips =
TRUE}, a \code{tip_states} matrix of joint probabilities of the state
combination at the terminal nodes.
}
\description{
Calculate all joint probabilities from a pair of ancestral state
reconstructions
}
