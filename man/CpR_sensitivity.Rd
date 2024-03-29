% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CpR_sensitivity.R
\name{CpR_sensitivity}
\alias{CpR_sensitivity}
\title{Runs a sensitivity analysis for rates of accumulation of a given phylogenetic index}
\usage{
CpR_sensitivity(tree, vec, mat, adj, rate, samp, comp, method, criterion, ncor)
}
\arguments{
\item{tree}{phylo. An ultrametric phylogenetic tree in the "phylo" format.}

\item{vec}{numeric vector. A numeric vector containing a series of numbers of slices.}

\item{mat}{matrix. A presence/absence matrix containing all studied species and sites.}

\item{adj}{matrix. A square adjacency matrix containing the presence/absence information of all sites and their spatially adjacent ones.}

\item{rate}{character string. The desired cumulative phylogenetic rate to be assessed, which can be the phylogenetic diversity (CpD), phylogenetic endemism (CpE), phylogenetic B-diversity (CpB), or phylogenetic B-diversity range-weighted (CpB_RW). Default is NULL, but must be filled with "CpD", "CPE", "CpB_RW", or "CpB".}

\item{samp}{numeric. The number of assemblages, or sites, to be sampled to make the sensitivity analysis.}

\item{comp}{character string. The component of beta-diversity that the user wants to calculate the CpB. It can be either "sorensen", turnover" or "nestedness". This argument works only when "rate = CpB". Default is "sorensen".}

\item{method}{character string. The method for calculating the CpB-rate. It can be either "pairwise" or "multisite". This argument works only when the argument "rate" is set to run for "CpB" or "CpB_RW". Default is "multisite".}

\item{criterion}{character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".}

\item{ncor}{numeric. A value indicating the number of cores the user wants to parallelize. Default is 0.}
}
\value{
This function returns a data frame containing the sensitivity analysis for a given rate of accumulation of a phylogenetic index. This outputted data frame contains, for each row or assemblage, a column with the rate value assessed for each inputted number of slices.
}
\description{
This function allows the evaluation of the sensitivity of the estimated rates of accumulation of a given phylogenetic index (e.g., \code{\link[=CpD]{CpD()}}, \code{\link[=CpE]{CpE()}}, \code{\link[=CpB]{CpB()}}, \code{\link[=CpB_RW]{CpB_RW()}}) to the number of slices inputted by the user.
}
\details{
\bold{Parallelization}

Users are advised to check the number of available cores within their machines before running parallel programming.

\bold{Plotting}

For plotting the sensitivity analysis output users can use \code{\link[=CpR_sensitivity_plot]{CpR_sensitivity_plot()}}.
}
\examples{
# Generate a random tree
tree <- ape::rcoal(20)

# Create a presence-absence matrix
mat <- matrix(sample(c(1,0), 20*10, replace = TRUE), ncol = 20, nrow = 10)
colnames(mat) <- tree$tip.label

# Calculate the CpD for 100 tree slices
CpD(tree, n = 100, mat = mat)

# Create a vector of number of slices
vec <- c(25, 50, 75, 100, 125, 150)

# Calculate the sensitivity of the CpD
CpR_sensitivity(tree, vec, mat, rate = "CpD", samp = 5)

}
\references{
See the tutorial on how to use this function on our \href{https://araujomat.github.io/treesliceR/articles/Passeriformes-diversification.html}{website}.
}
\seealso{
Other cumulative phylogenetic index rate analysis: \code{\link[=CpD]{CpD()}}, \code{\link[=CpE]{CpE()}}, \code{\link[=CpB]{CpB()}}, \code{\link[=CpB_RW]{CpB_RW()}}
}
\author{
Matheus Lima de Araujo \href{mailto:matheusaraujolima@live.com}{matheusaraujolima@live.com}
}
