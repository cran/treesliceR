% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CpE.R
\name{CpE}
\alias{CpE}
\title{Calculates the rate of accumulation of phylogenetic endemism (CpE) over time slices}
\usage{
CpE(tree, n, mat, criterion = "my", pEO = 5, ncor = 0)
}
\arguments{
\item{tree}{phylo. An ultrametric phylogenetic tree in the "phylo" format.}

\item{n}{numeric. A numeric value indicating the number of temporal slices (method = 1) or the time interval in million years (or phylogenetic diversity) among the tree slices (method = 2). Default is 1.}

\item{mat}{matrix. A presence/absence matrix containing all studied species and sites.}

\item{criterion}{character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".}

\item{pEO}{numeric. A value indicating the numeric proportion to define the temporal origin at which the phylogenetic endemism (PE) started to accumulate in a given assemblage. Default is 5\%.}

\item{ncor}{numeric. A value indicating the number of cores the user wants to parallelize. Default is 0.}
}
\value{
The function returns a data frame containing the assemblages' rates of cumulative phylogenetic endemism (CpE), their total phylogenetic endemism (PE), and their PE origin (pEO).
}
\description{
This function estimates the rates of accumulation of phylogenetic endemism (CpE) over time for inputted assemblages.
}
\details{
\bold{Parallelization}

Users are advised to check the number of available cores within their machines before running parallel programming.
}
\examples{
# Generate a random tree
tree <- ape::rcoal(20)

# Create a presence-absence matrix
mat <- matrix(sample(c(1,0), 20*10, replace = TRUE), ncol = 20, nrow = 10)
colnames(mat) <- tree$tip.label

# Calculate the CpE for 100 tree slices
CpE(tree, n = 100, mat = mat)

}
\seealso{
Other cumulative phylogenetic index analysis: \code{\link[=CpD]{CpD()}}, \code{\link[=CpB]{CpB()}}, \code{\link[=CpB_RW]{CpB_RW()}}
}
\author{
Matheus Lima de Araujo \href{mailto:matheusaraujolima@live.com}{matheusaraujolima@live.com}
}
