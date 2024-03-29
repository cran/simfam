% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/admix_fam.R
\name{admix_fam}
\alias{admix_fam}
\title{Calculate admixture matrix of a pedigree with known admixture of founders}
\usage{
admix_fam(admix, fam, missing_vals = c("", 0))
}
\arguments{
\item{admix}{The admixture proportions matrix of the founders (individuals along rows and ancestries along columns).
This matrix must have row names that identify each founder (matching codes in \code{fam$id}).
Individuals may be in a different order than \code{fam$id}.
Extra individuals in \code{admix} but absent in \code{fam$id} will be silently ignored.
All values should be non-negative and each row of \code{admix} should sum to one; for speed, this code does not check that \code{admix} is valid, just averages data as-is.}

\item{fam}{The pedigree data.frame, in plink FAM format.
Only columns \code{id}, \code{pat}, and \code{mat} are required.
\code{id} must be unique and non-missing.
Founders must be present, and their \code{pat} and \code{mat} values must be missing (see below).
Non-founders must have both their parents be non-missing.
Parents must appear earlier than their children in the table.}

\item{missing_vals}{The list of ID values treated as missing.
\code{NA} is always treated as missing.
By default, the empty string ('') and zero (0) are also treated as missing (remove values from here if this is a problem).}
}
\value{
The admixture proportions matrix of the entire \code{fam} table, based on the admixture of the founders.
These are expectations, calculated for each individual as the average ancestry proportion of the parents.
The rows of this admixture matrix correspond to \code{fam$id} in that order.
The columns (ancestries) are the same as in the input \code{admix}.
}
\description{
Calculates a full admixture proportions matrix (for all individuals in the provided pedigree FAM table) starting from the admixture proportions of the founders as provided.
}
\examples{
# The smallest pedigree, two parents and a child.
# A minimal fam table with the three required columns.
# Note "mother" and "father" have missing parent IDs, while "child" does not
library(tibble)
fam <- tibble(
  id = c('father', 'mother', 'child'),
  pat = c(NA, NA, 'father'),
  mat = c(NA, NA, 'mother')
)

# admixture proportions of the parents
admix <- rbind( c(0.3, 0.3, 0.4), c(0.5, 0.25, 0.25) )
# Name the parents with same codes as in `fam`
# (order can be different)
rownames( admix ) <- c('mother', 'father')
# name ancestries too
colnames( admix ) <- c('African', 'European', 'Asian')

# Calculate the full admixture proportions matrix
admix_all <- admix_fam( admix, fam )

# This is a 3x3 matrix with row names matching fam$id.
# The parent submatrix equals the input (reordered),
# but now there's admixture to the child too (averages of parents)
admix_all

}
\seealso{
Plink FAM format reference:
\url{https://www.cog-genomics.org/plink/1.9/formats#fam}
}
