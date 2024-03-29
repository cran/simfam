% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{recomb_map_hg}
\alias{recomb_map_hg}
\alias{recomb_map_hg38}
\alias{recomb_map_hg37}
\title{Simplified recombination maps for human genomes}
\format{
A list with 22 elements (autosomes, not named), each a tibble with two columns defining the recombination map at that chromosome:
\itemize{
\item \code{pos}: position in base pairs
\item \code{posg}: position in centiMorgans (cM)
}

An object of class \code{list} of length 22.
}
\source{
Raw genetic maps downloaded from this location prior to above processing:
\url{https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/}

Chromosome lengths from:
\url{https://www.ncbi.nlm.nih.gov/grc/human/data}
}
\usage{
recomb_map_hg38

recomb_map_hg37
}
\description{
Human genetic recombination maps for builds 38 (GRCh38/hg38) and 37 (GRCh37/hg19, below suffixed as hg37 for simplicity although technically incorrect).
Processed each first with \code{\link[=recomb_map_fix_ends_chr]{recomb_map_fix_ends_chr()}} to shift and extrapolate to sequence ends, then simplified with \code{\link[=recomb_map_simplify_chr]{recomb_map_simplify_chr()}} to remove all values that can be extrapolated with an error of up to \code{tol = 0.1}, in order to reduce their sizes and interpolation runtime.
Defaults were used, which resulted in extrapolated recombination rates close to and centered around the average of 1e-6 cM/base).
Autosomes only.
}
\keyword{datasets}
