## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- cache = FALSE, include = FALSE------------------------------------------
## copied from examples from the "simmer" R package
## after: https://www.enchufa2.es/archives/suggests-and-vignettes.html
## by Iñaki Úcar
required <- c("bnpsd", "popkin", "kinship2", "RColorBrewer")

if (!all(sapply(required, requireNamespace, quietly = TRUE)))
  knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
library(simfam)

# data dimensions
n <- 15 # number of individuals per generation
G <- 3 # number of generations

# the result changes in every run unless we set the seed
set.seed(1)

# draw the random pedigree with those properties
data <- sim_pedigree( n, G )
# save these objects separately
fam <- data$fam
ids <- data$ids
kinship_local_G <- data$kinship_local

## -----------------------------------------------------------------------------
# the top of the table
head( fam )

# and the bottom
tail( fam )

## -----------------------------------------------------------------------------
kinship_local_G

## ---- fig.width = 7, fig.height = 4, fig.align = 'center'---------------------
# load the package
library(kinship2)

# Convert the `fam` table into a `pedigree` object required by the `kinship2` package.
obj <- pedigree(
    id = fam$id,
    dadid = fam$pat,
    momid = fam$mat,
    sex = fam$sex,
    affected = fam$pheno
)

# plot the pedigree
plot( obj, cex = 0.7 )

## ---- fig.width = 7, fig.height = 4, fig.align = 'center'---------------------
# replace original with pruned FAM table
# containing only ancestors of last generation
fam <- prune_fam( fam, ids[[ G ]] )
# inspect
fam
# and plot pedigree
obj <- pedigree(
    id = fam$id,
    dadid = fam$pat,
    momid = fam$mat,
    sex = fam$sex,
    affected = fam$pheno
)
plot( obj, cex = 0.7 )

## -----------------------------------------------------------------------------
# number of loci to draw
m <- 5000
# draw uniform ancestral allele frequencies
p_anc <- runif( m )
# draw independent (unstructured) genotypes for founders (includes founders without descendants)
X_1 <- matrix(
    rbinom( n * m, 2, p_anc ),
    nrow = m,
    ncol = n
)
# `simfam` requires names for the founders on `X_1`, to validate identities and align
colnames( X_1 ) <- ids[[ 1 ]]

## ---- fig.width = 7, fig.height = 3.4, fig.align = 'center'-------------------
library(popkin)
# estimate kinship
kinship_est_1 <- popkin( X_1 )
# true/expected kinship
kinship_1 <- diag( n ) / 2
# assign same IDs as `X_1` and `fam`; `simfam` generally requires these to validate/align
colnames( kinship_1 ) <- colnames( X_1 )
rownames( kinship_1 ) <- colnames( X_1 )
# plot them together for comparison
plot_popkin(
    list( kinship_1, kinship_est_1 ),
    titles = c('Expected', 'Estimated'),
    names = TRUE,
    leg_width = 0.2,
    mar = c(3, 2)
)

## -----------------------------------------------------------------------------
X <- geno_fam( X_1, fam )
# expected kinship of entire pruned pedigree, starting from true kinship of founders.
# Called "local" kinship because it ignores the ancestral population structure (will be reused)
# NOTE: agrees with `kinship2::kinship( fam )`, but only because founders are unstructured
kinship_local <- kinship_fam( kinship_1, fam )

## ---- fig.width = 7, fig.height = 3.4, fig.align = 'center'-------------------
# estimate kinship
kinship_local_est <- popkin( X )
# plot them together for comparison
plot_popkin(
    list( kinship_local, kinship_local_est ),
    titles = c('Expected', 'Estimated'),
    names = TRUE,
    names_cex = 0.8,
    leg_width = 0.2,
    mar = c(3, 2)
)

## ---- fig.width = 7, fig.height = 2.3, fig.align = 'center'-------------------
# indexes (really logicals) marking individuals in last generation,
# to subset data (which is aligned with `fam`)
indexes_G <- fam$id %in% ids[[ G ]]
# estimate kinship
kinship_local_est <- popkin( X )
# plot them together for comparison
plot_popkin(
    list(
        kinship_local[ indexes_G, indexes_G ],
        kinship_local_est[ indexes_G, indexes_G ],
        kinship_local_G
    ),
    titles = c('Expected', 'Estimated', 'Local'),
    names = TRUE,
    mar = c(3, 2)
)

## -----------------------------------------------------------------------------
library(bnpsd)
# additional admixture parameters:
# number of ancestries
K <- 3
# define population structure
# FST values for k=3 subpopulations
inbr_subpops <- c(0.2, 0.4, 0.6)
# admixture proportions from 1D geography (founders)
admix_proportions_1 <- admix_prop_1d_linear( n, K, sigma = 1 )
# add founder names to `admix_proportions_1` (will propagate to other vars of interest)
rownames( admix_proportions_1 ) <- ids[[ 1 ]]
# also name subpopulations simply as S1, S2, S3
colnames( admix_proportions_1 ) <- paste0( 'S', 1:K )

# calculate true kinship matrix of the admixed founders
# NOTE: overwrites `kinship` for unstructured founders
kinship_1 <- coanc_to_kinship( coanc_admix( admix_proportions_1, inbr_subpops ) )

# draw genotypes from admixture model
# NOTE: overwrites `X_1` for unstructured founders
X_1 <- draw_all_admix( admix_proportions_1, inbr_subpops, m )$X

## ---- fig.width = 7, fig.height = 3.4, fig.align = 'center'-------------------
# estimate kinship (NOTE: overwrites `kinship_est_1` for unstructured founders)
kinship_est_1 <- popkin( X_1 )
# plot them together for comparison
plot_popkin(
    list( kinship_1, kinship_est_1 ),
    titles = c('Expected', 'Estimated'),
    names = TRUE,
    leg_width = 0.2,
    mar = c(3, 2)
)

## -----------------------------------------------------------------------------
X <- geno_fam( X_1, fam )
# expected kinship of entire pruned pedigree, starting from true kinship of founders
# NOTE: DOESN'T agree with `kinship2::kinship( fam )` anymore because founders are structured!
kinship_total <- kinship_fam( kinship_1, fam )

## ---- fig.width = 7, fig.height = 3.4, fig.align = 'center'-------------------
# estimate kinship
kinship_total_est <- popkin( X )
# plot them together for comparison
plot_popkin(
    list( kinship_total, kinship_total_est ),
    titles = c('Expected', 'Estimated'),
    names = TRUE,
    names_cex = 0.8,
    leg_width = 0.2,
    mar = c(3, 2)
)

## ---- fig.width = 7, fig.height = 2.3, fig.align = 'center'-------------------
plot_popkin(
    list(
        kinship_total[ indexes_G, indexes_G ],
        kinship_total_est[ indexes_G, indexes_G ],
        kinship_local_G
    ),
    titles = c('Total Expected', 'Total Estimated', 'Local'),
    names = TRUE,
    mar = c(3, 2)
)

## ---- fig.width = 7, fig.height = 2.5, fig.align = 'center'-------------------
# admixture proportions of pruned family
admix_proportions <- admix_fam( admix_proportions_1, fam )

# visualize as bar plot

# for nice colors
library(RColorBrewer)
# colors for independent subpopulations
col_subpops <- brewer.pal( K, "Paired" )
# shrink default margins
par_orig <- par(mar = c(4, 4, 0.5, 0) + 0.2)
barplot(
    t( admix_proportions ),
    col = col_subpops,
    legend.text = TRUE,
    args.legend = list( cex = 0.7, bty = 'n' ),
    border = NA,
    space = 0,
    ylab = 'Admixture prop.',
    las = 2
)
mtext('Individuals', side = 1, line = 3)
par( par_orig ) # reset `par`

## -----------------------------------------------------------------------------
# this is the admixture-only relatedness
kinship_admix <- coanc_to_kinship( coanc_admix( admix_proportions, inbr_subpops ) )

## ---- fig.width = 7, fig.height = 6.8, fig.align = 'center'-------------------
# this is the approximation of the total kinship:
# NOTE: a bit complicated because diagonal has to be transformed.
# Approximation operates on inbreeding coefficients, not self-kinship,
# so `popkin::inbr_diag` converts self-kinship to inbreeding,
# `bnpsd::coanc_to_kinship` converts inbreeding back to self-kinship.
# Off-diagonal values are unmodified by both functions.
kinship_total_approx <- coanc_to_kinship(
    inbr_diag( kinship_local ) +
    inbr_diag( kinship_admix ) -
    inbr_diag( kinship_admix ) * inbr_diag( kinship_local )
)

# plot all model kinship matrices (no estimates from genotypes this time)
plot_popkin(
    list( kinship_admix, kinship_local, kinship_total, kinship_total_approx ),
    titles = c('Admixture only', 'Family only', 'Total', 'Total approx.'),
    names = TRUE,
    names_cex = 0.8,
    leg_width = 0.2,
    layout_rows = 2,
    mar = c(3, 2)
)

## -----------------------------------------------------------------------------
library(simfam)

# data dimensions
n <- 500 # number of individuals per generation
G <- 10 # number of generations
K <- 3 # number of ancestries
m <- 5000 # number of loci

# draw the random pedigree with those properties.
data <- sim_pedigree( n, G )
fam <- data$fam
ids <- data$ids
# use this local kinship calculation in plot
kinship_local_G <- data$kinship_local

# prune pedigree to speed up simulations/etc
fam <- prune_fam( fam, ids[[ G ]] )

### admixture model

# define population structure
# FST values for k=3 subpopulations
inbr_subpops <- c(0.2, 0.4, 0.6)
# admixture proportions from 1D geography (founders)
admix_proportions_1 <- admix_prop_1d_linear( n, K, sigma = 1 )
# add founder names to `admix_proportions_1` (will propagate to other vars of interest)
rownames( admix_proportions_1 ) <- ids[[ 1 ]]
# also name subpopulations simply as S1, S2, S3
colnames( admix_proportions_1 ) <- paste0( 'S', 1:K )

# draw genotypes from admixture model
# admixed founders
X_1 <- draw_all_admix( admix_proportions_1, inbr_subpops, m )$X
# draw genotypes, returning last generation only
X_G <- geno_last_gen( X_1, fam, ids )
# total kinship estimated from genotypes
kinship_total_G_est <- popkin( X_G )

# calculate true total kinship matrix
kinship_total_1 <- coanc_to_kinship( coanc_admix( admix_proportions_1, inbr_subpops ) )
kinship_total_G <- kinship_last_gen( kinship_total_1, fam, ids )

# kinship due to admixture only
admix_proportions_G <- admix_last_gen( admix_proportions_1, fam, ids )
kinship_admix_G <- coanc_to_kinship( coanc_admix( admix_proportions_G, inbr_subpops ) )

# total kinship approximation from components
kinship_total_approx_G <- coanc_to_kinship(
    inbr_diag( kinship_local_G ) +
    inbr_diag( kinship_admix_G ) -
    inbr_diag( kinship_admix_G ) * inbr_diag( kinship_local_G )
)

## ---- fig.width = 7, fig.height = 3.4, fig.align = 'center'-------------------
plot_popkin(
    list( kinship_total_G, kinship_total_G_est ),
    titles = c('Expected', 'Estimated'),
    leg_width = 0.2,
    mar = c(3, 2)
)

## ---- fig.width = 7, fig.height = 6.8, fig.align = 'center'-------------------
# plot all model kinship matrices (no estimates from genotypes this time)
plot_popkin(
    list( kinship_admix_G, kinship_local_G, kinship_total_G, kinship_total_approx_G ),
    titles = c('Admixture only', 'Family only', 'Total', 'Total approx.'),
    leg_width = 0.2,
    layout_rows = 2,
    mar = c(3, 2)
)

## ---- fig.width = 7, fig.height = 4, fig.align = 'center'---------------------
# the IDs in `fam` (the pruned FAM table) contain the individuals of interest
frac_per_gen <- vector('numeric', G)
for ( g in 1 : G ) {
    # this is the desired fraction
    frac_per_gen[g] <- sum( fam$id %in% ids[[ g ]] ) / n
}
# visualize as a barplot
barplot(
    frac_per_gen,
    names.arg = 1 : G,
    xlab = 'Generation',
    ylab = 'Frac. individuals with descendants'
)
# add a line marking the maximum fraction of individuals per generation
abline( h = 1, lty = 2, col = 'red' )

