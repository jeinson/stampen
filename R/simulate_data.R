#' A function for generating fake data
#'
#' This function can simulate some haplotype combinations based on allele
#' frequencies sampled from beta distributions. The qtl frequencies are assumed
#' to represent the higher expressed allele. qtl frequencies are drawn from a
#' beta distribution with a = b = 1
#'
#' @param n_indvs The number of individuals to simulate
#' @param n_genes The number of genes to simulate
#' @param c_alpha shape parameter for the csnp allele frequency. Default = 4
#' @param c_beta shape parameter beta for the csnp allele frequency. Default = 3000
#' @param qtl_alpha shape parameter alpha for the qtl allele frequency. Default = 1
#' @param qtl_beta shape paramteter beta for the qtl allele frequency. Default = 1
#'
#' @export
#'
#' @examples
#' simulate_haplotype_counts(n_indvs = 100, n_genes = 5)
#'
#' @export
simulate_haplotype_counts <- function(
  n_indvs,
  n_genes,
  c_alpha = 0.656,
  c_beta = 437.47,
  qtl_alpha = 1,
  qtl_beta = 1)
{

  # Randomly draw the QTL allele frequencies for all genes simulated
  qtl_af <- rbeta(n_genes, qtl_alpha, qtl_beta)
  qtl_af <- .9 * qtl_af + .05 # Use a little y = mx+b to get the haplotypes to the right place

  # Simulate the allele frequency of coding variants corresponding to each gene
  csnp_af <- rbeta(n_genes, c_alpha, c_beta)

  # Compile the output
  haplotype_configs <- data.frame()
  for(i in 1:n_genes){

    # Generate sQTL haplotypes for all individuals for gene i
    ssnp_haps <-
      t(
        sapply(1:n_indvs, function(x){
          sample(c(0,1), 2, replace = T, prob = c(qtl_af[i], 1-qtl_af[i]))
        })
      )

    # Randomly draw which haplotype the coding SNP is on.
    csnp_haps <- t(sapply(1:n_indvs, function(x){
      sample(c(0,1), 2, replace = T, prob = c(1-csnp_af[i], csnp_af[i]))
    }))

    # Mark each type of haplotype combination
    haplotypes <- cbind(ssnp_haps, csnp_haps)

    abAB = 0
    ABab = 0
    abaB = 0
    aBab = 0

    AbaB = 0
    aBAb = 0
    AbAB = 0
    ABAb = 0

    for(j in 1:n_indvs){
      x <- haplotypes[j,]


      if(identical(x, c(0, 1, 0, 1))) abAB <- abAB + 1
      if(identical(x, c(1, 0, 1, 0))) ABab <- ABab + 1
      if(identical(x, c(0, 0, 0, 1))) abaB <- abaB + 1
      if(identical(x, c(0, 0, 1, 0))) aBab <- aBab + 1
      if(identical(x, c(1, 0, 0, 1))) AbaB <- AbaB + 1
      if(identical(x, c(0, 1, 1, 0))) aBAb <- aBAb + 1
      if(identical(x, c(1, 1, 0, 1))) AbAB <- AbAB + 1
      if(identical(x, c(1, 1, 1, 0))) ABAb <- ABAb + 1
    }

    out <- data.frame(gene = sprintf("gene%i", i), sqtl_af = qtl_af[i], AbaB, aBAb, AbAB, ABAb, abAB, ABab, abaB, aBab)
    haplotype_configs <- rbind(haplotype_configs, out)

  }
  return(haplotype_configs)
}

# Simulate the data that's included in this package
# set.seed(1234321)
# test_data <- simulate_haplotype_counts(1000, 5000)
# save(test_data, file = "data/test_data.rda")
#
# Save the haplotype configurations to the data directory
# beta_config_sqtl <-
#   c(abAB = 1,
#     ABab = 1,
#     abaB = 1,
#     aBab = 1,
#     AbaB = 0,
#     aBAb = 0,
#     AbAB = 0,
#     ABAb = 0
#   )
# save(beta_config_sqtl, file = "data/beta_config_sqtl.rda")
#
# beta_config_eqtl <-
#   c(abAB = 0,
#     ABab = 0,
#     abaB = 1,
#     aBab = 1,
#     AbaB = 1,
#     aBAb = 1,
#     AbAB = 0,
#     ABAb = 0
#   )
# save(beta_config_eqtl, file = "data/beta_config_eqtl.rda")

beta_plotter <- function(a,b, ...){
  x <- seq(0, 1, by = .0005)
  y <- dbeta(x, a, b)
  plot(x,y, type = "l", ...)
}
