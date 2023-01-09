#' A function for generating fake data
#'
#' This function can simulate some haplotype combinations based on allele
#' frequencies sampled from beta distributions. The qtl frequencies are assumed
#' to represent the higher expressed allele. qtl frequencies are drawn from a
#' beta distribution, defined by the parameters `qtl_alpha` and `qtl_beta`
#'
#' @param n_indvs The number of individuals to simulate
#' @param n_genes The number of genes to simulate
#' @param c_alpha shape parameter for the csnp allele frequency. Default = 4
#' @param c_beta shape parameter beta for the csnp allele frequency. Default = 3000
#' @param qtl_alpha shape parameter alpha for the qtl allele frequency. Default = 1
#' @param qtl_beta shape paramteter beta for the qtl allele frequency. Default = 1
#' @param maf_cutoff cutoff for the minimum qtl allele frequency when simulating data. This will force the beta distribution for the qtl allele frequency to these limits.
#'
#'
#' @examples
#' simulate_null_haplotypes(n_indvs = 100, n_genes = 5)
#'
#' @export
simulate_null_haplotypes <- function(
  n_indvs,
  n_genes,
  c_alpha = 0.656,
  c_beta = 437.47,
  qtl_alpha = 1,
  qtl_beta = 1,
  maf_cutoff = .05)
{

  # Randomly draw the QTL allele frequencies for all genes simulated
  qtl_af <- stats::rbeta(n_genes, qtl_alpha, qtl_beta)
  qtl_af <- (1-2*maf_cutoff) * qtl_af + maf_cutoff # Use a little y = mx+b to get the haplotypes to the right place

  # Simulate the allele frequency of coding variants corresponding to each gene
  csnp_af <- rbeta(n_genes, c_alpha, c_beta)

  # Compile the output
  haplotype_configs <- list()
  for(i in 1:n_genes){

    # Generate sQTL haplotypes for all individuals for gene i
    ssnp_haps <-
      t(
        sapply(1:n_indvs, function(x){
          sample(c(0,1), 2, replace = T, prob = c(qtl_af[i], 1-qtl_af[i]))
        })
      )

    # Recompute the QTL allele freqency based on what we got
    qtl_af[i] <- 1-sum(ssnp_haps)/(n_indvs*2)

    # Randomly draw which haplotype the coding SNP is on.
    csnp_haps <- t(sapply(1:n_indvs, function(x){
      sample(c(0,1), 2, replace = T, prob = c(1-csnp_af[i], csnp_af[i]))
    }))

    # Mark each type of haplotype combination
    haplotypes <- cbind(ssnp_haps, csnp_haps)

    hap_map <- function(numerical_hap){
      qtl1 = ifelse(numerical_hap[1], 'a', 'A')
      qtl2 = ifelse(numerical_hap[2], 'a', 'A')
      crv1 = ifelse(numerical_hap[3], 'b', 'B')
      crv2 = ifelse(numerical_hap[4], 'b', 'B')

      character_hap <- paste0(qtl1, crv1, qtl2, crv2)

      return(character_hap)
    }

    hap_symbols <- apply(haplotypes, 1, hap_map)

    output <- data.frame(indv = paste0("indv", 1:n_indvs),
                         haplotype = hap_symbols,
                         qtl_af = qtl_af[i])

    # Get rid of the cases of homozygous for the non-mutated cSNP, which will
    # be most of them.
    output <- output[!(output$haplotype %in% c("ABAB", "aBAB", "ABaB", "aBaB")),]

    # Append to the output list
    haplotype_configs[[paste0("gene", i)]] <- output

  }
  return(dplyr::bind_rows(haplotype_configs, .id = "gene"))
}

# Simulate the data that's included in this package
# set.seed(1234321)
# test_data <- simulate_haplotype_counts(1000, 5000)
# save(test_data, file = "data/test_data.rda", compress = "xz")

#Save the haplotype configurations to the data directory
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

#' Beta plotting function
#'
#' This function is simply used to plot a beta distribution, to give the user a
#' sense for what the simulated allele frequency distribution looks like.
#'
#' @param a The alpha parameter of a beta distribution
#' @param b The beta parameter of a beta distribution
#' @param ... Additional arguments to the plot function
#'
#' @export
beta_plotter <- function(a,b, ...){
  x <- seq(0, 1, by = .0005)
  y <- stats::dbeta(x, a, b)
  plot(x,y, type = "l", ...)
}


#' Beta Method-of-Moments Estimator
#'
#' A useful function for estimating the parameters for a beta distribution that
#' generated some given data. Returns estimated alpha and beta values
#'
#' @param x A vector of beta-distributed values (between 0 and 1)
#'
#' @export
beta_mom_estimator <- function(x){
  xbar <- stats::mean(x)
  vbar <- stats::var(x)

  ahat <- xbar*((xbar*(1-xbar))/vbar - 1)
  bhat <- (1-xbar)*(xbar*(1-xbar)/vbar - 1)
  return(c(ahat, bhat))
}
