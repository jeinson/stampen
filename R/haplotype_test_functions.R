#' Reshape data and classify the haplotypes according to the model
#'
#' This function takes the haplotype configurations in the standard input format
#' and reshapes them for more intuitive processing. See the sample dataset
#'
#' @param dataset A data.frame or tibble. The input dataset. See sample data for format
#' @param beta_config Vector. The values assigned to beta for each haplotype possibility. This can be changed to tweak the hypothesis.
#' Values are provided in the package as `beta_config_eqtl` and `beta_config_sqtl`, depending on the hypothesis
#' you are testing.
#' @param qtl_allele_frequency_column String. The name of the column in the input data with QTL allele frequencies. Default: `sqtl_af`
#' the modified penetrance hypothesis.
#' @param gene_name_column String. The name of the column in the input data with the gene name. Default: `gene`
#'
#' @export
#'
characterize_haplotypes <- function(dataset, beta_config,
                                    qtl_allele_frequency_column = "sqtl_af",
                                    gene_name_column = "gene"){
  x <- dataset

  # Make sure to use the correct beta configurations
  beta_config <- beta_config[names(beta_config) %in% names(x)]

  # Rename the input data
  x <- dplyr::rename(x, "sqtl_af" = qtl_allele_frequency_column)
  x <- dplyr::rename(x, "gene" = gene_name_column)
  x <- dplyr::select(x, "sqtl_af", "gene", names(beta_config))

  # Reshape and spread haplotypes out to one per line
  x <- tidyr::gather(x, hap_config, num, -"gene", -"sqtl_af", convert = F, factor_key = F)
  x <- x[x$num != 0,]
  x <- as.data.frame(lapply(x, rep, x$num))
  x <- x[,-4]
  x$beta <- beta_config[as.character(x$hap_config)] # Fuck you R
  x$hom <- grepl("A.A.", x$hap_config) | grepl("a.a.", x$hap_config)

  x$exp_beta <- (x$sqtl_af)^2 / (x$sqtl_af^2 + (1-x$sqtl_af)^2)

  # Switch out exp_beta for heterozygotes
  for(i in 1:nrow(x)){
    if(!x$hom[i]){
      x$exp_beta[i] <- .5
    }
  }

  return(x)
}


#' Calculate epsilon
#'
#' This function calculates the value of epsilon, dependent on the value of
#' beta and expected value of beta given each haplotype combination observed in
#' a set of haplotypes based on phased genetic data.
#'
#' @param dataset ...
#' @param af_weighted logical. If TRUE, weight the epsilon calculation by the sQTL minor allele frequency
#'
#' @export
#'
calculate_epsilon <- function(dataset, beta_config, af_weighted = FALSE, ...){
  x <- characterize_haplotypes(dataset, beta_config, ...)

  # Do the math
  ifelse(af_weighted,
         mean( (x$beta - x$exp_beta) / x$exp_beta) * x$sqtl_af^2,
         mean( (x$beta - x$exp_beta) / x$exp_beta)
  )
}


#' Calculate a bootstrapped p-value
#'
#' This function takes the input data and returns a p-value for the null
#' hypothesis that haplotype configurations which influence coding variant
#' penetrance are not significantly enriched in the population. The p-value is
#' estimated through a bootstrapping procedure, and a confidence interval
#' is returned.
#'
#' @param haps The output of `characterize_haplotypes`, or another haplotype calling tool. There is one line per individual
#' @param beta_config The values assigned to beta for each haplotype possibility. This can be changed to tweak
#' @param af_weighted logical. If TRUE, weight the epsilon calculation by the sQTL minor allele frequency
#' @param B The number of bootstraps to perform. Default: 1000
#'
#' @export
#'
bootstrap_test <- function(haps, exp_beta = "exp_beta", af_weighted = FALSE,
                           B = 1000, ...){
  x <- haps
  epsilon <-
    ifelse(af_weighted,
           mean((as.numeric(x$beta) - x[[exp_beta]]) / x[[exp_beta]]) * x$sqtl_af^2,
           mean((as.numeric(x$beta) - x[[exp_beta]]) / x[[exp_beta]])
    )

  p_b = rep(0, B)
  for(b in 1:B){
    ix <- sample(1:nrow(x), nrow(x), replace = T)
    x_b <- x[ix,]

    # Do the math
    p_b[b] <-
      ifelse(af_weighted,
             mean((as.numeric(x_b$beta) - x_b[[exp_beta]]) / x_b[[exp_beta]]) * x_b$sqtl_af^2,
             mean((as.numeric(x_b$beta) - x_b[[exp_beta]]) / x_b[[exp_beta]])
      )
  }

  p_b <- sort(p_b)
  p_Ho <- 2 * min(sum(p_b < 0) / B, sum(p_b > 0) / B)
  lower <- p_b[ceiling(B * .05)]
  upper <- p_b[ceiling(B * .95)]

  n_haplotypes <- nrow(x)

  return(c(bootstrap_p = p_Ho,
           epsilon = epsilon,
           lower = lower,
           upper = upper,
           n_haplotypes = n_haplotypes))

}

#' Calculate a poisson-binomial p-value
#'
#' This function takes the input data and returns a p-value for the null
#' hypothesis that haplotype configurations which influence coding variant
#' penetrance are not significangtly enriched in the population. The p-value
#' is calculated analytically using the poisson-binomial distribuion, which
#' is the distribution over the sum of n independent and NOT identically
#' distributed bernoulli random variables.
#'
#' @param dataset
#' @param beta_config
#'
#' @export
#'
poison_binomial_test <- function(dataset, beta_config, ...){
  x <- characterize_haplotypes(dataset, beta_config, ...)
  epsilon <- mean((as.numeric(x$beta) - x$exp_beta) / x$exp_beta)

  p1 <- poisbinom::ppoisbinom(sum(x$beta), x$exp_beta, lower_tail = T)
  p2 <- poisbinom::ppoisbinom(sum(x$beta), x$exp_beta, lower_tail = F)
  p <- 2 * min(p1, p2)

  return(c(poison_binomial_p = p,
           epsilon = epsilon,
           n_haplotypes = nrow(x)
  )
  )
}
