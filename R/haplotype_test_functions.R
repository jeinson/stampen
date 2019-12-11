#' Reshape data and classify the haplotypes according to the model
#'
#' This function takes the haplotype configurations in the standard input format
#' and reshapes them for more intuitive processing. See the sample dataset
#'
#' @param dataset The input dataset. See sample data for format
#' @param beta_config The values assigned to beta for each haplotype possibility. This can be changed to tweak
#' the modified penetrance hypothesis. d
#'
#' Internal function, not exported
characterize_haplotypes <- function(dataset, beta_config){
  x <- tidyr::gather(dataset, hap_config, num, -gene, -sqtl_af, convert = F, factor_key = F)
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
#' @param dataset
calculate_epsilon <- function(dataset, beta_config){
  x <- characterize_haplotypes(dataset, beta_config)

  # Do the math
  mean( (x$beta - x$exp_beta) / x$exp_beta)
}

#' Calculate epsilon (alternative method)
#'
#' This function calculates the value of epsilon, dependent on the value of
#' beta and expected value of beta given each haplotype combination observed in
#' a set of haplotypes based on phased genetic data. In this function,
#' subtract 0.5 from beta, instaed of E[beta(g)]
#'
#' @param dataset
calculate_epsilon_alt <- function(dataset, beta_config){
  x <- characterize_haplotypes(dataset, beta_config)

  # Do the math
  mean( (x$beta - 0.5) / abs(x$beta - x$exp_beta) )
}



#' Calculate a bootstrapped p-value
#'
#' This function takes the input data and returns a p-value for the null
#' hypothesis that haplotype configuations which influence coding variant
#' penetrance are not significantly enriched in the population.
#'
#' @param dataset The input dataset. See sample data for format
#' @param beta_config The values assigned to beta for each haplotype possibility. This can be changed to tweak
#' @param B The number of bootstraps to perform. Default: 1000

bootstrap_test <- function(dataset, beta_config, B = 1000){
  x <- characterize_haplotypes(dataset, beta_config)

  p_b = rep(0, B)
  for(b in 1:B){
    ix <- sample(1:nrow(x), nrow(x), replace = T)
    x_b <- x[ix,]

    # Do the math
    p_b[b] <- mean((as.numeric(x_b$beta) - x_b$exp_beta) / x_b$exp_beta)
  }

  p_b <- sort(p_b)
  p_Ho <- 2 * min(sum(p_b < 0) / B, sum(p_b > 0) / B)
  lower <- p_b[ceiling(B * .05)]
  upper <- p_b[ceiling(B * .95)]

  return(c(bootstrap_p = p_Ho, lower = lower, upper = upper))

}

#' Calculate a bootstrapped p-value (alternative)
#'
#' This function takes the input data and returns a p-value for the null
#' hypothesis that haplotype configuations which influence coding variant
#' penetrance are not significantly enriched in the population.
#'
#' Here, subtract 0.5 from beta instead of e_beta
#'
#' @param dataset The input dataset. See sample data for format
#' @param beta_config The values assigned to beta for each haplotype possibility. This can be changed to tweak
#' @param B The number of bootstraps to perform. Default: 1000

bootstrap_test_alt <- function(dataset, beta_config, B = 1000){
  x <- characterize_haplotypes(dataset, beta_config)

  p_b = rep(0, B)
  for(b in 1:B){
    ix <- sample(1:nrow(x), nrow(x), replace = T)
    x_b <- x[ix,]

    # Do the math
    p_b[b] <- mean( (x_b$beta - 0.5) / abs(x_b$beta - x_b$exp_beta) )
  }

  p_b <- sort(p_b)
  p_Ho <- 2 * min(sum(p_b < 0) / B, sum(p_b > 0) / B)
  lower <- p_b[ceiling(B * .05)]
  upper <- p_b[ceiling(B * .95)]

  return(c(bootstrap_p = p_Ho, lower = lower, upper = upper))

}

