#' Some simulated haplotype count data from 1000 individuals and 500 genes
#'
#' This file is the result of the simulate_data.R script, saved for faster
#' bootstrap testing. When simulated, there is no dependence between the coding
#' allele and the sqtl allele, so the result of the test should be null
#'
#' @format A data frame with about 15,000 rows and 4 columns
#' \describe{
#'     \item{gene}{simulated gene name}
#'     \item{indv}{simulated individual name}
#'     \item{haplotype}{The simulated haplotype with 'a' and 'A' being high and low included alleles, and 'b' and 'B' being rare and common coding variant alleles.}
#'     \item{qtl_af}{The frequency of the 'penetrance driving' QTL qllele. i.e. lower splicing or higher expression allele.}
#'     ...
#'}
"test_data"

#' beta_config_eqtl
#'
#' ...
#'
#' @format ...
#'
"beta_config_eqtl"

#' beta_config_sqtl
#'
#' ...
#'
#' @format ...
#'
"beta_config_sqtl"
