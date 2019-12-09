#' Some simulated haplotype count data from 1000 individuals and 500 genes
#'
#' This file is the result of the simulate_data.R script, saved for faster
#' bootstrap testing. When simulated, there is no dependence between the coding
#' allele and the sqtl allele, so the result of the test should be null
#'
#' @format A data frame with 500 rows and 10 columns
#' \describe{
#'     \item{gene}{simulated gene name}
#'     \item{sqtl_af}{The frequency of the higher included sQTL haplotype}
#'     \item{[Aa][Bb][Aa][Bb]}{The counts of each haplotype configuration across the 1000 individuals}
#'     ...
#'}
"test_data"
