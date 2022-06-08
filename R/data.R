#' Extensive dictionary for the example dataset.
#'
#' A List containing mapping of HLA ambiguity codes to all possible alleles.
#' Some ambiguities are situational and do not correspond to official MAC;
#' those are denoted by UNK#### with 4 digits numeric identifiers.
#'
#' @format A List with a 17930 entries:
#'
#'
"dict"


#' List of data.frames of allele frequencies for the example dataset.
#'
#' A List with each locus as a list element. For each locus a data.frame with
#' the alleles as rows and different regions or populations as columns.
#'
#' @format A List with 4 loci ("A", "B", "C", "DRB1").
"egf"


#' Haplotypes frequency for the example dataset.
#'
#' A dataset with two different populations (full and simplified) according to six loci haplotypes.
#'
#' @format A data frame with 4381 rows and 3 variables:
#' \describe{
#'   \item{alleles}{the haplotypes with alleles separated by "~"}
#'   \item{NA.}{frequency for the example data}
#'   \item{NA.smpl}{frequency for the simplified example data}
#'   ...
#' }
#'
"egh"


#' Sample genotypes for the example dataset.
#'
#' A dataset with sample genotypes (4 loci).
#'
#' @format A data frame with 100 rows and 10 variables:
#' \describe{
#'   \item{IND}{the individual identifier}
#'   \item{A_1}{first HLA-A allele}
#'   \item{A_2}{second HLA-A allele}
#'   \item{B_1}{first HLA-B allele}
#'   \item{B_2}{second HLA-B allele}
#'   \item{C_1}{first HLA-C allele}
#'   \item{C_2}{second HLA-C allele}
#'   \item{DRB1_1}{first HLA-DRB1 allele}
#'   \item{DRB1_2}{second HLA-DRB1 allele}
#'   \item{CRA}{region/population of the individual}
#'   ...
#' }
#'
"egt"


#' Simple dictionary for the test dataset.
#'
#' A List containing mapping of HLA ambiguity codes to all possible alleles.
#'
#' @format A List with a single entry:
#'
#'
"d_test"


#' List of data.frames of allele frequencies for the test dataset.
#'
#' A List with each locus as a list element. For each locus a data.frame with
#' the alleles as rows and different regions or populations as columns.
#'
#' @format A List with 4 loci ("A", "B", "C", "DRB1").
"f_test"


#' Haplotypes frequency for the test dataset.
#'
#' A dataset with two different populations (full and simplified) according to four loci haplotypes.
#'
#' @format A data frame with 7 rows and 3 variables:
#' \describe{
#'   \item{alleles}{the haplotypes with alleles separated by "~"}
#'   \item{01}{frequency for the test data}
#'   \item{02}{frequency for the simplified test data}
#'   ...
#' }
#'
"h_test"


#' Sample genotypes for the test dataset.
#'
#' A dataset with sample genotypes (4 loci).
#'
#' @format A data frame with 4 rows and 10 variables:
#' \describe{
#'   \item{IND}{the individual identifier}
#'   \item{A_1}{first HLA-A allele}
#'   \item{A_2}{second HLA-A allele}
#'   \item{B_1}{first HLA-B allele}
#'   \item{B_2}{second HLA-B allele}
#'   \item{C_1}{first HLA-C allele}
#'   \item{C_2}{second HLA-C allele}
#'   \item{DRB1_1}{first HLA-DRB1 allele}
#'   \item{DRB1_2}{second HLA-DRB1 allele}
#'   \item{CRA}{region/population of the individual}
#'   ...
#' }
#'
"t_test"
