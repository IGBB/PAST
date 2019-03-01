#' Genes data
#'
#' A dataset containing the assigned genes
#'
#' @format A data frame with 183 rows and 5 variables:
#' \describe{
#'   \item{Chromosome}{The chromosome on which the SNP associated with the gene is found.}
#'   \item{Position}{The physical location, in base pairs from the start of the chromosome, where the SNP associated with the gene is found. This may vary depending on the reference sequence version which was used to align the SNPs.}
#'   \item{Gene}{The name of the gene}
#'   \item{Effect}{The phenotypic effect that the allele had on the trait. The strongest effect is assigned to the gene.}
#'   \item{P_value}{The p-value from the T-test of each individual marker-trait association for the SNP associated with the gene.}
#'   ...
#' }
"genes"

#' A dataframe containing the necessary information
#' to create the rugplots in the next step
#'
#' @format A data frame with 17 rows and 9 variables:
#' \describe{
#'   \item{pathway_id}{The pathway ID from the pathways database}
#'   \item{NESrank}{A variable that should be renamed to pathway number}
#'   \item{gene_id}{The physical location, in base pairs from the start of the chromosome, where the SNP is found. This may vary depending on the reference sequence version which was used to align the SNPs.}
#'   \item{rank}{The rank of the gene in the data by its effect on the trait.}
#'   \item{running_enrichment_score}{enrichment score (ES) for the pathway is the maximum deviation from zero}
#'   \item{pathway_name}{The pathway name from the pathways file}
#'   \item{pvalue}{he probability of seeing the test statistic (in this case, the peak running enrichment score) by chance}
#'   \item{FDR}{The FDR of the pathway (not really)}
#'   \item{qvalue}{The FDR adjusted p-value is the q-value. If set at 0.05, this implies that 5% of significant tests will result in false positives}
#'   ...
#' }
"rugplots_data"