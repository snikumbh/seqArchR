#' seqArchR: A package for de novo discovery of different sequence
#' architectures
#'
#' Given a set of DNA sequences, \code{seqArchR} enables unsupervised
#' discovery of _de novo_ clusters with characteristic sequence
#' architectures characterized by position-specific motifs or composition
#' of stretches of nucleotides, e.g., CG-richness, etc.
#'
#' The seqArchR package provides three categories of important functions:
#' related to data preparation and manipulation, performing non-negative
#' matrix factorization, performing clustering, and visualization-related
#' functions.
#'
#' @importFrom methods is
#'
#' @section Functions for data preparation and manipulation:
#' \itemize{
#' \item \code{\link{prepare_data_from_FASTA}}
#' \item \code{\link{get_one_hot_encoded_seqs}}
#' }
#'
#'
#' @section Functions for visualizations:
#' \itemize{
#' \item \code{\link{plot_arch_for_clusters}}
#' \item \code{\link{plot_ggseqlogo_of_seqs}}
#' \item \code{\link{viz_bas_vec}}
#' \item \code{\link{viz_seqs_acgt_mat}}
#' \item \code{\link{viz_pwm}}
#' }
#'
#' @docType package
#' @name seqArchR
NULL


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# if(getRversion() >= "2.15.1")  utils::globalVariables(c("seqArchRconfig"))
