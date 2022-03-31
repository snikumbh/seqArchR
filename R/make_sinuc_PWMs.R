#' @title Make a PWM-resembling matrix out of a given n-vector
#'
#' @description The given matrix (or simply a vector) is reshaped to have four
#' rows for four nucleotides and a relevant number of columns.
#'
#' @param vec A vector that will be reshaped into a PWM matrix of DNA
#' sequences. Note that the matrix is formed by row.
#' @param add_pseudo_counts Logical, taking values TRUE or FALSE, specifying
#' whether or not pseudocounts are added to the matrix.
#' @param scale Logical, taking values TRUE or FALSE, specifying whether or
#' not the matrix is scaled column-wise, i.e., all columns summed to 1.
#' @param sinuc Logical. Specify TRUE for mononucleotides (default), FALSE to
#' for dinucleotides.
#'
#' @return A PWM. If sinuc is `TRUE`, the PWM has 4 rows corresponding to the
#' 4 nucleotides (A, C, G, T) and the relevant number of columns (i.e.,
#' number of elements in given vector/4).
#' If dinucleotide is selected, by setting `sinuc` to `FALSE`, the PWM has
#' 16 rows corresponding to the dinucleotide combinations of the four
#' nucleotides (A, C, G, T) and the relevant number of columns (i.e.,
#' number of elements in given vector/16).
#'
#' @examples
#'
#' ## Mononucleotides case
#' ## Make a dummy PWM of dimensions 4 * 10 from a vector
#' vec <- runif(4*10)
#' pwm <- seqArchR::make_PWMs(vec = vec, add_pseudo_counts = FALSE)
#'
#' ## Dinucleotides case
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' pwm <- seqArchR::make_PWMs(get_clBasVec_m(res,iter=1)[,1],
#'                         add_pseudo_counts = FALSE, sinuc = FALSE)
#'
#' @export
#'
make_PWMs <- function(vec, add_pseudo_counts = TRUE, scale = TRUE,
                      sinuc = TRUE) {
    ##
    if(add_pseudo_counts) vec <- add_pseudo_count(vec)
    ##
    rnames <- Biostrings::DNA_BASES
    if(!sinuc) rnames <- get_dimers_from_alphabet(Biostrings::DNA_BASES)
    ##
    this_mat <- make_matrix(vec, use_rnames = rnames)
    ##
    if (scale) {
        this_mat <- scale_matrix(this_mat)
    }
    ##
    if(sinuc) return(this_mat)
    ## sinuc == FALSE reaches here
    return_mat <- collapse_into_sinuc_matrix(dinuc_mat = this_mat,
                                             feat_names = rnames)
    return(return_mat)
}


add_pseudo_count <- function(vec, pseudo_count = 10^-5){
    vec <- vec + pseudo_count
}

make_matrix <- function(vec, use_rnames){
    this_mat <- t(matrix(vec, ncol = length(use_rnames), byrow = FALSE))
    rownames(this_mat) <- use_rnames
    this_mat
}

scale_matrix <- function(this_mat){
    scaled <- base::sweep(this_mat, 2, colSums(this_mat), "/")
    scaled
}

# @title Similarly to the PWM-like matrix for mononucleotides, make one for
#  dinucleotides
#
# @description This function converts the basis matrix with basis vectors
# of dinucleotide information into matrix of dimension
# 16 x (sequence_length) for visualization.
#
# @param vec Column vector of the basis matrix
# @param add_pseudo_counts Whether pesudocounts are to be added. TRUE or FALSE.
# @param scale Whether to perform per position scaling of the matrix. TRUE or
# FALSE
#
# @return A (PWM) matrix with 16 rows corresponding to the dinucleotide
# combinations of the four nucleotides (A, C, G, T) and the relevant number
# of columns (i.e., number of elements in given vector/16)
#
# @examples
#
# res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#          package = "seqArchR", mustWork = TRUE))
#
# pwm <- seqArchR::make_dinuc_PWMs(get_clBasVec_m(res,iter=1)[,1],
#                         add_pseudo_counts = FALSE)
#
# @export
#
# make_dinuc_PWMs <- function(vec, add_pseudo_counts = TRUE, scale = TRUE) {
#     # column vector is expected as input
#     dinuc <- get_dimers_from_alphabet(Biostrings::DNA_BASES)
#     if(add_pseudo_counts) vec <- add_pseudo_count(vec)
#     this_mat <- make_matrix(vec, use_rnames = dinuc)
#     if (scale) this_mat <- scale_matrix(this_mat)
#     return_mat <- collapse_into_sinuc_matrix(dinuc_mat = this_mat,
#                                             feat_names = dinuc)
#     return(return_mat)
# }





## Collapse the dinucleotide matrix to single nucleotide
collapse_into_sinuc_matrix <- function(dinuc_mat, feat_names){
    # Collapse the matrix of dinuc features into sinuc features
    # When V is (features x #seqs)

    sinuc_feat_mat <-  matrix(rep(0, times=4*ncol(dinuc_mat)),
                            nrow=4, ncol=ncol(dinuc_mat))
    rownames(sinuc_feat_mat) <- c('A', 'C', 'G', 'T')

    useMat <- dinuc_mat
    splitted_names <- strsplit(feat_names, "")

    for(i in seq_along(splitted_names)){
        for (j in seq_len(ncol(dinuc_mat)-1)){
            sinuc_feat_mat[splitted_names[[i]][1], j] <-
                sinuc_feat_mat[splitted_names[[i]][1], j] + 10*useMat[i,j]
            ##
            ##
            sinuc_feat_mat[splitted_names[[i]][2], j+1] <-
                sinuc_feat_mat[splitted_names[[i]][2], j+1] + 10*useMat[i,j]
        }
    }
    scaled_mat <- base::sweep(sinuc_feat_mat, 2, colSums(sinuc_feat_mat),
        "/")
    return(scaled_mat)
}
