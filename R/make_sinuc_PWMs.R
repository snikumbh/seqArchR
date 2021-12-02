#' @title Make a PWM-resembling matrix out of a given n-vector
#'
#' @description The given matrix (or simply a vector) is reshaped to have four
#' rows for four nucleotides and a relevant number of columns.
#'
#' @param mat Actually a vector that will be reshaped into a (PWM)
#' matrix of DNA sequences.
#' @param add_pseudo_counts Logical, taking values TRUE or FALSE, specifying
#' whether or not pseudocounts are added to the matrix.
#' @param scale Logical, taking values TRUE or FALSE, specifying whether or
#' not the matrix is scaled column-wise, i.e., all columns summed to 1.
#'
#' @return A (PWM) matrix with 4 rows corresponding to the 4 nucleotides (A, C,
#' G, T) and the relevant number of columns (i.e., number of elements in given
#' vector/4)
#'
#' @export
#' 
make_sinuc_PWMs <- function(mat, add_pseudo_counts = TRUE, scale = TRUE) {
    ##
    sinuc <- c("A", "C", "G", "T")
    if (add_pseudo_counts) {
        mat <- mat + 10^-5
    }
    this_mat <- t(matrix(mat, ncol = length(sinuc),
                                    byrow = FALSE))
    rownames(this_mat) <- sinuc
    if (scale) {
        scaled <- this_mat
        scaled <- base::sweep(this_mat, 2, colSums(this_mat),
                                "/")
        return(scaled)
    } else {
        return(this_mat)
    }
}


#' @title Similarly to the PWM-like matrix for mononucleotides, make one for 
#'  dinucleotides 
#' 
#' @description This function converts the basis matrix with basis vectors 
#' of dinucleotide information into matrix of dimension 
#' 16 x (sequence_length) for visualization.
#' 
#' @param vec Column vector of the basis matrix
#' @param add_pseudo_counts Whether pesudocounts are to be added. TRUE or FALSE.
#' @param scale Whether to perform per position scaling of the matrix. TRUE or 
#' FALSE
#' 
#' @return A (PWM) matrix with 16 rows corresponding to the dinucleotide 
#' combinations of the four nucleotides (A, C, G, T) and the relevant number 
#' of columns (i.e., number of elements in given vector/16)
#' 
#' @export
#' 
make_dinuc_PWMs <- function(vec, add_pseudo_counts = TRUE, scale = TRUE) {
    
    # column vector is expected as input
    dna_alphabet <- c("A", "C", "G", "T")
    dinuc <- do.call(paste0, expand.grid(dna_alphabet, dna_alphabet))
    if (add_pseudo_counts) {
        vec <- vec + 10^-5
    }
    this_mat <- t(matrix(vec, ncol = length(dinuc),
                                    byrow = FALSE))
    rownames(this_mat) <- dinuc
    if (scale) {
        this_mat <- base::sweep(this_mat, 2, colSums(this_mat), "/")
    }
    return_mat <- collapse_into_sinuc_matrix(dinuc_mat = this_mat, 
                                            feat_names = dinuc)
    return(return_mat)
}


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
