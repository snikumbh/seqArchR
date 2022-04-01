# @title
# Convert DNA sequences from DNAStringSet object to a numeric matrix
#
# @description This function converts a DNAStringSet object into a numeric
# matrix
# @param seqs The sequences to converted to matrix
# @param position_labels Labels for the axis denoting sequence positions
# @param annClusters For later, when clusters are to be annotated
# @param sinuc_or_dinuc For later, maybe
.seqs_to_mat <- function(seqs, pos_lab, annClusters = NULL,
                            sinuc_or_dinuc = "sinuc") {
    nSeqs <- length(seqs)
    nPos <- length(pos_lab)
    ## handling single- and di-nucleotides separately
    if (sinuc_or_dinuc == "sinuc") {
        nuc_list <- unlist(lapply(
                    seq_along(seqs),
                    function(x) {str_seq <- seqs[x]
                        nucleotides <- unlist(strsplit(str_seq, split = NULL))
                    }))
        nuc_list[which(nuc_list == "A")] <- 1
        nuc_list[which(nuc_list == "C")] <- 2
        nuc_list[which(nuc_list == "G")] <- 3
        nuc_list[which(nuc_list == "T")] <- 4
        nuc_list_num <- as.numeric(nuc_list)
        nuc_mat <- matrix(nuc_list_num, byrow = TRUE, ncol = nPos, nrow = nSeqs)
    } else if (sinuc_or_dinuc == "dinuc") {
        ## handling dinucleotides
        ## not needed any more?!
        stop("Nothing for dinuc here. Unrelated!")
    }
    return(nuc_mat)
}
## =============================================================================


#' @title
#' Visualize raw DNA sequences as an image
#'
#' @description This function plots the collection of sequences as an image
#' matrix.
#'
#' @param seqs The sequences as a DNAStringSet object.
#' @param pos_lab The labels to be used for the sequence positions.
#' Default: Sequence positions are labeled from 1 to the length of the
#' sequences.
#' @param xt_freq The x-axis tick frequency. Expects a positive integer less
#' than the length of the sequences. Default is 5.
#' @param yt_freq The y-axis tick frequency. Expects a positive integer less
#' than number of sequences. Default is 100.
#' @param use_col A vector of four colors used for the DNA bases A, C, G,
#' and T (in that order).
#' @param add_legend Logical. Whether legend should be added to the plot.
#' Default is TRUE.
#' @param use_legend A character vector of letters in the input sequences.
#' Default is \code{\link[Biostrings]{DNA_BASES}}, used for DNA sequences.
#' @param save_fname Specify the filename (with extension) for saving the
#' plot to disk.
#' @param file_type Specify the file type, namely PNG, JPEG, TIFF.
#' @param f_width Specify the width for the plot. This depends on the length of
#' sequences.
#' @param f_height Specify the height for the plot. This depends on the number
#' of sequences.
#' @param f_units Specify the units in which the height and width are given.
#'
#' @importFrom Biostrings width DNA_BASES
#'
#' @return Nothing returned to the R interpreter.
#' @family visualization functions
#' @importFrom grDevices png dev.off
#' @importFrom graphics axis image legend par
#'
#' @examples
#'
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' # Image matrix of sequences in the input order
#' viz_seqs_acgt_mat(seqs = seqs_str(res))
#'
#' # Image matrix of sequences ordered by the clustering from seqArchR
#' use_seqs <- seqs_str(res, iter = NULL, cl = NULL, ord = TRUE)
#' viz_seqs_acgt_mat(seqs = use_seqs)
#'
#' # Image matrix of sequences belonging to a single cluster
#' use_seqs <- seqs_str(res, iter = 2, cl = 2)
#' viz_seqs_acgt_mat(seqs = use_seqs)
#'
#' @export
viz_seqs_acgt_mat <- function(seqs, pos_lab = NULL,
                                    xt_freq = min(length(pos_lab), 5),
                                    yt_freq = min(length(seqs), 100),
                                    use_col = c("darkgreen", "blue",
                                            "orange", "red"),
                                    add_legend = TRUE,
                                    use_legend = Biostrings::DNA_BASES,
                                    save_fname = NULL,
                                    file_type = "PNG",
                                    f_width = 450, f_height = 900,
                                    f_units = "px"
                                    ) {
    ##
    pos_lab <- set_default_pos_lab2(Biostrings::DNAString(seqs[1]), pos_lab)
    if(xt_freq <= 0 || xt_freq > length(pos_lab)){
        warning("Expected positive integer (< length of sequences) for",
            "xt_freq. Reverting to default value")
        xt_freq <- min(length(pos_lab), 5)
    }
    if(yt_freq <= 0 || yt_freq > length(seqs)){
        warning("Expected positive integer (< number of sequences) for",
            "yt_freq. Reverting to default value")
        yt_freq <- min(length(seqs), 100)
    }
    ##
    nSeqs <- length(seqs)
    nPos <- length(pos_lab)
    seq_mat <- .seqs_to_mat(seqs = seqs, pos_lab = pos_lab)
    seq_mat <- seq_mat[rev(seq_len(nSeqs)),]
    ##
    if(!is.null(save_fname)){
        if(file_type == "PNG" || file_type == "png"){
            grDevices::png(filename = save_fname, width = f_width,
                height = f_height, units = f_units, bg = "white")
        }else{
            ## TODO
        }
    }else{ ## do nothing
    }

    use_xtick_labs <- set_xtick_labels(pos_lab, xt_freq)
    ##
    ytick_names <- rev(seq(yt_freq, nSeqs, by = yt_freq))
    ytick_loc <- 1 + nSeqs - c(rev(seq(yt_freq, nSeqs, by = yt_freq)))

    par(xpd = TRUE, mar = par()$mar + c(0,0,0,2.5))
    graphics::image(x = seq_len(nPos), y = seq_len(nSeqs),
                z = t(seq_mat), col = use_col, useRaster = TRUE,
                ylab = paste0("Sequences (n = ", nSeqs, ")"),
                xlab = "Positions", axes = FALSE)
    axis(side = 1, at = use_xtick_labs$breaks,
            labels = use_xtick_labs$labels, las = 2)
    axis(side = 2, at = ytick_loc, labels = ytick_names, las = 2)
    ##
    if(add_legend){
        legend("topright", inset = c(-0.07, 0), col = use_col,
            legend = use_legend, horiz = FALSE,
            pch = 15, pt.cex = 1.2, bty = "n")
    }

    par(mar=c(5, 4, 4, 2) + 0.1)
    if(!is.null(save_fname)){
        dev.off()
    }
}
