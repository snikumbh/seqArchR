#' @title Plot a given position weight matrix as a heatmap using ggplot2
#'
#' @description The given matrix is plotted as a heatmap using \code{ggplot2}'s
#' \code{geom_tile}.
#'
#' @param pwm_mat Matrix (usually a PWM, but can be non-normalized/any matrix)
#' to be represented as a heatmap.
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the
#' positions are labeled from 1 to the length of the sequences.
#' @param pdf_name Name of the file which will be saved as PDF.
#' @param fixed_coord Set this to TRUE to use a fixed aspect ratio for the
#' plot. Default is FALSE.
#'
#' @return A ggplot object so you can simply call \code{print} or \code{save}
#' on it later. If \code{pdf_name} is given, it is also saved and the
#' \code{ggplot} object returned.
#'
#' @export
#'
#' @family visualization functions
#'
#' @seealso \code{\link{plot_ggseqlogo}} for plotting PWMs as sequence logos
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import ggseqlogo
#'
#' @examples
#'
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' pwm <- seqArchR::make_dinuc_PWMs(get_clBasVec_m(res,iter=1)[,1],
#'                         add_pseudo_counts = FALSE)
#'
#' plot_ggheatmap(pwm_mat = pwm, fixed_coord = TRUE)
#'
plot_ggheatmap <- function(pwm_mat, pos_lab = NULL, pdf_name = NULL,
                            fixed_coord = FALSE){
    if(is.null(pos_lab)) pos_lab <- set_default_pos_lab2(pwm_mat)

    check_vars(pwm_mat, pos_lab)

    ##
    ## Convert pwm_mat to df, heatmap by ggplot-way
    pwm_mat_df <- as.data.frame(pwm_mat)
    pwm_mat_df$Nucleotides <- rownames(pwm_mat_df)
    colnames(pwm_mat_df) <- c(pos_lab, "Nucleotides")
    pwm_mat_df_for_ggheatmap <- melt(pwm_mat_df, id.vars = c("Nucleotides"),
                                    variable.name = "positions")

    ##
    p1 <- ggplot2::ggplot(data = pwm_mat_df_for_ggheatmap, mapping = aes(
        x = positions,
        ## Here, 'positions' is the column_name, see previous statement.
        ## Do not change it to position_labels
        y = Nucleotides,
        fill = value
    )) +
        ggplot2::geom_tile() +
        ggplot2::theme_bw() +
        ggplot2::xlab(label = element_blank()) +
        ggplot2::scale_fill_gradient2(name = "", low = "white",
                                mid = "white", high = "#012345") +
        ggplot2::theme(legend.position = "top",
            legend.justification = "center",
            legend.margin = margin(0,-1,0,0),
            axis.text.x = element_text(size = rel(0.8), angle = 90,
                            hjust = 1, vjust = 0.5),
            plot.margin = margin(0,0,0,0)
        )

    p1 <- fix_coord(p1, nPos = length(pos_lab), method = "heatmap",
        fixed_coord = fixed_coord)

    ##
    if (!is.null(pdf_name)) {
        if (file.exists(pdf_name)) {
            warning("File exists, will overwrite", immediate. = TRUE)
        }
        ggplot2::ggsave(filename = pdf_name, plot = p1, device = "pdf",
                        width = 20, height = 2.5)
    }

    return(p1)
}
## =============================================================================

#' @title Visualize a given (PWM) matrix as a sequence logo.
#'
#' @param pwm_mat Matrix (usually a PWM, but can be any non-normalized matrix)
#' to be represented as a sequence logo. Rownames must be letters.
#'
#' @param method For \code{ggseqlogo}; either of 'custom', 'bits', or
#' 'probability'. Default is 'bits'.
#'
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the
#' positions are labeled from 1 to the length of the sequences.
#'
#' @param pdf_name Name of the file which will be saved as PDF.
#'
#' @param bits_yax Specify 'full' if the information content y-axis limits
#' should be 0-2 or 'auto' for a suitable limit. The 'auto' setting adjusts
#' the y-axis limits according to the maximum information content of the
#' sequence logo. Default is 'full'.
#'
#' @param fixed_coord Set this to TRUE to use a fixed aspect ratio for the
#' plot. Default is FALSE.
#'
#' @return A ggplot2 object so you can simply call \code{print} or \code{save}
#' on it later. If \code{pdf_name} is given, it is also saved in addition to
#' returning the ggplot object.
#'
#' @export
#'
#' @family visualization functions
#'
#' @seealso \code{\link{plot_ggheatmap}} for plotting PWMs as heatmaps,
#' \code{\link{plot_ggseqlogo_of_seqs}} for visualizing a collection of
#' sequences by their sequence logo.
#'
#' @import ggplot2
#' @import ggseqlogo
#'
#' @examples
#'
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' pwm <- seqArchR::make_dinuc_PWMs(get_clBasVec_m(res,iter=1)[,1],
#'                         add_pseudo_counts = FALSE)
#'
#' plot_ggseqlogo(pwm_mat = pwm, fixed_coord = TRUE)
#'
plot_ggseqlogo <- function(pwm_mat, method = "bits", pos_lab = NULL,
    pdf_name = NULL, bits_yax = "full", fixed_coord = FALSE){
    ##
    if(is.null(pos_lab)) pos_lab <- set_default_pos_lab2(pwm_mat)
    check_vars(pwm_mat, pos_lab)
    ##
    p1 <- ggplot() +
        ggseqlogo::geom_logo(pwm_mat, method = method, seq_type = "dna") +
        ggplot2::theme_linedraw() +
        ggplot2::theme(axis.text.x = element_text(size = rel(0.8),
            angle = 90, hjust = 1, vjust=0.5),
            axis.text.y = element_text(size = rel(0.8)),
            panel.grid = element_blank()) +
        ggplot2::xlab(label = "Positions") +
        ##
        ggplot2::scale_x_continuous(breaks = seq_len(ncol(pwm_mat)),
            labels = pos_lab,
            expand = expansion(mult = c(0, 0))) #+
        ##
        # ggplot2::theme(axis.text.x = element_text(size = rel(0.9),
        #     angle = 90, hjust = 1),
        #     axis.text.y = element_text(size = rel(0.9)))
    ##
    p1 <- add_lims_lab(p1, method, bits_yax)
    p1 <- fix_coord(p1, nPos = length(pos_lab),method = method,
                    fixed_coord = fixed_coord)
    ##
    if (!is.null(pdf_name)) {
        if (file.exists(pdf_name)) {
            warning("File exists, will overwrite", immediate. = TRUE)
        }
        ggsave(filename = pdf_name, plot = p1, device = "pdf",
                width = 25, height = 2.5)
    }
    return(p1)
}
## =============================================================================


set_default_pos_lab2 <- function(pwm_mat){
    pos_lab <- seq_len(ncol(pwm_mat))
}

add_lims_lab <- function(p1, method, yax){
    ## add limits to ggplot2 object
    if(method == "prob"){
        p1 <- p1 + ggplot2::ylab(label = "Probability")
    }else if(method == "bits"){
        p1 <- p1 + ggplot2::ylab(label = "Bits")
        if(yax == "full") p1 <- p1 + ggplot2::ylim(0, 2)
    }
    return(p1)
}

check_vars <- function(pwm_mat, pos_lab){
    if (!is.matrix(pwm_mat)) {
        stop("Expecting a matrix with 4 rows")
    }
    if (sum(dim(pwm_mat)) == 2 && is.na(pwm_mat)) {
        stop("Empty matrix")
    }
    if (!(nrow(pwm_mat) == 4)) {
        stop("Expecting a matrix with 4 rows corresponding to DNA chars ",
            "'A', 'C', 'G', 'T'")
    }
    ##
    if (length(pos_lab) < ncol(pwm_mat)) {
        stop("Inadequate position labels supplied",
            ncol(pwm_mat) - length(pos_lab)
        )
    }
    ##
    if (length(pos_lab) > ncol(pwm_mat)) {
        stop("Overabundant position labels supplied",
            length(pos_lab) - ncol(pwm_mat)
        )
    }
}

