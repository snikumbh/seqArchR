#' @title Visualize a position weight matrix as a heatmap or sequence logo
#'
#' @description The given position weight matrix is plotted as a heatmap
#' or sequence logo
#'
#' @param pwm_mat Matrix (usually a PWM, but can be any non-normalized matrix)
#' to be visualized. Rownames must be letters.
#' @param method Character. Set this to 'heatmap' when plotting a heatmap, else
#' you can set it to either of 'custom', 'bits', or 'probability' when you
#' wish to visualize it as a sequence logo. Default is 'heatmap'.
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the
#' positions are labeled from 1 to the length of the sequences.
#' @param pdf_name Name of the file which will be saved as PDF.
#' @param fixed_coord Set this to TRUE to use a fixed aspect ratio for the
#' plot. Default is FALSE.
#' @param bits_yax Specify 'full' if the information content y-axis limits
#' should be 0-2 or 'auto' for a suitable limit. The 'auto' setting adjusts
#' the y-axis limits according to the maximum information content of the
#' sequence logo. Default is 'full'.
#'
#' @return A ggplot object so you can simply call \code{print} or \code{save}
#' on it later. If \code{pdf_name} is given, it is also saved and the
#' \code{ggplot2} object returned.
#'
#' @export
#'
#' @family visualization functions
#'
#' @seealso \code{\link{plot_ggseqlogo_of_seqs}} for visualizing a collection of
#' sequences by their sequence logo.
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot theme_bw geom_tile xlab scale_fill_gradient2
#' @importFrom ggplot2 theme element_blank element_text margin ggsave aes rel
#' @importFrom ggplot2 theme_linedraw theme element_text
#' @importFrom ggplot2 expansion scale_x_continuous aes_string
#' @importFrom ggseqlogo geom_logo
#'
#' @examples
#'
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' pwm <- seqArchR::make_PWMs(get_clBasVec_m(res,iter=1)[,1],
#'                         add_pseudo_counts = FALSE, sinuc = FALSE)
#'
#' viz_pwm(pwm_mat = pwm, method = "heatmap", fixed_coord = TRUE)
#'
#' viz_pwm(pwm_mat = pwm, method = "bits", fixed_coord = TRUE)
#'
viz_pwm <- function(pwm_mat, method = "heatmap", pos_lab = NULL,
                     pdf_name = NULL, fixed_coord = FALSE, bits_yax = "full"){

    pos_lab <- set_default_pos_lab2(pwm_mat, pos_lab)
    check_vars(pwm_mat, pos_lab)
    ##
    if(method == "heatmap"){
        ## Convert pwm_mat to df, heatmap by ggplot-way
        pwm_mat_df <- as.data.frame(pwm_mat)
        pwm_mat_df$Nucleotides <- rownames(pwm_mat_df)
        colnames(pwm_mat_df) <- c(pos_lab, "Nucleotides")
        pwm_mat_df_for_ggheatmap <- melt(pwm_mat_df, id.vars = c("Nucleotides"),
                                         variable.name = "positions")
        p1 <- get_ggheatmap(pwm_mat_df = pwm_mat_df_for_ggheatmap)
    }else{
        p1 <- get_ggseqlogo(pwm_mat, method = method, pos_lab = pos_lab)
        p1 <- add_lims_lab(p1, method, bits_yax)
    }
    ##
    p1 <- fix_coord(p1, nPos = length(pos_lab), method = method,
                    fixed_coord = fixed_coord)
    ##
    handle_pdffile_plotting(pdf_name, pl = p1, plw = 20, plh = 2.5)
    ##
    return(p1)

}
## =============================================================================

get_ggheatmap <- function(pwm_mat_df){
    ##
    p1 <- ggplot2::ggplot(data = pwm_mat_df, mapping = aes(
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
    ##
    p1
}
## =============================================================================


get_ggseqlogo <- function(pwm_mat, method, pos_lab){
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
    p1
}
## =============================================================================

handle_pdffile_plotting <- function(pdf_name, pl, plw, plh){
    if (!is.null(pdf_name)) {
        if (file.exists(pdf_name)) {
            warning("File exists, will overwrite", immediate. = TRUE)
        }
        ggplot2::ggsave(filename = pdf_name, plot = pl, device = "pdf",
                        width = plw, height = plh)
    }
}
## =============================================================================

set_default_pos_lab2 <- function(obj, pos_lab){
    if(!is.null(pos_lab)) return(pos_lab)
    if(!any(is.matrix(obj) || is(obj, "DNAString") || is(obj, "XString")))
        stop("Only handles matrix or DNAString object")
    if(is.matrix(obj)){
        pos_lab <- seq_len(ncol(obj))
        return(pos_lab)
    }
    pos_lab <- seq_len(Biostrings::nchar(obj))
}
## =============================================================================

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
## =============================================================================

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
## =============================================================================
