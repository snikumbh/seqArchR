


#' @title Visualize the NMF basis vectors
#'
#' @description The given features matrix is visualized as a paired heatmap
#' and sequence logo where the positions are aligned for better
#' visualization., or as a single heatmap or as a single sequence logo.
#'
#' @param feat_mat The features matrix (basis vectors matrix) from seqArchR.
#' @param ptype Character vector of length one or two. Specify just one of
#' "heatmap" or "seqlogo" to visualize the basis vectors as such, or specify
#' a vector of length two for plotting both, heatmap and seqlogo. These are
#' then arranged one below the other, the first on top and the second under it.
#' @param method Specify either of "custom", "bits", or "probability" for
#' plotting sequence logo. Default is "bits".
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the
#' positions are labeled from 1 to the length of the sequences.
#' @param pdf_name Filename to save the plot, also provide the extension.
#' @param add_pseudo_counts Logical, taking values TRUE or FALSE, default
#' set to FALSE. Setting it to TRUE will enable adding pseudo-counts to the
#' features matrix.
#' @param sinuc_or_dinuc "sinuc" or "dinuc" for choosing between mono- and
#' dinucleotide profiles respectively.
#' @param fixed_coord Set this to TRUE to use a fixed aspect ratio for the
#' plot irrestive of the width and height of the PDF. Default is FALSE.
#'
#' @return nothing
#'
#' @export
#' @family visualization functions
#'
#' @importFrom ggplot2 theme margin
#'
#'
#' @examples
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' # Visualize basis vectors at iteration 1 of seqArchR result as heatmap and
#' # sequence logo
#' viz_bas_vec(feat_mat = get_clBasVec_m(res,iter=1), sinuc_or_dinuc = "dinuc",
#'                 ptype = c("heatmap", "seqlogo"))
#'
#'
#' # Visualize basis vectors at iteration 1 of seqArchR result as sequence logos
#' viz_bas_vec(feat_mat = get_clBasVec_m(res,iter=1), ptype = "seqlogo",
#'                  sinuc_or_dinuc = "dinuc")
#'
#'
#' # Visualizing basis vector for a single cluster as a heatmap
#' viz_bas_vec(feat_mat = as.matrix(get_clBasVec_m(res,iter=1)[,3]),
#'                  ptype = "heatmap", sinuc_or_dinuc = "dinuc")
#'
viz_bas_vec <- function(feat_mat, ptype = c("heatmap", "seqlogo"),
                        method = "bits", pos_lab = NULL, pdf_name = NULL,
                        add_pseudo_counts = FALSE, sinuc_or_dinuc = "sinuc",
                        fixed_coord = FALSE){
    ## Visualize all basis factors (expected as columns of the given features
    ## matrix) as heatmaps or seqlogos or both combined
    check_vars2(feat_mat)
    pos_lab <- set_default_pos_lab(feat_mat, sinuc_or_dinuc, pos_lab)
    expType <- c("seqlogo", "heatmap")
    if(!(length(match(ptype, expType)) == length(ptype))){
        stop("Expected values for arg 'ptype' are 'seqlogo' and 'heatmap'")
    }

    pl_list <- apply(feat_mat, MARGIN = 2, function(x) {
        set_sinuc <- TRUE
        if (sinuc_or_dinuc == "dinuc") {
            set_sinuc <- FALSE
        }
        pwm <- make_PWMs(x, add_pseudo_counts = FALSE, sinuc = set_sinuc)
        pl <- vector("list", length(expType))
        ##
        if("heatmap" %in% ptype){
            p1 <- viz_pwm(pwm_mat = pwm, method = "heatmap", pos_lab = pos_lab)
            # p1 <- p1 + theme(plot.margin = margin(0,0,0,0, unit="cm"))
            if(length(ptype) == 1) return(p1)
            pl[[match("heatmap", ptype)]] <- p1
        }
        if("seqlogo" %in% ptype){
            p2 <- viz_pwm(pwm_mat = pwm, method = method, pdf_name = NULL,
                            pos_lab = pos_lab, fixed_coord = fixed_coord)
            # p2 <- p2 + theme(plot.margin = margin(0,0,0,0, unit="cm"))
            if(length(ptype) == 1) return(p2)
            pl[[match("seqlogo", ptype)]] <- p2
        }
        if(length(ptype) == length(expType)){
            check_cowplot()
            pl <- cowplot::plot_grid(plotlist = pl, nrow = 2, align="v")
            return(pl)
        }
        pl
    })

    handle_pdffile_plotlist(pdf_name, pl_list)
    pl_list
}
## =============================================================================

handle_pdffile_plotlist <- function(pdf_name, pl_list, pw = 20, ph = 4){
    if(!is.null(pdf_name)){
        if (file.exists(pdf_name)) {
            warning("File exists, will overwrite", immediate. = TRUE)
        }
        grDevices::pdf(file=pdf_name, width=pw, height=ph)
        lapply(pl_list, print)
        dev.off()
        return(invisible(NULL))
    }
}
## =============================================================================

set_default_pos_lab <- function(feat_mat, sinuc_or_dinuc, pos_lab){
    if(!is.null(pos_lab)) return(pos_lab)
    deno <- 4
    if(sinuc_or_dinuc == "dinuc"){
        deno <- 16
    }
    pos_lab <- seq_len(nrow(feat_mat)/deno)
}
## =============================================================================

check_cowplot <- function(){
    if(!requireNamespace("cowplot", quietly = TRUE)){
        stop("Please install the R package 'cowplot' to use this function")
    }
}
## =============================================================================

check_vars2 <- function(feat_mat){
    if (!is.matrix(feat_mat)) {
        stop("feat_mat not of type matrix")
    }
    if (sum(dim(feat_mat)) == 2 && is.na(feat_mat)) {
        stop("Empty feat_mat")
    }
}
## =============================================================================
