#' @title Plot cluster architectures as sequence logos.
#'
#' @description Given a collection of FASTA sequences as a DNAStringSet object,
#' and the clusters information, this function plots the architectures for all
#' clusters. If a name for the PDF file is provided, the resulting set of
#' architecture sequence logos are saved as a multi-page PDF.
#'
#' @param seqs Sequences as a \code{\link[Biostrings]{DNAStringSet}}.
#'
#' @param clust_list Clusters as a list of sequence IDs in each cluster.
#'
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the
#' positions are labeled from 1 to the length of the sequences.
#'
#' @param xt_freq Frequency of x-axis ticks.
#'
#' @param set_titles Specify TRUE if titles are to be written for the plots.
#' With FALSE, there are no titles for the plots. The title for each plot
#' includes the current cluster number, total number of clusters, start and
#' end sequence numbers in the collection.
#'
#' @param pdf_width,pdf_height Width and height in inches of the PDF file.
#' Default values are 11 and 2.
#'
#' @param pdf_name Specify the PDF filename.
#'
#' @param show Set TRUE if plot should be immediately shown/plotted. Default is
#' TRUE. By setting FALSE, one can simply collect the list of plots and use
#' any other approach to arrange/display them. See examples.
#'
#' @param ... Additional args passed to \code{\link{plot_ggseqlogo_of_seqs}}.
#'
#' @return A list of (ggplot2-based) sequence logo plots is returned. When a
#' valid file name is specified, the list of plots is also written to the PDF
#' file (one plot per page).
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom utils capture.output
#'
#' @examples
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' # Default position labels 1 to length of the sequences.
#' # Can also set pos_lab based on biology, e.g., use -50 to 49 denoting
#' # 50 basepairs upstream and 49 downstream of the transcription start site
#' # located at position 0.
#' arch_pl <- plot_arch_for_clusters(seqs = seqs_str(res),
#'                                   clust_list = res$clustSol$clusters,
#'                                   pos_lab = NULL,
#'                                   pdf_name = NULL,
#'                                   fixed_coord = TRUE)
#'
#'
#' # Using cowplot::plot_grid
#' arch_pl <- plot_arch_for_clusters(seqs = seqs_str(res),
#'                                   clust_list = res$clustSol$clusters,
#'                                   pos_lab = seq(100),
#'                                   method = "bits",
#'                                   pdf_name = NULL, show = FALSE)
#' cowplot::plot_grid(plotlist = arch_pl, ncol=1)
#'
#' # Plotting architecture sequence logos with probability instead of
#' # information content
#' arch_pl <- plot_arch_for_clusters(seqs = seqs_str(res),
#'                                   clust_list = res$clustSol$clusters,
#'                                   pos_lab = seq(100),
#'                                   method = "prob",
#'                                   pdf_name = NULL, show = FALSE)
#' cowplot::plot_grid(plotlist = arch_pl, ncol=1)
#'
#' @export
plot_arch_for_clusters <- function(seqs,
                                clust_list,
                                pos_lab = NULL,
                                xt_freq = 5,
                                set_titles = TRUE,
                                pdf_width = 11,
                                pdf_height = 2,
                                pdf_name = NULL,
                                show = FALSE,
                                ...){
    ##
    stopifnot(!is.null(clust_list))
    stopifnot(!is.null(seqs))
    ## seqs can be a DNAStringSet object or a character vector
    if(is(seqs, "character") && length(seqs) > 0){
        ## OK
    }else if(is(seqs, "DNAStringSet")){
        ## OK
    }else{
        stop("Expecting either a DNAStringSet object or a character vector")
    }
    ##
    if(is.null(pos_lab)){
        pos_lab <- seq_len(Biostrings::width(seqs[1]))
    }
    ##
    plot_titles <- make_plot_titles(clust_list, set_titles)
    plot_list <- lapply(seq_along(clust_list), function(x){
        pl <- plot_ggseqlogo_of_seqs(seqs=seqs[clust_list[[x]]],
                pos_lab = pos_lab,
                xt_freq = xt_freq,
                title = plot_titles[[x]],
                ...)
    })
    ##
    if(show){
        invisible(capture.output(plot_list))
    }
    ##
    if(!is.null(pdf_name)){
        grDevices::pdf(file=pdf_name, width=pdf_width, height=pdf_height)
        lapply(plot_list, print)
        dev.off()
        return(invisible(NULL))
    }
    return(plot_list)
}
## =============================================================================

make_plot_titles <- function(clust_list, set_titles){
    nClust <- length(clust_list)
    clust_lens <- unlist(lapply(clust_list, length))
    cumsums_clust_lens <- cumsum(clust_lens)
    clust_names <- sort(as.character(seq_along(clust_list)))
    ##
    clust_starts <- c(1, 1+cumsums_clust_lens[seq_len(nClust-1)])
    clust_ends <- cumsums_clust_lens
    if(set_titles){
        ##
        pl_titles <- lapply(seq_along(clust_starts), function(x){
            make_plot_title_str(x, nClust, clust_names[x], clust_lens[x],
                clust_starts[x], clust_ends[x])
        })
    }else{
        pl_titles <- lapply(seq_along(clust_starts), function(x) NULL)
    }
    pl_titles
}
## =============================================================================

make_plot_title_str <- function(i, n, name, this_size, st, ed){
    paste0("(", i , "/", n, ") Arch '",
        name, "': ", this_size, " sequences (",  st, "-",  ed, ")")
}
## =============================================================================

#' @title Plot sequence logo of a collection of sequences
#'
#' @description A wrapper to ggseqlogo plotting. Given a collection of
#' sequences, this function plots the sequence logo.
#'
#' @param seqs Collection of sequences as a
#' \code{\link[Biostrings]{DNAStringSet}} object.
#'
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the
#' positions are labeled from 1 to the length of the sequences.
#'
#' @param xt_freq Specify the frequency of the x-axis ticks.
#'
#' @param method Specify either 'bits' for information content or
#' 'prob' for probability.
#'
#' @param title The title for the plot. Deafult is NULL.
#'
#' @param bits_yax Specify 'full' if the information content y-axis limits
#' should be 0-2 or 'auto' for a suitable limit. The 'auto' setting adjusts
#' the y-axis limits according to the maximum information content of the
#' sequence logo. Default is 'full'.
#'
#' @param fixed_coord Specify TRUE if the aspect ratio of the plot should be
#' fixed, FALSE otherwise. Default is TRUE. When `method` argument is set to
#' 'bits', ratio is 4, when 'prob', ratio is 6.
#'
#' @return A sequence logo plot of the given DNA sequences.
#'
#' @seealso \code{\link{plot_arch_for_clusters}} for obtaining multiple
#' sequence logo plots as a list.
#'
#' @importFrom Biostrings width
#'
#' @examples
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' # Default, using information content on y-axis
#' pl <- plot_ggseqlogo_of_seqs(seqs = seqs_str(res, iter=1, cl=3),
#'                              pos_lab = seq_len(100), title = NULL,
#'                              fixed_coord = TRUE)
#' pl
#'
#' # Using probability instead of information content
#' pl <- plot_ggseqlogo_of_seqs(seqs = seqs_str(res, iter=1, cl=3),
#'                              pos_lab = seq_len(100), title = "",
#'                              method = "prob", fixed_coord = TRUE)
#' pl
#'
#' @export
plot_ggseqlogo_of_seqs <- function(seqs, pos_lab = NULL, xt_freq = 5,
                                    method = "bits", title = NULL,
                                    bits_yax = "full", fixed_coord = FALSE){
    ##
    if(is.null(pos_lab)){
        pos_lab <- seq_len(Biostrings::width(seqs[1]))
    }
    ##
    if(xt_freq > Biostrings::width(seqs[1])){
        xt_freq <- 5
    }
    ##
    # nPos <- length(pos_lab)
    # xtick_cal <- seq(0, nPos, by = xt_freq)
    # xtick_cal[1] <- 1
    # xtick_cal[length(xtick_cal)] <- nPos
    ##
    use_xtick_labs <- set_xtick_labels(pos_lab, xt_freq)
    ##
    foo_p <-
        ggseqlogo::ggseqlogo(
            as.character(seqs),
            seq_type = "dna",
            method = method
        ) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(axis.text.x = element_text(size = rel(0.9),
                                        angle = 90, hjust = 1, vjust=0.5),
                        axis.text.y = element_text(size = rel(0.9)),
                        panel.grid = element_blank()
        ) +
        ## Add additional bold tick labels
        ggplot2::scale_x_continuous(breaks = use_xtick_labs$breaks,
                                    labels = use_xtick_labs$labels,
                                    expand = expansion(mult = c(0, 0)))
    ##
    foo_p <- add_lims_lab(foo_p, method, bits_yax)
    foo_p <- fix_coord(foo_p, nPos=length(pos_lab), method, fixed_coord)
    ##
    if(!is.null(title)){
        foo_p <- foo_p + ggplot2::ggtitle(title)
        .msg_pstr("Plot title:", title, flg=TRUE)
    }
    ##

    foo_p
}
## =============================================================================

handle_only_pos_ticks <- function(xtick_cal, nPos){
    use_tick_pos <- c(1+xtick_cal[1], xtick_cal[-1])
    use_tick_pos <- use_tick_pos[use_tick_pos <= nPos]
}

handle_only_neg_ticks <- function(xtick_cal, nPos){
    use_tick_pos <- 1+xtick_cal
    use_tick_pos <- use_tick_pos[use_tick_pos <= nPos]
}


set_xtick_labels <- function(pos_lab, xt_freq){
    ## Note: This works in cases where the position labels are
    ## multiples of 5. Because one doesn't generally choose sequences of length,
    ## say, 923. One would select either 925/920.
    nPos <- length(pos_lab)
    xtick_cal <- seq(0, nPos, by = xt_freq)

    has_zero <- intersect(0, pos_lab)
    has_neg <- any(pos_lab < 0)
    has_pos <- any(pos_lab > 0)
    if(length(has_zero) > 0){
        ## has a zero
        use_tick_pos <- handle_only_neg_ticks(xtick_cal, nPos)
    }
    else if(length(has_zero) == 0 && !has_neg){
        ## has only negative values with no zero
        use_tick_pos <- handle_only_pos_ticks(xtick_cal, nPos)
    }
    else if(length(has_zero) == 0 && !has_pos){
        ## has only negative values with no zero
        use_tick_pos <- handle_only_neg_ticks(xtick_cal, nPos)
    }
    else{
        pos_idx <- which(pos_lab > 0)
        pos_lab_p <- pos_lab[pos_idx]
        pos_lab_n <- pos_lab[-pos_idx]
        ##
        nPos <- length(pos_lab_n)
        xtick_cal <- seq(0, nPos, by = xt_freq)
        utn <- handle_only_neg_ticks(xtick_cal, nPos)
        ##
        nPos <- length(pos_lab_p)
        xtick_cal <- seq(0, nPos, by = xt_freq)
        utp <- handle_only_pos_ticks(xtick_cal, nPos)
        ##
        use_tick_pos <- c(utn, length(pos_lab_n) + utp)
    }
    use_labs <- pos_lab[use_tick_pos]
    use_breaks <- use_tick_pos
    use_labs <- list(breaks = use_breaks, labels = use_labs)
}
## =============================================================================

set_ratio <- function(nPos, for100 = 4){
    ## ratio of 4 for 100 positions looks good.
    return(nPos*for100/100)
}
## =============================================================================


fix_coord <- function(p1, nPos, method, fixed_coord){
    if(fixed_coord){
        if(method == "bits") use_ratio <- set_ratio(nPos, 4)
        if((method == "prob" || method == "custom")){
            use_ratio <- set_ratio(nPos, 8)
        }
        if(method == "heatmap") use_ratio <- set_ratio(nPos, 2)
        p1 <- p1 + ggplot2::coord_fixed(ratio = use_ratio)
    }
    return(p1)
}
## =============================================================================
