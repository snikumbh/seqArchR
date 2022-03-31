## Getter functions ============================================================

## get functions for components of seqArchR's result object
## These keep such accesses in other functions agnostic to how they are
## maintained in the seqArchR result object

#' @title Get functions for seqArchR result object
#'
#' @description Basis vectors' information for the selected iteration.
#'
#' @param res seqArchR result object.
#' @param iter Choose the iteration of seqArchR result to get from.
#'
#' @return get_clBasVec A list with two elements `nBasisVectors` (integer) and
#' `basisVectors` (matrix).
#'
#' @examples
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' k <- get_clBasVec_k(res=res, iter=2)
#'
#' bMat <- get_clBasVec_m(res=res, iter=2)
#'
#' ## cluster labels of sequences from final clustering
#' scLab <- get_seqClLab(res=res, iter=2)
#'
#' @export
get_clBasVec <- function(res, iter){
    return(res$clustBasisVectors[[iter]])
}


#' @describeIn get_clBasVec Get the number of basis vectors (clusters)
#' at the  selected iteration.
#'
#' @return get_clBasVec_k The number of basis vectors (integer).
#'
#' @export
get_clBasVec_k <- function(res, iter){
    return(res$clustBasisVectors[[iter]]$nBasisVectors)
}

#' @describeIn get_clBasVec The basis vectors matrix at the selected
#' iteration. Note that eatures along rows.
#'
#' @return get_clBasVec_m The basis vectors' matrix with features along the
#' rows of the matrix.
#'
#' @export
get_clBasVec_m <- function(res, iter){
    return(res$clustBasisVectors[[iter]]$basisVectors)
}

#' @describeIn get_clBasVec Get the cluster IDs for each sequence. Note that
#' order of sequences here is as per the input.
#'
#' @return get_seqClLab A character vector denoting the cluster IDs for
#' each sequence.
#'
#' @seealso \code{\link{seqs_str}}
#' @export
get_seqClLab <- function(res, iter = NULL){
    if(is.null(iter)){
        return(res$clustSol$seqsClustLabels)
    }else{
        return(res$seqsClustLabels[[iter]])
    }
}



#' @title Get sequences from the seqArchR result object
#' @description Wrapper to fetch sequences from the seqArchR result object as
#' character
#'
#' @param res seqArchR result object
#'
#' @param iter Specify the iteration of seqArchR result. If set to NULL
#' (the default), the original set of sequences (`seqArchRresult$rawSeqs`) is
#' returned.
#'
#' @param cl Specify the cluster number. Sequences belonging to this cluster in
#' iteration `iter` of seqArchR result are returned as character. When `iter` is
#' NULL, this is treated as denoting the cluster number in seqArchR's final
#' clustering solution (`seqArchRresult$clustSol$clusters`).
#'
#' @param ord Specify TRUE if sequences are ordered by clusters. The original
#' ordering of the sequences can be fetched by setting `iter` to NULL and `ord`
#' to FALSE.
#'
#' @details Setting iter to NULL will fetch sequences as per the final
#' clustering solution of seqArchR (`clustSol$clusters`). When `iter` is not
#' NULL, use `cl` to further choose a particular cluster. When `cl` is NULL,
#' the set of sequences returned can be ordered by clusters with `ord = TRUE`.
#' Using `ord = FALSE` fetches the sequences by their original order.
#'
#' @return The selected DNA sequences from the DNAStringSet object as a
#' character vector.
#'
#' @export
#'
#' @examples
#'
#' res <- system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE)
#'
#' # Fetch sequences from 2nd cluster of seqArchR's final solution
#' ans <- seqArchR::seqs_str(readRDS(res), iter=NULL, cl=2)
#'
#' # Fetch all sequences ordered by the final clustering
#' ans <- seqArchR::seqs_str(readRDS(res), iter=NULL, cl=NULL, ord=TRUE)
#'
#' # Fetch sequences belonging to first cluster in seqArchR's first iteration
#' ans <- seqArchR::seqs_str(readRDS(res), iter=1, cl=1)
#'
#'
seqs_str <- function(res, iter = NULL, cl = NULL, ord = FALSE){

    ## when iter is NULL, return sequences belonging to a cluster (specified
    ## by cl) in the final solution
    if(is.null(iter) && is.null(cl) && !ord){
        return(as.character(res$rawSeqs))
    }
    if(is.null(iter) && is.null(cl) && ord){
        use_ord <- unlist(res$clustSol$clusters)
        return(as.character(res$rawSeqs[use_ord]))
    }
    if(is.null(iter) && !is.null(cl)){
        cl_mem <- res$clustSol$clusters[[cl]]
        return(as.character(res$rawSeqs[cl_mem]))
    }
    if(!is.null(iter) && !is.null(cl)){
        clust_list <- get_seqs_clust_list(get_seqClLab(res, iter))
        cl_mem <- clust_list[[cl]]
        return(as.character(res$rawSeqs[cl_mem]))
    }
    if(!is.null(iter) && is.null(cl) && ord){
        clust_list <- get_seqs_clust_list(get_seqClLab(res, iter))
        use_ord <- unlist(clust_list)
        return(as.character(res$rawSeqs[use_ord]))
    }

}


## =============================================================================

# @title Assign samples to clusters
#
# @param samplesMatrix
#
# @return A list that can be assigned as an element in globClustAssignments
#
.get_cluster_memberships_per_run <- function(samplesMatrix, iChunksColl = NA,
                                        iChunkIdx = NA, test_itr = NA,
                                        oChunkIdx = NA) {
    .assert_seqArchR_samplesMatrix(samplesMatrix)
    nClusters <- nrow(samplesMatrix)
    returnClusterAsList <- vector("list", nClusters)
    ## Fetch the cluster memberships for each sample (along the columns)
    clustMemberships <- apply(samplesMatrix, 2, which.max)
    ##
    if (length(unique(clustMemberships)) != nClusters) {
        cli::cli_alert_warning(
            c("Basis vector that got no sequences assigned: ",
            setdiff( unique(clustMemberships), seq_len(nClusters))))
    }

    return(clustMemberships)
}
## =============================================================================

.adjustSampleMemberships <- function(old_mem, samplesMatrix, has_overfit){

    ## Assign all overfitted samples back to cluster 1
    new_mem <- old_mem
    ##
    unoverfit <- setdiff(unique(old_mem), has_overfit)
    toassign <- utils::tail(sort(unoverfit, decreasing = FALSE), 1)
    for(i in has_overfit){
        ## simply assigning 1 is problematic if/when 1 is among the overfit
        ## and/or the only one overfit. This has earlier resulted in an error.
        ## Instead assign the lowest cluster ID which is not overfit
        new_mem[which(old_mem == i)] <- toassign
    }
    names(new_mem) <- NULL

    ## At this point new memberships for samples can contain clust ids
    ## such 1, 4 and 5 because 2 and 3 are overfit. (i.e. they are not
    ## consecutive). We need to get them to be consecutive
    new_mem_elem <- sort(unique(new_mem))
    for(i in seq_along(new_mem_elem)){
        new_mem[new_mem == new_mem_elem[i]] <- i
    }
    ##

    new_mem_labs <- unique(new_mem)

    stopifnot((length(unique(old_mem)) - length(has_overfit))
                == length(new_mem_labs))
    if(length(new_mem_labs) == 1) stopifnot(unique(new_mem_labs) == 1)
    ##
    new_mem
}
## =============================================================================

## This function detects if, for a given number of factors, overfitting has
## occurred.
##
## Return value: the outlier clust ids or numeric(0) if no outlier
##
## Input args:
## -- samplesMat
## -- clusterMemberships
## -- minSeqs (This can be a different value than that passed via
## seqArchR config)
## --
.detect_overfitting <- function(samplesMat, clustMemberships, minSeqs,
                                sizeThreshold = 10){
    # put factors along columns
    samplesMat <- t(samplesMat)
    len_clusts <- unlist(lapply(seq_len(ncol(samplesMat)), function(x){
        length(which(clustMemberships == x))
        }))

    out_cl_idx_size <- which(len_clusts <= sizeThreshold)

    qual_cl_idx <- which(len_clusts < minSeqs & len_clusts > sizeThreshold)

    ## Update:
    ## -- Keep only IQR overfitted outliers when the qualified clusters are
    ##    among the ones in between.
    ## -- For any in-between outliers, assign them back to the first cluster,
    ##    instead ofthe next maximum loading cluster
    ## -- If the m qualifying clusters are from n-m to n, they can be
    ##    considered for either IQR or range-outlier case
    ## --
    out_cl_idx_other <- NULL
    out_cl_idx <- out_cl_idx_size

    if(length(qual_cl_idx) > 0){
        out_cl_idx_other <- .compare_clusters(samplesMat, clustMemberships,
                                        qual_cl_idx, zscore_thresh = 5)
        out_cl_idx <- unique(union(out_cl_idx_other, out_cl_idx_size))
    }

    ## Account for when all available clusters qualify to be checked for
    ## overfitting cases
    if(length(out_cl_idx) == length(len_clusts)){
        out_cl_idx <- setdiff(out_cl_idx, 1)
    }

    if(length(out_cl_idx) > 0) return(out_cl_idx)
    ##

    return(intersect(1,2))

}
## =============================================================================

.compare_clusters <- function(samplesMat, clustMemberships, qual_cl_idx,
                                zscore_thresh = 5){

    clustwise_matlist <- lapply(seq_len(ncol(samplesMat)), function(x){
        temp_mat <- matrix(samplesMat[clustMemberships == x,],
            ncol = ncol(samplesMat), byrow = FALSE)
        temp_mat
    })

    # Adjust qualifying IDs based on Update in compare_clusters func
    # -- Allow any ID for IQR check

    out_clust_iqr <- .compare_iqr(clustwise_matlist, qual_cl_idx,
                                    zscore_thresh = zscore_thresh)
    # out_clust_range <- NULL
    out_clust_range <- .compare_range(clustwise_matlist, qual_cl_idx,
                                    zscore_thresh = zscore_thresh)

    ##
    ## Those that qualify as outliers by IQR comparison, will also
    ## qualify by range comparison. But for those which do not qualify
    ## by IQR but by range, we need further checks for them.
    ## -- where does their box lie,
    ##      -- near the maximum? Not an outlier
    ##      -- near others? outlier
    ##
    check_these <- setdiff(out_clust_range, out_clust_iqr)
    if(length(check_these) > 0){
        out_clust_range <- .add_check_range_outliers(check_these,
            clustwise_matlist)
    }
    ##
    out_clust_size <- NULL
    ## This is commented and out_clust_size is set NULL because this is
    ## already handled above. See detect_overfitting/sizeThreshold
    ##
    return(unique(c(out_clust_iqr, out_clust_range, out_clust_size)))

}
## =============================================================================

.add_check_range_outliers <- function(check_these, clustwise_matlist){
    ## get boxplot.stats
    all_stats_max <- lapply(clustwise_matlist, function(x){
        tempx <- lapply(seq_len(ncol(x)), function(y){
            max(as.vector(x[,y]))
        })
        unlist(tempx)
    })

    all_stats_midpoints <- lapply(clustwise_matlist, function(x){
        tempx <- lapply(seq_len(ncol(x)), function(y){
            tempy <- grDevices::boxplot.stats(as.vector(x[,y]))
            sum(tempy$stats[c(2,4)])/2
        })
        unlist(tempx)
    })

    #####
    all_stats_midpoint_nl <- unlist(all_stats_midpoints)
    #####
    out_cl_idx <- NULL
    for(i in check_these){
        ## use midpoints to compare locations of boxes
        check_midpoint <- all_stats_midpoints[[i]][i]
        distWithMax <- abs(check_midpoint - all_stats_max[[i]][i])
        distWithMedianMidpoints <- check_midpoint -
            stats::median(all_stats_midpoint_nl)

        if(distWithMax > distWithMedianMidpoints){
            ## outlier
            out_cl_idx <- c(out_cl_idx, i)
        }else{
            # do nothing
        }
    }
    return(out_cl_idx)
}
## =============================================================================

## Discard clusters with length smaller than 10 sequences
## This func is called when IQR and range fail to identify any outliers
## In which case if the IQRs are all zeros (or very small), that is sure shot
##  overfitting occurrence
.compare_size <- function(clustwise_matlist, qual_cl_idx, threshold = 0.01){

    all_iqr <- lapply(clustwise_matlist, function(x){
        ## get column-wise IQRs
        matrixStats::colIQRs(x)
    })

    ## detect if all IQRs for a cluster are 0
    if_all_zeros <- unlist(lapply(all_iqr, function(x){
        ifelse(all(!(x > threshold)), TRUE, FALSE)
        ## enter here if all zeros is detected
    }))


    ## Alternatively, simply check the number of samples in this cluster,
    ## if < 10, discard.
    ## -- This is already handled in a function hierarchical above this func.
    ##    See detect_overfitting/sizeThreshold

    out_cl_idx <- intersect(which(if_all_zeros), qual_cl_idx)
    return(out_cl_idx)
}
## =============================================================================

## return indices of outliers by comparing IQRs
.compare_iqr <- function(clustwise_matlist, qual_cl_idx, zscore_thresh = 5){

    ncl <- ncol(clustwise_matlist[[1]])

    all_iqr <- lapply(clustwise_matlist, function(x){
        ## get column-wise IQRs
        matrixStats::colIQRs(x)
    })
    all_iqr_vec <- unlist(all_iqr)
    ##
    all_iqr_mad <- stats::mad(all_iqr_vec)
    iqr_zscore <- (all_iqr_vec - stats::median(all_iqr_vec))/all_iqr_mad
    ##
    out_idx <- which(iqr_zscore > zscore_thresh)
    if(length(out_idx) > 0){
        clust_id <- ceiling(out_idx / ncl)
        return(intersect(clust_id, qual_cl_idx))
    }else{
        # do nothing
    }
    return(NULL)
}
## =============================================================================

## return indices of outliers by comparing ranges
.compare_range <- function(clustwise_matlist, qual_cl_idx, zscore_thresh = 5){

    ncl <- ncol(clustwise_matlist[[1]])
    all_range <- lapply(clustwise_matlist, function(x){
        apply(apply(as.matrix(x), MARGIN = 2, range), MARGIN = 2, diff)
    })
    all_range_vec <- unlist(all_range)
    ##
    all_range_mad <- stats::mad(all_range_vec)
    range_zscore <- (all_range_vec - stats::median(all_range_vec))/all_range_mad
    ##
    out_idx <- which(range_zscore > zscore_thresh)
    if(length(out_idx) > 0){
        clust_id <- ceiling(out_idx / ncl)
        return(intersect(clust_id, qual_cl_idx))
    }else{
        # do nothing
    }
    return(NULL)
}
## =============================================================================

.assign_samples_to_clusters <- function(clusterMembershipsVec, nClusters,
                                        iChunksColl, iChunkIdx) {
    returnClusterAsList <- vector("list", nClusters)
    for (i in seq_along(returnClusterAsList)) {
        returnClusterAsList[[i]] <-
            iChunksColl[[iChunkIdx]][clusterMembershipsVec == i]
    }
    return(returnClusterAsList)
}
## =============================================================================


#' @title Collate sequence IDs from an existing clustering according to a new,
#' given clustering of the existing clusters
#'
#' @description Collate sequences original divided across n clusters into a
#'  new set of m clusters. These m clusters obtained by clustering the original
#'  n clusters.
#' Assume a collection of 100 sequences across seven existing
#' clusters. These seven clusters are collated to obtain three resulting
#' clusters. Collating 100 sequences distributed across the seven clusters into
#' the resulting three clusters can be achieved with collate clusters
#'
#' @param to_clust A list giving clustering of factors. In other words this is
#' the clustering of clusters of sequences given in orig_clust
#' @param orig_clust  A list of sequence IDs in existing clusters
#'
#' @return A list with sequence IDs collated by the specified clustering
#' @export
#'
#' @examples
#'
#'
#'
#' set.seed(123)
#' n <- 7; nn <- 100
#' orig_clust_labels <- ceiling(n * runif(nn))
#' orig_clust <- seqArchR::get_seqs_clust_list(orig_clust_labels)
#'
#' to_clust <- list(c(1,4), c(2,3,5), c(6,7))
#'
#' collate_clusters(to_clust = to_clust, orig_clust = orig_clust)
#'
collate_clusters <- function(to_clust, orig_clust) {
    if (is.null(to_clust)) {
        ## When/if factor to_clust is not performed, the
        ## orig_clust variable is directly assigned
        return(orig_clust)
    }else{
        .assert_seqArchR_globClustAssignments(orig_clust)
        ##
        nClusters <- length(to_clust)
        clustSizes <- unlist(lapply(to_clust, length))
        ##
        collClAssign <- vector("list", nClusters)
        for(clustIdx in seq_along(collClAssign)){
            temp <- unlist(lapply(to_clust[[clustIdx]], function(x){
                orig_clust[[x]]
            }))
            collClAssign[[clustIdx]] <- temp
        }
        ##
        return(collClAssign)
    }
}
## =============================================================================

.msg_print <- function(vec_or_list){
    if(is.list(vec_or_list)){
        return(paste(vec_or_list, collapse = "\n"))
    }
    if(is.vector(vec_or_list)){
        return(paste(vec_or_list, collapse = ", "))
    }
}
## =============================================================================

.msg_pstr <- function(..., flg){
    if(flg){
        message(...)
    }
}
## =============================================================================

## @title Update labels of sequences in a cluster
## @param oldSeqsClustLabels
##
## @param collatedClustAssignments
## @param flags List. Flags variable from seqArchR config
##
## @return newSeqsClustLabels Vector os clustLabels
.update_cluster_labels <- function(oldSeqsClustLabels,
                                collatedClustAssignments) {
    ## Important: The collatedClustAssignments variable is a list where each
    ## element holds the indices of sequences falling in the same cluster,
    ## like so: For a set of 200 sequences, collatedClustAssignments is
    ## [[1]]
    ## [1]   1   2   6  11  20  23  32  46  50  52  63  68  71  72  75
    ## 80  82  83  85  86  87  88  92  99 102 104 105
    ##
    ## [[2]]
    ## [1]   3   5   7   8  10  12  13  16  17  18  19  22  25  30  31
    ## 36  37  38  39  41  43  45  47  49  51  53  54
    ## --
    ## This is different than the globClustAssignment variable which will hold
    ## the indices of the factors that are combined into one cluster.
    ## For instance, when there are 25 factors/clusters which have combined
    ## to give 5 clusters, globClustAssignments will hold something like this:
    ## [[1]]
    ## [1] "3"  "13" "11" "16" "6"
    ##
    ## [[2]]
    ## [1] "15" "1"  "8"  "25" "18"
    ##
    .assert_seqArchR_seqsClustLabels(oldSeqsClustLabels)
    .assert_seqArchR_globClustAssignments(collatedClustAssignments)
    nClusters <- length(collatedClustAssignments)
    ## Use numerics w/ as.character and pre-sort to have them in the
    ## order that will be returned by levels (in get_seqs_clusters_in_list fn)
    candidateClustLabels <- sort(as.character(seq_len(nClusters)))
    newSeqsClustLabels <- oldSeqsClustLabels
    for (i in seq_len(nClusters)) {
        needUpdateIdx <- collatedClustAssignments[[i]]
        newSeqsClustLabels[needUpdateIdx] <-
            vapply(newSeqsClustLabels[needUpdateIdx], function(x) {
                        paste0(candidateClustLabels[i])}, character(1)
                )
    }
    .assert_seqArchR_seqsClustLabels(newSeqsClustLabels)
    return(newSeqsClustLabels)
}
## =============================================================================

# @title Prepare chunks out of given set of sequences (sequence IDs)
# @param total_set set of sequences to be chunked
# @param reqdChunkSize the given chunk size
# @param flags Specify the flags from the config param set for seqArchR
#
# @return A list with chunks of sequences (sequence IDs) as its elements
.prepare_chunks <- function(total_set, reqdChunkSize) {
    ## total_set is the set of seq_ids to be chunked (not array indices)
    if (is.null(total_set)) {
        stop("Preparing chunks, 'total_set' is NULL")
    }
    chunkLength <- length(total_set)
    if (chunkLength == 0) {
        stop("Preparing chunks, length of 'total_set' is 0")
    }
    .assert_seqArchR_chunkSize_independent(reqdChunkSize)
    ##
    ## When chunkLength (i.e., total sequences) < reqdChunkSize, the else
    ## condition, return the totalSet as the only chunk.
    if (chunkLength > reqdChunkSize) {
        chunkStarts <- seq(1, chunkLength, by = reqdChunkSize)
        chunkEnds <- seq(reqdChunkSize, chunkLength, by = reqdChunkSize)
        if (length(chunkStarts) > length(chunkEnds)) {
            if ((chunkLength - chunkStarts[length(chunkStarts)]) >
                round(0.5*reqdChunkSize)) {
                chunkEnds <- append(chunkEnds, chunkLength)
            } else {
                chunkStarts <- chunkStarts[-length(chunkStarts)]
                chunkEnds[length(chunkEnds)] <- chunkLength
            }
        }
        ##
        preparedChunks <- vector("list", length(chunkStarts))
        for (i in seq_along(chunkStarts)) {
            preparedChunks[[i]] <-
                total_set[seq(chunkStarts[i],chunkEnds[i],by=1)]
        }
        ##
    } else {
        preparedChunks <- vector("list", 1)
        preparedChunks[[1]] <- total_set
    }
    if (length(preparedChunks) < 1) {
        stop("Preparing chunks, length was 0")
    }
    return(preparedChunks)

}
## =============================================================================

# @title Handle clustering of NMF factors
#
# @param globFactorsMat Specify the NMF factors as a matrix with the individual
# factors along the columns.
#
# @param distMethod Specify the distance measure to be used. Default is cosine
# measure.
#
# @param clustMethod Specify hc. This is the only option available.
#
# @param linkage Specify one of linkage options for hclust.
#
# @param flags Specify the flags from the config param set for seqArchR.
#
# @param returnOrder
#
# @param useCutree By setting this arg to TRUE, the cutree version of
# get_clusters func is used. This func uses silhouette and/or
# Calinski-Harabasz index. See additional details there. Note: currently, this
# is the only avaiable func so useCutree should beset to TRUE always, FALSE
# will not work.
#
# @param minClusters Passed to get_clusters function. See explanation there.
#
# @param parentChunks
#
# @return List with clustered factors
#
# Change on 2020-12-28:
# - change default value for distMethod to 'euclid',
# - change default value for linkage to 'ward.D',
# - remove returnOrder arg. not needed; this was useful during development
# -
# Change on 2021-06-03:
# - Added back returnOrder argument
#
.handle_clustering_of_factors <- function(globFactorsMat,
                                        clustMethod = "hc",
                                        linkage = "average",
                                        distMethod = "cosangle",
                                        flags = list(debugFlag = FALSE,
                                                    verboseFlag = FALSE,
                                                    plotVerboseFlag = FALSE,
                                                    timeFlag = FALSE),
                                        useCutree = TRUE,
                                        minClusters = 2,
                                        parentChunks = NULL,
                                        returnOrder = FALSE,
                                        ...) {
    ##
    .assert_seqArchR_featuresMatrix(globFactorsMat)
    .assert_seqArchR_flags(flags)
    dbg <- flags$debugFlag
    vrbs <- flags$verboseFlag
    ##
    gfDisMat <- .compute_factor_distances(globFactorsMat, distMethod)
    if(clustMethod == "hc"){
        if(distMethod == "euclid") use_threshold <- 3
        if(distMethod == "cor") use_threshold <- 0.75

        as_dist_mat <- stats::as.dist(gfDisMat)
        temp_hclust <- stats::hclust(as_dist_mat, method = linkage)
        ord <- temp_hclust$order
        .msg_pstr("Cluster order by ", linkage," linkage w/ ",
            distMethod, " distance:", flg=dbg)
        .msg_pstr("New order: ", .msg_print(ord), flg=dbg)
        for(k in seq_len(length(ord)-1)){
            .msg_pstr(paste(paste(ord[k], ord[k+1], sep=", "), "=",
                    gfDisMat[ord[k], ord[k+1]], collapse="\n"), flg=dbg)
        }
        .msg_pstr("New order: ", .msg_print(ord), flg=dbg)
        if(returnOrder) return(temp_hclust)
        .msg_pstr("=== Fetching clusters === ", flg=dbg)
        if(useCutree){
            clustList <- .get_clusters_from_hc_using_cutree(
                            hcObj = temp_hclust, distMat = as_dist_mat,
                            hStep = 0.05, parentChunks = parentChunks,
                            minClusters = minClusters,
                            distThreshold = use_threshold,
                            verbose = flags$debugFlag, ...)
        }

        if(length(clustList) == 1){
            clustList <- lapply(seq_len(max(clustList[[1]])),
                                function(x){x})
        }
        .msg_pstr("== DONE ==", "\nClustList: ", .msg_print(clustList), flg=dbg)
        return(clustList)
    }
}
## =============================================================================

.unfurl_nodeList <- function(nodeList, vec_ver = FALSE, vrbs=FALSE){
    ##
    returnVal <- .assert_seqArchR_list_properties(nodeList)
    if(returnVal != "FOO") stop(returnVal)
    ##
    if(vec_ver){
        element_lengths <- unlist(lapply(nodeList, length))
        if(any(element_lengths != 1)){
            .msg_pstr("Needs unfurling...", flg=vrbs)
            new_list <- vector("list", sum(element_lengths))
            iter1_new <- 0
            iter2_old <- 1
            while(iter2_old <= length(nodeList)){
                if(length(nodeList[[iter2_old]]) == 1){
                    iter1_new <- iter1_new + 1
                    new_list[[iter1_new]] <- nodeList[[iter2_old]]
                }else{
                    for(i in seq_along(nodeList[[iter2_old]])){
                        iter1_new <- iter1_new + 1
                        new_list[[iter1_new]] <- nodeList[[iter2_old]][[i]]
                    }
                }
                iter2_old <- iter2_old + 1
            }
            return(new_list)

        }else{
            .msg_pstr("No unfurling...", flg=vrbs)
            return(nodeList)
        }
    }else{
        element_lengths <- unlist(lapply(nodeList, function(elem){
            ifelse(is.list(elem), length(elem), 1)
            # else 1
        }))
        if(any(element_lengths != 1)){
            .msg_pstr("Needs unfurling...", flg=vrbs)
            new_list <- vector("list", sum(element_lengths))
            iter1 <- 0
            iter2 <- 1
            while(iter2 <= length(nodeList)){
                if(!is.list(nodeList[[iter2]])){
                    iter1 <- iter1 + 1
                    new_list[[iter1]] <- nodeList[[iter2]]
                }else{
                    for(i in seq_along(nodeList[[iter2]])){
                        iter1 <- iter1 + 1
                        new_list[[iter1]] <- nodeList[[iter2]][[i]]
                    }
                }
                iter2 <- iter2 + 1
            }
            return(new_list)
        }else{
            .msg_pstr("No unfurling...", flg=vrbs)
            return(nodeList)
        }
    }

}
## =============================================================================


# @title Order cutree result by hierarchical clustering order
#
# @description Make clusters obtained by cutree to follow the
# hierarchical clustering order
#
# @param cutree_res Result object of cutree
# @param hcorder The ordering from hierarchical clustering that is accessible
# in its result object
#
# @return A list of clusters, including each clusters' elements, ordered by
# their appearance in the hierarchical clustering order
#
# @examples
#
# @export
fetch_cutree_by_hc_order <- function(clust_list, cutree_res = NULL, hcorder){

    if(!is.null(cutree_res)){
        names(cutree_res) <- NULL
        clust_list <- lapply(unique(cutree_res),
                                function(x){which(cutree_res == x)})
    }

    ## Arrange clusters themselves by the hc_order
    starters <- lapply(clust_list, function(x){x[[1]]})
    starters_pos <- unlist(lapply(unlist(starters), function(x){
        which(x == hcorder)
    }))

    use_idx <- sort(starters_pos, index.return = TRUE)

    clust_list <- lapply(seq_along(use_idx$ix), function(x){
        clust_list[[use_idx$ix[x]]]
    })

    # For each ID in list elements, get their position in the hcorder vector
    cl_order_idx <- unlist(lapply(unlist(clust_list), function(x){
        which(x == hcorder)
    }))

    ## Arrange each cluster elements by hc_order

    cl_lens <- unlist(lapply(clust_list, length))

    cl_ends <- cumsum(cl_lens)
    cl_starts <- c(1, 1+cl_ends[seq(length(cl_ends)-1)])

    stopifnot(tail(cl_ends, 1) == length(hcorder))

    clust_list_new <- vector("list", length(cl_lens))
    for(x in seq_along(cl_starts)){
        adjust_idx <- cl_order_idx[cl_starts[x]:cl_ends[x]]
        if(x > 1){
            adjust_idx <- cl_order_idx[ cl_starts[x]:cl_ends[x] ] - cl_ends[x-1]
        }
        clust_list_new[[x]][ adjust_idx ] <- clust_list[[x]]
    }

    stopifnot(length(unlist(clust_list_new)) == length(hcorder))
    stopifnot(length(unique(unlist(clust_list_new))) == length(hcorder))

    return_vec <- rep(0, length(hcorder))

    return_vec <- lapply(seq_along(clust_list_new), function(x){
        return_vec
    })

    return(clust_list_new)
}



############################## CUTREE VERSION ##################################
## This works best/is best referred with ward.D linkage
# @title Get clusters out of hierarchical clustering object using cutree.
# Used internally by seqArchR.
#
# @description Clusters from a hierarchical clustering object are obtained
# by using cutree at different heights of the tree. The optimum number of
# clusters are decided based on the average silhouette value of the clustering.
#
# @param hcObj The hierarchical clustering object as returned by
# \code{\link[stats]{hclust}}.
# @param distMat The distance matrix that was used by hclust, a
# \code{\link[stats]{dist}} object.
# @param hStep Numeric. The step size used to increment height values for
# cutree. Default value is 0.05.
# @param parentChunks List. Specify the factor numbers in the previous
# iteration of seqArchR that factors in the current iteration resulted from.
# Default value is NULL.
# @param keepSiblingsUncollated Logical. Specify TRUE if all clusters from a
# single parent chunk should not be collated, even when recommended so based
# on silhouette value computation. Default value FALSE. i.e. if they get
# collated, it is left as collated.
# @param enableSwitchSilToCH Logical. When and if the clustering result has
# singletons, using the average silhouette value may not be reliable to decide
# on the optimum number of clusters. In this case, if this argument is TRUE,
# the Calinski-Harabasz index is used instead. Setting it to FALSE, the
# switch is disabled, i.e., average silhouette value is used. Default is FALSE.
# @param minClusters Integer. Specify the minimum number of clusters to be
# fetched from the HAC solution.
# @param distThreshold Numeric. Specify a threshold of units of distance for
# considering any two elements as close enough to be merged in the same
# cluster. The default value is specified for Euclidean distance.
# @param verbose Logical. Specify TRUE for verbose output.
#
# @importFrom stats cutree
# @importFrom fpc calinhara
# @importFrom cluster silhouette
#
## Addtional comments:
## Keeping the min to 0 (as below) leads to as many clusters as
## the number of elements. This errors when at mean[, "sil_width"]
## cut_heights <- seq(0, max(hcObj$height), by = hStep)
##
## Important: Address the question of whether any merging is required?
## -- We could use information from seqArchR iteration
##    -- Obtain the cutree-using-h clustering here.
##    -- Is it combining consecutive factors, like 1-2-3, 3-4 etc.?
##          a. This could mean that it is undoing the clusters
##          identified by seqArchR in a previous iteration.
##          b. But, if 1 came from chunk 1 and 2-3 came from chunk 2,
##          and if the clustering results in combining 1-2, then this
##          is a plausible merge and should be left alone.
##          c. So, we use a parentChunks variable that notes the parent
##          chunk for each of the current factors. This info can be
##          used to make this decision unambiguously.
##
## Update: 2020-12-13: We currently let these get combined if
## HAC+cutree deems it fine, and we recommend that the clusters at the
## iteration of seqArchR can be left uncollated for reference.
## The final stage will then hold the best possible set of clusters
## guided by silhouette value. Now handled by logical argument
## keepSiblingsUncollated.
##
## Update: 2020-12-24: With singletons in a clustering, using
## silhouette values could be problematic, because it arbitrarily
## assigns s(i) = 0 for a singleton cluster. This can lead to
## smaller average silhouette value, artificially making that
## clustering look bad. Thus, if and when we detect such a case,
## test using the Calinski-Harabasz (CH) index automatically.
## See argument `enableSwitchSilToCH`.
.get_clusters_from_hc_using_cutree <- function(hcObj, distMat, hStep = 0.05,
                                            parentChunks = NULL,
                                            keepSiblingsUncollated = FALSE,
                                            enableSwitchSilToCH = FALSE,
                                            minClusters = 2,
                                            ## number of clusters in a
                                            ## previous iteration on
                                            ## which collation was performed
                                            distThreshold = 0.75,
                                            ## the default value for
                                            ## - euclidean distance: 3
                                            ## - correlation distance: 0.75
                                            verbose = FALSE){
    ##
    if(minClusters < 2){
        warning("minClusters < 2. Setting it to 2")
        minClusters <- 2
    }
    .msg_pstr("minClusters is ", minClusters, flg=verbose)
    ##
    if(min(distMat) > distThreshold){
        clust_list <- lapply(hcObj$order, function(x) x)
        .msg_pstr("No element pairs close enough by given dist threshold: ",
                distThreshold, flg=verbose)
        .msg_pstr("#Clusters: ", length(clust_list), flg=verbose)
        return(clust_list)
    }else{
        ##
        cut_heights <- seq(min(hcObj$height), max(hcObj$height), by = hStep)
        clust_list <- .get_clusts_sil_or_ch(cut_heights, hcObj, distMat,
                        minClusters, use_sil = TRUE, verbose = verbose)
        ##
        clust_list <- fetch_cutree_by_hc_order(clust_list = clust_list,
                                                hcorder = hcObj$order)
        ##
        .msg_pstr("#Clusts using sil.vals:", length(clust_list), flg=verbose)
        .msg_pstr(.msg_print(clust_list), flg=verbose)
        clust_list_lengths <- unlist(lapply(clust_list, length))
        if(enableSwitchSilToCH && any(clust_list_lengths == 1)){
            .msg_pstr("But singleton(s) present: clusters ",
            paste(which(clust_list_lengths == 1), collapse= " "), flg=verbose)
            .msg_pstr("Using Calinski-Harabasz index instead", flg=verbose)
            ### Calinski-Harabasz Index
            clust_list <- .get_clusts_sil_or_ch(cut_heights, hcObj, distMat,
                minClusters, use_sil = FALSE, verbose = verbose)
        }
        ## number of clusters decided, get final clustering result
        .msg_pstr("Final #Clusters: ", length(clust_list), flg=verbose)
        .msg_pstr(.msg_print(clust_list), flg=verbose)
        ##
        if(!is.null(parentChunks) && keepSiblingsUncollated){
            clust_list <- .check_and_uncollate_siblings(clust_list,
                            parentChunks, verbose)
        }
        .msg_pstr("#Clusters: ", length(clust_list), flg=verbose)
        .msg_pstr(paste(clust_list), flg=verbose)
        return(clust_list)
    }
}
## =============================================================================

.detect_just_for_sake_clust <- function(cheight_idx, cl_list, vrbs=FALSE){
    upd_cl_list <- cl_list
    cl_lens <- lapply(cl_list, length)
    if(cheight_idx == 1 && length(which(cl_lens == 2)) == 1){
        .msg_pstr("just for sake clustering detected!", flg=vrbs)
        upd_cl_list <- .unfurl_nodeList(cl_list, vec_ver=TRUE, vrbs)
    }
    upd_cl_list
}
## =============================================================================

.check_and_uncollate_siblings <- function(clust_list, parentChunks, verbose){
    .msg_pstr("Uncollating siblings, in case...", flg=verbose)
    .msg_pstr("Checking parent chunks", flg=verbose)
    childrenPerParent <- lapply(unique(parentChunks), function(x){
        which(parentChunks == x)
    })
    updated_clust_list <- lapply(seq_along(clust_list), function(x){
        parents <- unique(parentChunks[clust_list[[x]]])
        thisX <- clust_list[[x]]
        if(length(parents) == 1 && length(thisX) > 1){
            ## if cluster elements come from same parent chunk of previous
            ## iteration + now check if all elements in this cluster
            ## are the only children of that parent chunk (in other words
            ## nothing got separated and combined with other factors)
            childrenThisParent <- childrenPerParent[[as.integer(parents)]]
            if(identical(as.numeric(childrenThisParent), as.numeric(thisX))){
                .msg_pstr("Parents are identical", flg=verbose)
                ## Update: 2020-12-13
                ## See update above.
                ## With keepSiblingsUncollated as FALSE, this check is not
                ## performed.
                ##
                ## Split the merge back into separate clusters
                return(lapply(childrenThisParent, function(y){y}))
            }else{
                .msg_pstr("returning thisX", flg=verbose)
                return(thisX)
            }
        }
        else{
            ## cluster elements come from different parent chunks
            ## (of previous iteration)
            return(thisX)
        }
    })
    clust_list <- .unfurl_nodeList(updated_clust_list, vec_ver=FALSE)
    clust_list
}
## =============================================================================

## This works best/is best referred with ward.D linkage
# @title Get clusters out of hierarchical clustering using either silhouette
# value or Calinski-Harabasz index.
# Used internally by seqArchR.
#
# @description Clusters from a hierarchical clustering object are obtained
# by using cutree at different heights of the tree. The optimum number of
# clusters are decided based on the average silhouette value of the clustering.
#
# @param hcObj The hierarchical clustering object as returned by
# \code{\link[stats]{hclust}}.
# @param distMat The distance matrix that was used by hclust, a
# \code{\link[stats]{dist}} object.
# @param hStep Numeric. The step size used to increment height values for
# cutree. Default value is 0.05.
# @param parentChunks List. Specify the factor numbers in the previous
# iteration of seqArchR that factors in the current iteration resulted from.
# Default value is NULL.
# @param keepSiblingsUncollated Logical. Specify TRUE if all clusters from a
# single parent chunk should not be collated, even when recommended so based
# on silhouette value computation. Default value FALSE. i.e. if they get
# collated, it is left as collated.
# @param enableSwitchSilToCH Logical. When and if the clustering result has
# singletons, using the average silhouette value may not be reliable to decide
# on the optimum number of clusters. In this case, if this argument is TRUE,
# the Calinski-Harabasz index is used instead. Setting it to FALSE, the
# switch is disabled, i.e., average silhouette value is used. Default is FALSE.
# @param minClusters Integer. Specify the minimum number of clusters to be
# fetched from the HAC solution.
# @param distThreshold Numeric. Specify a threshold of units of distance for
# considering any two elements as close enough to be merged in the same
# cluster. The default value is specified for Euclidean distance.
# @param verbose Logical. Specify TRUE for verbose output.
#
# @importFrom stats cutree
# @importFrom fpc calinhara
# @importFrom cluster silhouette
#
.get_clusts_sil_or_ch <- function(cut_heights, hcObj, distMat, minClusters,
                            use_sil = TRUE, verbose = FALSE){
    measure_cut_h <- unlist(lapply(cut_heights, function(x){
        foo_try <- stats::cutree(hcObj, h = x)
        names(foo_try) <- NULL
        if(length(unique(foo_try)) >= minClusters){
            if(use_sil){
                sils <- cluster::silhouette(foo_try, dist = distMat)
                retVal <- base::mean(sils[, "sil_width"])
            }else{
                retVal <- fpc::calinhara(distMat, foo_try)
            }
            retVal
        }else return(-100)
    }))

    ## multiple matches, first match index is returned with which.max
    cheight_idx <- which.max(measure_cut_h)
    score_str <- ifelse(use_sil, "sil score", "CH index")
    .msg_pstr("Max.", score_str, "at index =", cheight_idx, ", h:",
        cut_heights[cheight_idx], flg=verbose)

    ## number of clusters decided, get final clustering result
    cut_result <- stats::cutree(hcObj, h = cut_heights[cheight_idx])
    names(cut_result) <- NULL
    clust_list <- lapply(seq_along(unique(cut_result)),
        function(x) which(cut_result == x))
    clust_list
}
## =============================================================================


.regularizeMat <- function(basisMat, topN = 10){
    basisMat2 <- basisMat
    for(i in seq_len(ncol(basisMat))){
        asVec <- as.vector(basisMat[,i])
        threshold <- utils::tail(utils::head(
            sort(asVec, decreasing = TRUE), topN),1)
        basisMat2[(basisMat[,i] < threshold), i] <- 0.0
    }
    return(basisMat2)
}
## =============================================================================

#' @title Collate raw clusters at the chosen iteration of seqArchR result
#'
#' @description We use hierarchical clustering for reordering/collating raw
#' clusters from seqArchR's given iteration.
#'
#' @param result The seqArchR result object.
#'
#' @param iter Specify clusters at which iteration of seqArchR are to be
#' reordered/collated. Default is the last iteration of the seqArchR result
#' object.
#'
#' @param clust_method Specify 'hc' for hierarchical clustering. Currently, only
#' hierarchical clustering is supported.
#'
#' @param aggl_method One of linkage values as specified for hierarchical
#' clustering with \code{\link[stats]{hclust}}. Default is 'ward.D'.
#'
#' @param dist_method Distance measure to be used with hierarchical clustering.
#' Available options are "euclid" (default), "cor" for correlation, "cosangle"
#' for cosine angle, "modNW" for modified Needleman-Wunsch similarity (see
#' \code{\link[TFBSTools]{PFMSimilarity}}).
#'
#' @param regularize Logical. Specify TRUE if regularization is to be performed
#' before comparison. Default is FALSE. Also see argument 'topN'.
#'
#' @param topn Use only the top N dimensions of each basis vector for
#' comparing them. Note that since each basis vector has 4L or 16L (mono- or
#' dinucleotides) dimensions, each dimension is a combination of nucleotide and
#' its position in the sequence. This argument selects the top N dimensions of
#' the basis vector. This is ignored when argument 'regularize' is FALSE.
#'
#' @param collate Logical. Specify TRUE if collation using
#' hierarchical agglomerative clustering is to be performed, otherwise FALSE.
#'
#' @param return_order Logical. Use this argument when you want hierarchical
#'  clustering to be performed but not collation of clusters. Therefore,
#'  setting return_order to TRUE will return the hierarchical clustering
#'  object itself. This enables custom downstream processing/analysis.
#'
#' @param flags Pass the flags object similar to the flags in configuration
#' of the seqArchR result object.
#'
#' @param ... ignored
#'
#' @return When `collate` is TRUE, a list with the following elements is
#' returned:
#' \describe{
#' \item{basisVectorsCLust}{A list storing collation information of the
#' basis vectors, i.e, IDs of basis vectors that were collated into one.}
#' \item{clusters}{A list of sequences in each collated cluster.}
#' \item{seqClustLabels}{Cluster labels for all sequences according to the
#' collated clustering.}
#' }
#'
#' When 'collate' is FALSE,
#' it returns the already existing basis vectors, each as singleton clusters.
#' The sequence cluster labels and sequence clusters are also handled
#' accordingly. All are available as part of the same list as the earlier case.
#'
#' When 'return_order' is set to TRUE, the hierarchical clustering result is
#' returned instead.
#'
#' @examples
#' res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
#'          package = "seqArchR", mustWork = TRUE))
#'
#' # While the default settings for collation use Euclidean distance and
#' # ward.D agglomeration, one can choose to use different settings, say,
#' # correlation distance and complete linkage, and also regularizing to use
#' # only top 50 dimensions (nucleotide-positions combinations)
#' collated_res <- collate_seqArchR_result(result = res, iter = 2,
#'                         aggl_method = "complete", dist_method = "cor",
#'                         regularize = TRUE, topn = 50)
#'
#' names(collated_res)
#'
#' @importFrom stats hclust dist
#' @importFrom utils tail
#' @export
collate_seqArchR_result <- function(result,
                            iter = length(result$seqsClustLabels),
                            clust_method = "hc", aggl_method = "ward.D",
                            dist_method = "euclid", regularize = FALSE,
                            topn = 50, collate = TRUE, return_order = FALSE,
                            flags = list(debugFlag = FALSE, verboseFlag = TRUE),
                            ...) {
    ##
    dbg <- flags$debugFlag
    vrbs <- flags$verboseFlag

    stopifnot(iter > 0)
    basisMat <- get_clBasVec_m(result, iter=iter)
    clust_lab <- get_seqClLab(result, iter=iter)
    ##
    ## assert that topn value supplied is in the valid range
    if(!(topn > 0 && topn <= nrow(basisMat))){
        stop("Selected topn value ", topn," is outside expected range [",
                .msg_print(c(1, nrow(basisMat))), "]")
    }
    if(!is.logical(regularize) && !is.logical(collate)){
        stop("Arguments 'regularize' and 'collate' expect a logical ",
                " value (TRUE/FALSE), found otherwise")
    }
    ##
    if(regularize) basisMat <- .regularizeMat(basisMat, topn)
    ##
    clust_lab_pre <- NULL
    if(iter > 1) clust_lab_pre <- get_seqClLab(result, iter=iter-1)
    parentChunks <- .get_parent_chunks(clust_lab, clust_lab_pre, iter, dbg)
    ##
    if(collate){
        .msg_pstr("Collating clusters", flg=dbg)
        factorsClustering <-
            .handle_clustering_of_factors(basisMat, clustMethod = clust_method,
                    linkage = aggl_method, distMethod = dist_method,
                    flags = flags, parentChunks = parentChunks, ...)
        .msg_pstr("Factor clustering done, returning:", flg=dbg)
        .msg_pstr(paste("Returning: ", .msg_print(factorsClustering)), flg=dbg)
    }else{
        ## Do not collate, but prepare the return object
        nFactors <- get_clBasVec_k(result, iter)
        factorsClustering <- vector("list", nFactors)
        factorsClustering <- lapply(seq_len(nFactors), function(x){x})
        .msg_pstr("No factor clustering, returning:", flg=dbg)
        .msg_pstr(.msg_print(factorsClustering), flg=dbg)
    }
    if(return_order){
        .msg_pstr("Reordering clusters", flg=dbg)
        ## Do not collate just reorder, but prepare the return object
        nFactors <- get_clBasVec_k(result, iter)
        factorsClustering <- vector("list", nFactors)
        factorsClustering <- lapply(seq_len(nFactors), function(x){x})
        temp <- .handle_clustering_of_factors(basisMat,
                                        clustMethod = clust_method,
                                        linkage = aggl_method,
                                        distMethod = dist_method,
                                        returnOrder = return_order,
                                        flags = flags,
                                        parentChunks = parentChunks, ...)
        .msg_pstr("No collation of clusters, returning order:", flg=dbg)
        .msg_pstr(.msg_print(factorsClustering), flg=dbg)
        return(temp)
    }
    ##
    seqClusters <- get_seqs_clust_list(clust_lab)
    clusters <- collate_clusters(factorsClustering, seqClusters)
    clustLabels <- .update_cluster_labels(oldSeqsClustLabels = clust_lab,
                                collatedClustAssignments = clusters)

    cluster_sol <- list(basisVectorsClust = factorsClustering,
                        clusters = clusters, seqsClustLabels = clustLabels)
    return(cluster_sol)
}
## =============================================================================

.get_parent_chunks <- function(clust_lab, clust_lab_pre, iter, flg){
    ## Prepare info on parent chunks for each cluster/factor at given iteration
    ## We get this info from the seqsClustLabels
    parentChunks <- NULL
    if(iter > 1){
        clustsThisIter <- sort(unique(clust_lab))
        parentChunks <- unlist(lapply(clustsThisIter, function(x){
            relSeqsIds <- which(clust_lab == x)
            unique(clust_lab_pre[relSeqsIds])
        }))
    }
    parentChunks
}
## =============================================================================


#' @title Retrieve sequence clusters as a list from the sequence labels
#'
#' @description Given the sequence cluster labels from the seqArchR result
#' object, returns the clusters separated as a list.
#'
#' @param seqs_clust_lab Sequences with cluster labels as in the seqArchR
#' result object.
#'
#' @return A list holding sequence IDs belonging in each cluster.
#'
#' @examples
#'
#' clustLabels <- sample(seq_len(4), 50, replace = TRUE)
#' print(clustLabels)
#' get_seqs_clust_list(clustLabels)
#'
#' @export
get_seqs_clust_list <- function(seqs_clust_lab){
    ## check that labels are not empty/NULL
    .assert_seqArchR_seqsClustLabels(seqs_clust_lab)
    clusterLevels <- levels(as.factor(seqs_clust_lab))

    seqs_clusters_as_a_list <- lapply(seq_along(clusterLevels),
                                        function(x){
                                            which(seqs_clust_lab ==
                                                    clusterLevels[x])
                                    })

    return(seqs_clusters_as_a_list)
}
## =============================================================================
