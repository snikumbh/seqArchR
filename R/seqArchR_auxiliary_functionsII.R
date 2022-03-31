# @title Handle directory creation
#
# @description Given the output directory name with its complete path, this
# function checks if a directory of same name exists the given location. If
# yes, then it adds a suffix (a number) to the given directory name, and
# proceeds to create the directory.
#
# @param o_dir Specify the output directory name with its complete path.
#
# @param vrbs Set verbosity to TRUE or FALSE
#
# @return The (updated) dir name
#
#
handle_dir_creation <- function(o_dir, vrbs){
    ##
    cli::cli_alert_info(c("Output directory at path: ",
                    "{.emph {dirname(o_dir)}}"))
    if(dir.exists(o_dir)){
        # .msg_pstr("-- Directory exists: -- ", o_dir,
        #     "-- Changing name to: -- ", flg=vrbs)
        cli::cli_alert_warning("Directory exists: {.emph {basename(o_dir)}}")

        allExistingDirs <- list.dirs(path = dirname(o_dir),
                                        recursive = FALSE)
        dirsThatMatch <- grep(pattern = basename(o_dir), allExistingDirs,
                                value = TRUE)
        ## Suffix an integer, because directory with given name (oDir)
        ## exists for length(dirsThatMatch) times
        name_suffix <- length(dirsThatMatch)
        o_dir <- paste(o_dir, name_suffix, sep = "_")
        while(dir.exists(o_dir)){
            name_suffix <- name_suffix + 1
            o_dir <- paste(o_dir, name_suffix, sep = "_")
        }
    }
    retVal <- dir.create(o_dir, showWarnings = TRUE)
    stopifnot(retVal)
    cli::cli_alert_success("Writing to: {.emph {basename(o_dir)}}")
    o_dir
}
## =============================================================================


## Getter function to fetch the features matrix from NMF result object
## (from python)
##Dependency on python script perform_nmf.py
get_features_matrix <- function(nmfResultObj){
    returnVal <- .assert_seqArchR_list_properties(nmfResultObj)
    if (returnVal != "FOO") stop(returnVal)
    return(as.matrix(nmfResultObj[[1]]))
}
## =============================================================================

## Getter function to fetch the samples matrix from NMF result object
## (from python)
## Dependency on python script perform_nmf.py
get_samples_matrix <- function(nmfResultObj){
    returnVal <- .assert_seqArchR_list_properties(nmfResultObj)
    if (returnVal != "FOO") stop(returnVal)
    return(as.matrix(nmfResultObj[[2]]))
}
## =============================================================================

get_trimers_from_alphabet <- function(alph){
    if (is.null(alph)) stop("Expecting non-NULL alphabet")
    return(do.call(paste0, expand.grid(alph, alph, alph)))

}
## =============================================================================

get_dimers_from_alphabet <- function(alph){
    if (is.null(alph)) stop("Expecting non-NULL alphabet")
    return(do.call(paste0, expand.grid(alph, alph)))

}
## =============================================================================

manage_o_dir <- function(plt, o_dir){
    if(plt){
        if(is.null(o_dir)){
            stop("'plot' flag is TRUE but 'o_dir' is not provided. ",
                    "Did you forget to set 'o_dir'?")
        }
    }
}
## =============================================================================

# @importFrom Biostrings width
plot_all_seqs_logo <- function(seqs_raw, seqs_pos, dpath){
    if(is.null(seqs_raw)) stop("seqs_raw is NULL")
    if(is.null(dpath)) stop("directory path/name is NULL")
    if(is.null(seqs_pos)) seqs_pos <- seq(1, Biostrings::width(seqs_raw[1]))

    allSequencesLogo <- plot_ggseqlogo_of_seqs(
        seqs = seqs_raw,
        pos_lab = seqs_pos,
        title = paste("Sequence logo of all", length(seqs_raw),"sequences" ))
    ##

    ggsave(filename = file.path(dpath, "allSequencesLogo.pdf"),
    plot = allSequencesLogo,
    device = "pdf", width = 20, height = 2.5)

}
## =============================================================================


#' @title
#' Set seqArchR run configuration
#'
#' @description This function sets the configuration for `seqArchR`.
#'
#' @param chunk_size Numeric. Specify the size of the inner chunks of
#' sequences.
#' @param k_min Numeric. Specify the minimum of the range of values to be tested
#' for number of NMF basis vectors. Default is 1.
#' @param k_max Numeric. Specify the maximum of the range of values to be tested
#' for number of NMF basis vectors. Default is 50.
#' @param mod_sel_type Character. Specify the model selection strategy to
#' be used. Default is 'stability'. Another option is 'cv', short for
#' cross-validation. Warning: The cross-validation approach can be time
#' consuming and computationally expensive than the stability-based approach.
#' @param bound Numeric. Specify the lower bound value as criterion for choosing
#' the most appropriate number of NMF factors. Default is 1e-08.
#' @param cv_folds Numeric. Specify the number of cross-validation folds used
#' for model selection. Only used when mod_sel_type is set to 'cv'. Default
#' value is 5.
#' @param parallelize Logical. Specify whether to parallelize the procedure.
#' Note that running seqArchR serially can be time consuming, especially when
#' using cross-validation for model selection. See `n_cores`.
#' Consider parallelizing with at least 2 or 4 cores.
#' @param n_cores The number of cores to be used when `parallelize` is set
#' to TRUE. If `parallelize` is FALSE, nCores is ignored.
#' @param n_runs Numeric. Specify the number of bootstrapped runs
#' to be performed with NMF. Default value is 100. When using cross-validation
#' more than 100 iterations may be needed (upto 500).
#' @param alpha_base,alpha_pow Specify the base and the power for computing
#' 'alpha' in performing model selection for NMF. alpha = alpha_base^alpha_pow.
#' Alpha specifies the regularization for NMF. Default: 0 and 1 respectively.
#' _Warning_: Currently, not used (for future).
#' @param min_size Numeric. Specify the minimum number of sequences, such that
#' any cluster/chunk of size less than or equal to it will not be further
#' processed. Default is 25.
#' @param result_aggl Character. Specify the agglomeration method to be used
#' for final result collation with hierarchical clustering. Default is
#' 'complete' linkage. Possible values are those allowed with
#' \code{\link[stats]{hclust}}. Also see details below.
#' @param result_dist Character. Specify the distance method to be used for
#' final result collation with hierarchical clustering. Default is 'cor' for
#' correlation. Possible values are those allowed with
#' \code{\link[stats]{hclust}}. Also see details below.
#' @param checkpointing Logical. Specify whether to write intermediate
#' checkpoints to disk as RDS files. Checkpoints and the final result are
#' saved to disk provided the `o_dir` argument is set in \code{\link{seqArchR}}.
#' When `o_dir` argument is not provided or NULL, this is ignored.
#' Default is TRUE.
#' @param flags List with four logical elements as detailed.
#' \describe{
#'   \item{debug}{Whether debug information for the run is printed}
#'   \item{verbose}{Whether verbose information for the run is printed}
#'   \item{plot}{Whether verbose plotting is performed for the run}
#'   \item{time}{Whether timing information is printed for the run}
#' }
#'
#' @details Setting suitable values for the following parameters is dependent
#' on the data: 'inner_chunk_size', 'k_min', 'k_max', 'mod_sel_type',
#' 'min_size', 'result_aggl', 'result_dist'.
#'
#'
#' @return A list with all params for seqArchR set
#'
#' @examples
#' # Set seqArchR configuration
#' seqArchRconfig <- seqArchR::set_config(
#'     chunk_size = 100,
#'     parallelize = TRUE,
#'     n_cores = 2,
#'     n_runs = 100,
#'     k_min = 1,
#'     k_max = 20,
#'     mod_sel_type = "stability",
#'     bound = 10^-8,
#'     flags = list(debug = FALSE, time = TRUE, verbose = TRUE,
#'         plot = FALSE)
#' )
#'
#'
#' @export
set_config <- function(chunk_size = 500,
                            k_min = 1,
                            k_max = 50,
                            mod_sel_type = "stability",
                            bound = 10^-6,
                            cv_folds = 5,
                            parallelize = FALSE,
                            n_cores = NA,
                            n_runs = 100,
                            alpha_base = 0,
                            alpha_pow = 1,
                            min_size = 25,
                            result_aggl = "complete",
                            result_dist = "euclid",
                            checkpointing = TRUE,
                            flags = list(
                                debug = FALSE,
                                time = FALSE,
                                verbose = TRUE,
                                plot = FALSE)
                            ) {
    ## Configuration Params that can be set by user
    seqArchRconfig <- NULL
    ##
    if(is.null(flags)){
        useFlags <- list(
            debugFlag = FALSE,
            timeFlag = FALSE,
            verboseFlag = TRUE,
            plotVerboseFlag = FALSE)
    }else{
        useFlags <- list(
            debugFlag = flags$debug,
            timeFlag = flags$time,
            verboseFlag = flags$verbose,
            plotVerboseFlag = flags$plot)
    }
    ##
    seqArchRconfig <- list(modSelType = mod_sel_type,
                        # tol = tol,
                        bound = bound,
                        kFolds = cv_folds,
                        parallelize = parallelize,
                        nCoresUse = n_cores,
                        nRunsUse = n_runs,
                        paramRanges = list(
                            alphaBase = alpha_base,
                            alphaPow = alpha_pow,
                            k_vals = seq(k_min, k_max, by = 1)
                        ),
                        chunkSize = chunk_size,
                        result_aggl = result_aggl,
                        result_dist = result_dist,
                        checkpointing = checkpointing,
                        minSeqs = min_size,
                        flags = useFlags
                    )
    .assert_seqArchR_config(seqArchRconfig)
    return(seqArchRconfig)
}
## =============================================================================


## @title Decide processing of outer chunk based on its size
##
## @description Function to make the decision on whether the given (outer)
## chunk should be processed
##
## @param minThreshold Numeric. This is the minSeqs param from seqArchR config
## @param lengthOfOC Numeric. This is the length of the outer chunk
## @param kFoldsVal Numeric. This is the kFolds param in seqArchR config
##
## If the lengthOfOC == 0, STOP
## If the minThreshold < 4*kFolds, STOP
## If the lengthOfOC < minThreshold, DO NOT PROCESS/return TRUE
##
.decide_process_outer_chunk <- function(minThreshold, lengthOfOC, kFoldsVal) {
        stopifnot(lengthOfOC > 0)
        stopifnot(minThreshold > 0)
        ## Assert that minThreshold > 4*kFoldsVal
        nFoldsCondition <- 4 * kFoldsVal
        .assert_seqArchR_min_size_independent(minThreshold)
        if (minThreshold < nFoldsCondition) {
            stop("'min_size' should be at least 4 times 'kFolds'")
        }
        ##
        doNotProcess <- FALSE
        if (lengthOfOC <= minThreshold) {
            doNotProcess <- TRUE
            .msg_pstr("Sorry, will not process this small a chunk: ",
                    lengthOfOC, flg=TRUE)
        }
        ##
        return(doNotProcess)
    }
## =============================================================================

## @param factorsMat A matrix holding the factors along the columns
## @param distMethod character A string specifying the distance measure to
## be computed. Values are: 'modNW' for modified Needleman Wunsch, and any
## distance measure that is possible with HOPACH
##
## Default value 'modNW'
## Edit on 2021-01-02:
## Default value changed to euclid instead of modNW which is computed using
## a suggested package
##
## @return distance matrix from hopach (hdist object)
.compute_factor_distances <- function(factorsMat, distMethod = "euclid"){
    ## Assumption: Each column is a factor
    .assert_seqArchR_featuresMatrix(factorsMat)
    ## Since the default distMethod is euclid/euclidean, when the user wishes
    ## to use any other distance methods/metrics, we can check if hopach exists
    ## If not, we ask the user to install it.
    ## This lets move hopach to Suggests
    distMethods_hopach <- c("cosangle", "abscosangle",
                            "abseuclid", "cor", "abscor")
    distMethods_stats <- c("euclid")
    ##
    if(distMethod == "modNW"){
        distMat <- .get_modNW_dist(factorsMat)
        return(distMat)
    }
    if(any(distMethod == distMethods_hopach)){
        distMat <- .get_hopach_dist(factorsMat, distMethod)
        return(distMat)
    }
    if(distMethod == distMethods_stats){
        distMat <- .get_stats_dist(factorsMat)
        return(distMat)
    }
}
## =============================================================================

.get_hopach_dist <- function(factorsMat, distMethod){
    if(!requireNamespace("hopach", quietly = TRUE)){
        stop("Please install R package 'hopach' to use ", distMethod,
            " distance.")
    }else{
        ## - these distance metrics are available in hopach pkg
        ## - hopach::distancematrix func requires vectors along rows.
        ## - Distances are computed between row vectors
        if (nrow(factorsMat) > ncol(factorsMat)){
            factorsMat <- t(factorsMat)
        }
        hopachDist <- hopach::distancematrix(factorsMat, d = distMethod)
        ## hopachDistMat is a hopach hdist object
        stopifnot(length(hopachDist) == nrow(factorsMat))
        ## make as.matrix as done for dist object in the
        ## stats::dist case (see Else condition next)
        hopachDistMat <- hopach::as.matrix(hopachDist)

        return(hopachDistMat)
    }
}
## =============================================================================

.get_stats_dist <- function(factorsMat){
    ## dist method from stats // standard
    if (nrow(factorsMat) > ncol(factorsMat)) factorsMat <- t(factorsMat)
    as_dist <- stats::dist(factorsMat, method = "euclidean")
    distMat <- as.matrix(as_dist)
    return(distMat)
}

.get_modNW_dist <- function(factorsMat){
    if(!requireNamespace("TFBSTools", quietly = TRUE)){
        stop("Please install R package 'TFBSTools' for using modNW distance.")
    }else{
        ## Turn the factors which are vectors into a 2D matrix of
        ## dinucs x positions
        dim_names <- get_dimers_from_alphabet(c("A", "C", "G", "T"))
        nPositions <- nrow(factorsMat)/length(dim_names)
        ##
        # factorsMatList_as2D <- lapply(seq_len(ncol(factorsMat)),
        #     function(x){matrix(factorsMat[,x],
        #                     nrow = nrow(factorsMat)/nPositions,
        #                     byrow = TRUE,
        #                     dimnames = list(dim_names))
        #     })
        # ##
        # factorsMatList_asPFMs <- lapply(seq_len(length(factorsMatList_as2D)),
        #         function(x){
        #             sinucSparse <- collapse_into_sinuc_matrix(
        #                 given_feature_mat = as.matrix(factorsMat[,x]),
        #                 dinuc_mat = factorsMatList_as2D[[x]],
        #                 feature_names = dim_names)
        #             sinucSparseInt <- matrix(as.integer(round(sinucSparse)),
        #                 nrow = 4, byrow = FALSE,
        #                 dimnames = list(rownames(sinucSparse)))
        #         })
        ## After collapse_into_sinuc_matrix was updated, the above is changed
        ## to:
        factorsMatList_asPFMs <- lapply(seq_len(ncol(factorsMat)),
            function(x){
                ## A 16 x nPositions 2D matrix is obtained
                sinucSparse <- make_dinuc_PWMs(factorsMat[,x],
                                            add_pseudo_counts = FALSE,
                                            scale = FALSE)
                sinucSparseInt <- matrix(as.integer(round(sinucSparse)),
                                    nrow = 4, byrow = FALSE,
                                    dimnames = list(rownames(sinucSparse)))
            })

        ##
        lenPFMs <- length(factorsMatList_asPFMs)
        scoresMat <- matrix(rep(0, lenPFMs*lenPFMs), nrow = lenPFMs)
        rownames(scoresMat) <- seq(1,nrow(scoresMat),by=1)
        colnames(scoresMat) <- seq(1,ncol(scoresMat),by=1)

        for(i in seq_len(lenPFMs)){
            for(j in seq_len(lenPFMs)){
                temp <- TFBSTools::PFMSimilarity(
                    factorsMatList_asPFMs[[i]], factorsMatList_asPFMs[[j]])
                scoresMat[i,j] <- temp["score"]
                # relScoresMat[i,j] <- temp["relScore"]
            }
        }
        ## currently we use scoresMat, so we only return that
        distMat <- max(scoresMat) - scoresMat
    }
    return(distMat)
}
## =============================================================================


## For hierarchical clustering object, return the cluster medoids
## We currently use the first element of the cluster as its medoid
.get_factors_from_factor_clustering2 <- function(listObj, globFactorsMat){
    ##
    .assert_seqArchR_featuresMatrix(globFactorsMat)
    if (is.null(listObj)) {
        return(globFactorsMat)
    } else {
        medoids <- unlist(lapply(listObj, function(x){x[1]}))
        return(as.matrix(globFactorsMat[ , medoids]))
    }
}
## =============================================================================


## @title Process a chunk wih NMF
##
##
## On the given inner chunk,
## 1. perform model selection for #factors for NMF
## 2. Perform final NMF with chosen best_k (#Factors)
## 3. Store factors (globFactors)
## 4. Fetch clusters using k-means clustering
##      - assign clusters <--> factors
## 5. Store cluster assignments (globClustAssignments)
## 6. Return updated globFactors, globClustAssignments
##
## - innerChunkIdx is needed to appropriately index into the global variables:
##    globFactors and globClustAssignments
## - globClustAssignments variable is updated inside the function to
## additionally hold new ones
## - Similarly, globFactors variable is updated inside the function to
## additionally hold new ones
##
##
.handle_chunk_w_NMF2 <- function(innerChunkIdx,
                                    innerChunksColl,
                                    this_mat,
                                    cgfglinear = TRUE,
                                    coarse_step = 10,
                                    askParsimony = TRUE,
                                    config, test_itr, oChunkIdx,
                                    bpparam){
    .assert_seqArchR_flags(config$flags)
    dbg <- config$flags$debugFlag
    vrbs <- config$flags$verboseFlag
    tym <- config$flags$timeFlag
    ##
    if (is.null(this_mat) || !is.matrix(this_mat) &&
        !is(this_mat, "dgCMatrix")) {
        stop("Input matrix to model selection procedure is NULL or not a
            matrix")
    }
    #########################
    if(config$modSelType == "cv"){
        suffix_str <- ", without parsimony"
        if(askParsimony) suffix_str <- ", with parsimony"
        .msg_pstr("Performing cross validation-based model selection",
                suffix_str, flg=dbg)
        best_k <- .cv_model_select_pyNMF2(
                X = this_mat, param_ranges = config$paramRanges,
                kFolds = config$kFolds, nRuns = config$nRunsUse,
                # parallelDo = config$parallelize, nCores = config$nCoresUse,
                verboseFlag = config$flags$verboseFlag,
                debugFlag = config$flags$debugFlag,
                returnBestK = TRUE, cgfglinear = cgfglinear,
                coarse_step = coarse_step,
                askParsimony = askParsimony, bpparam = bpparam
            )
    }
    #########################
    if(config$modSelType == "stability"){
        .msg_pstr("Performing stability-based model selection", flg=dbg)
        best_k <- .stability_model_select_pyNMF2(
            X = this_mat, param_ranges = config$paramRanges,
            # parallelDo = config$parallelize, nCores = config$nCoresUse,
            nRuns = config$nRunsUse, bound = config$bound,
            flags = config$flags, returnBestK = TRUE, bootstrap = TRUE,
            bpparam = bpparam
        )
    }
    #########################
    if (best_k == max(config$paramRanges$k_vals)) {
        warning(c("Best K for this subset == 'k_max'. ",
                    "Consider selecting a larger 'k_max' value, or\n",
                    "smaller chunk_size, or\n",
                    "perhaps, further increasing 'n_runs'\n"),
                immediate. = TRUE)
    }
    cli::cli_alert_info("Best K for this chunk: {best_k}")

    ##
    if (best_k >= 1) {
        ## For fetching sequence clusters from samplesMat
        ## Cluster sequences
        ## New strategy, perform nRuns for bestK and use only the best one
        nRuns <- config$nRunsUse
        # .msg_pstr("Fetching ", best_k, " clusters", flg=(vrbs || dbg))

        featuresMatrixList <- vector("list", nRuns)
        samplesMatrixList <- vector("list", nRuns)
        new_ord <- vector("list", nRuns)
        nmf_nRuns_list <- .perform_multiple_NMF_runs(X = this_mat,
                                        kVal = best_k,
                                        alphaVal = 0,
                                        # parallelDo = config$parallelize,
                                        # nCores = config$nCoresUse,
                                        nRuns = nRuns,
                                        bootstrap = TRUE, bpparam = bpparam)
        featuresMatrixList <- lapply(nmf_nRuns_list$nmf_result_list,
                                        get_features_matrix)
        samplesMatrixList <- lapply(nmf_nRuns_list$nmf_result_list,
                                        get_samples_matrix)
        ##
        new_ord <- nmf_nRuns_list$new_ord
        ## Get reconstruction accuracies for them
        bestQ2 <- -1
        for (nR in seq_len(nRuns)){
            ##A <- this_mat[, new_ord[[nR]]]
            A <- this_mat
            recA <- as.matrix(featuresMatrixList[[nR]]) %*%
                as.matrix(samplesMatrixList[[nR]])
            this_q2 <- .compute_q2(as.matrix(A), recA)
            if(this_q2 > bestQ2){
                bestQ2 <- this_q2
                bestFeatMat <- featuresMatrixList[[nR]]
                bestSampMat <- samplesMatrixList[[nR]]
                bestOrd <- new_ord[[nR]]
            }
        }
        # .msg_pstr("Best Q2 giving run found: ", bestQ2, flg=dbg)
        # cli::cli_alert_info("Fetched {best_k} clusters")
        ##
        featuresMatrix <- bestFeatMat
        samplesMatrix <- bestSampMat
        ## When order was changing:
        ## put samples matrix back in order it should be
        tempM <- bestSampMat
        samplesMatrix <- matrix(rep(NA, length(tempM)), nrow = nrow(tempM))
        samplesMatrix[ ,bestOrd] <- tempM
        #####
        # .msg_pstr("Fetching ", best_k," cluster(s)", flg=dbg)
        clusterMembershipsForSamples <-
            .get_cluster_memberships_per_run(samplesMatrix = samplesMatrix,
                iChunksColl = innerChunksColl, iChunkIdx = innerChunkIdx,
                test_itr, oChunkIdx)
        ##
        ## Could handle overfitting here
        if(best_k > 1){
            has_overfit <- .detect_overfitting(samplesMatrix,
                                                clusterMembershipsForSamples,
                                                minSeqs = 50)
            ##
            ## -- Note which clusters are overfit
            ## -- Remove those columns from featuresMat
            ## -- Adjust clustMemberships
            if(length(has_overfit) > 0){
                clusterMembershipsForSamples <-
                    .adjustSampleMemberships(clusterMembershipsForSamples,
                                        samplesMatrix, has_overfit)
                featuresMatrix <- as.matrix(featuresMatrix[, -c(has_overfit)])
                best_k <- best_k - length(has_overfit)
                cli::cli_alert_info(c("Adjusting for overfitting, ",
                                    "fetched {best_k} cluster{?s}"))
                stopifnot(best_k == ncol(featuresMatrix))
                stopifnot(best_k ==
                            length(unique(clusterMembershipsForSamples)))
                if(best_k == 1){
                    stopifnot(unique(clusterMembershipsForSamples) == 1)
                }
            }

        }
        ##
        forGlobClustAssignments <- .assign_samples_to_clusters(
            clusterMembershipsVec = clusterMembershipsForSamples,
            nClusters = best_k, iChunkIdx = innerChunkIdx,
            iChunksColl = innerChunksColl)
        ##
    } else if (best_k < 1) {
        stop("Chosen number of factors: ", best_k)
    }
    .assert_seqArchR_featuresMatrix(featuresMatrix)
    .assert_seqArchR_globClustAssignments(forGlobClustAssignments)
    innerChunkNMFResult <- list(forGlobFactors = featuresMatrix,
                                forGlobClustAssignments =
                                    forGlobClustAssignments)
    ##
    return(innerChunkNMFResult)
}
## =============================================================================



.getMeanOfListOfMatrices <- function(listOfMats) {
    # Currently, assume all matrices have same dimensions
    # if discrepancy in dimensions, throws error "non-conformable arrays"
    meanMat <- Reduce("+", listOfMats)/length(listOfMats)
    return(meanMat)
}


## @title Setup the clustFactors list element for seqArchR result object
##
## @description Function to set up the clustFactors variable for seqArchR result
## object. Having a separate dedicated function enables seamless future changes.
##
## @param intClustFactors A matrix holding all factors from the just concluded
## iteration along the columns.
##
## @return A list with 2 elements having fixed element names. They are
## 'nBasisVectors':  This is the number of basis vectors (for easy info access)
##  and
## 'basisVectors': This is the matrix of basis vectors (intClustFactors itself)
.setup_clustFactors_for_seqArchR_result <- function(intClustFactors) {
    returnList <- list(nBasisVectors = ncol(intClustFactors),
                        basisVectors = intClustFactors)
    return(returnList)
}
## =============================================================================


intermediateResultsPlot <- function(seq_lab, seqs_raw = NULL,
                                pos_lab = NULL, iter = 0, fname = NULL,
                                name_suffix = NULL,
                                vrbs = TRUE){
    ## This function plots and prints resulting clusters -- the sequence image
    ## matrix (PNG file) and the sequence logos (PDF file).

    if(is.null(name_suffix)){
        if(is.numeric(iter)){
            name_suffix <- paste0("Iteration", iter)
            cli::cli_h2("Intermediate result")
        }else{
            name_suffix <- paste0("Final")
            cli::cli_rule(left="Final result")
        }
    }else{
        .msg_pstr("=== On-demand Result ===", flg=vrbs)
    }
    ##
    cli::cli_alert_info("Output directory: {.emph {fname}}")

    seqs_clust_list_ord <- get_seqs_clust_list(seq_lab)
    seqs_clust_vec_ord <- unlist(seqs_clust_list_ord)
    image_fname <- file.path(fname,
        paste0("ClusteringImage_", name_suffix, ".png"))
    cli::cli_alert_info(c("Sequence clustering image written to: ",
                            "{.emph {basename(image_fname)}}"))
    viz_seqs_acgt_mat(
        seqs =  as.character(seqs_raw[seqs_clust_vec_ord]),
                        pos_lab = pos_lab,
                        save_fname = image_fname,
                        f_width = 450,
                        f_height = 900,
                        xt_freq = 5,
                        yt_freq = 100)


    arch_fn <- file.path(fname,
        paste0("Architecture_SequenceLogos_", name_suffix, ".pdf"))
    cli::cli_alert_info(c("Architectures written to: ",
                            "{.emph {basename(arch_fn)}}"))

    plot_arch_for_clusters(seqs = seqs_raw,
                            clust_list = seqs_clust_list_ord,
                            pos_lab = pos_lab,
                            pdf_name = arch_fn)
}
## =============================================================================

save_final_result <- function(o_dir, temp_seqArchRresult){
    if(!is.null(o_dir)){
        rdsFilename <- file.path(o_dir, "seqArchRresult.rds")
        saveRDS(temp_seqArchRresult, file=rdsFilename)
        cli::cli_alert_info("Result saved to: {.emph {basename(rdsFilename)}}")
    }
}
## =============================================================================

save_checkpoint <- function(o_dir, test_itr, threshold_itr,
                            seqsClustLabelsList, clustFactors,
                            seqs_raw, config, call = NULL){
    ##
    if(!is.null(o_dir) && test_itr != threshold_itr){
        ##
        curr_seqArchRresult <- list(seqsClustLabels = seqsClustLabelsList,
                                clustBasisVectors = clustFactors,
                                rawSeqs = seqs_raw,
                                config = config,
                                call = call)
        rdsFilename <- file.path(o_dir, paste0("seqArchRresult_checkpoint",
                                test_itr, ".rds"))
        saveRDS(curr_seqArchRresult, file=rdsFilename)
        cli::cli_alert_info(c("Checkpointing at iteration {test_itr}: ",
            "{basename(rdsFilename)}"))
    }

}
## =============================================================================


show_ellapsed_time <- function(use_str = "Time ellapsed since start: ",
                                use_time){
    complTime1 <- as.difftime(Sys.time() - use_time)
    cli::cli_alert(
        c(use_str, "{prettyunits::pretty_dt(complTime1)}"))
    as.double.difftime(complTime1)
}
## =============================================================================

decisionToCollate <- function(clustFactors, dbg){
    decisionToCollate <- TRUE
    iterations <- length(clustFactors)
    if(iterations > 1){
        prevIterNFactors <- clustFactors[[iterations-1]]$nBasisVectors
        currIterNFactors <- clustFactors[[iterations]]$nBasisVectors
        if(currIterNFactors == prevIterNFactors){
            decisionToCollate <- FALSE
            .msg_pstr("Reordering decision: FALSE", flg=dbg)
        }
    }
    decisionToCollate
}
## =============================================================================

keepMinClusters <- function(set_ocollation, temp_res = NULL,
                            totOuterChunksColl = NULL, test_itr = 1,
                            nClustEachOC = NULL, nClustEachIC = NULL,
                            dbg, clustFactors = NULL, stage = "Final"){

    if(!is.null(stage) && stage == "Final" && test_itr > 1){
        setMinClustersFinal <- 2 ## the default value
        ## For setting minClusters, note last iteration collated
        if(any(set_ocollation)){
            lastItrC <- tail(which(set_ocollation), 1)
            setMinClustersFinal <- get_clBasVec_k(temp_res, lastItrC)
                # temp_res$clustBasisVectors[[lastItrC]]$nBasisVectors
        }
        return(setMinClustersFinal)
    }

    if(totOuterChunksColl > 1){
        .msg_pstr("meanClustersOC: ", ceiling(mean(nClustEachOC)), flg=dbg)
        ## For setting minClusters, note last iteration collated
        chkIdx <- seq_len(test_itr-1)
        if(any(set_ocollation[chkIdx])){
            lastItrC <- tail(which(set_ocollation[chkIdx]), 1)
            setMinClusters <- clustFactors[[lastItrC]]$nBasisVectors
        }else{
            ## average clusters identified in each chunk of the 1st iter
            setMinClusters <- max(ceiling(mean(nClustEachOC[1])), 2)
        }
    }else{
        ## When totOuterChunks is == 1, this is the first iteration
        ## Use the mean of nClustEachIC
        .msg_pstr("meanClustersIC: ", ceiling(mean(nClustEachIC)), flg=dbg)
        setMinClusters <- max(ceiling(mean(nClustEachIC)), 2)
    }
    return(setMinClusters)
}
## =============================================================================


perform_setup <- function(config, total_itr, o_dir, fresh,
                        seqs_pos, seqs_raw, seqs_ohe_mat,
                        set_ocollation){
    dbg <- config$flags$debugFlag
    vrbs <- config$flags$verboseFlag
    plt <- config$flags$plotVerboseFlag
    crs <- config$nCoresUse
    parallelize <- config$parallelize
    modSelType <- config$modSelType
    bound <- config$bound

    ## assert total_itr is a positive integer
    .assert_seqArchR_thresholdIteration(total_itr)
    ## TODO Provide a summary function

    if(!is.null(o_dir)){
        if(fresh){
            o_dir <- handle_dir_creation(o_dir, vrbs||dbg)
        }else{
            if(!dir.exists(o_dir)){
                stop(o_dir, " not found")
            }
        }
    }
    ##
    if(is.null(seqs_pos)){
        seqs_pos <- seq_len(Biostrings::width(seqs_raw[1]))
    }
    ##
    if(plt){
        manage_o_dir(plt, o_dir) # this will stop if o_dir is NULL
        plot_all_seqs_logo(seqs_raw, seqs_pos, dpath=o_dir)
    }
    ## Make checks for params in configuration
    .assert_seqArchR_config(config, seqs_size = ncol(seqs_ohe_mat))
    .assert_seqArchR_thresholdIteration(total_itr)
    .assert_seqArchR_collation(set_ocollation, total_itr)

    if(length(set_ocollation) > total_itr){
        set_ocollation <- set_ocollation[seq(1,total_itr)]
        cli::cli_alert_info(c("Changing length of 'set_ocollation' to ",
                                "same as that of total_itr"))
    }

    #### Start cluster only once -- using BiocParallel
    if(parallelize){
        if(.Platform$OS.type == "windows"){
            if (is.null(crs)) crs <- BiocParallel::multicoreWorkers()
            cl <- BiocParallel::SnowParam(workers = crs, tasks = crs,
                                            exportglobals = TRUE)
            cli::cli_alert_info("Parallelization: Yes")
        }else{
            if (is.null(crs)) crs <- BiocParallel::multicoreWorkers()
            cl <- BiocParallel::MulticoreParam(workers = crs, tasks = crs)
            cli::cli_alert_info("Parallelization: Yes")
        }
    }else{
        cl <- BiocParallel::SerialParam()
        cli::cli_alert_info("Parallelization: No")
    }

    # #### Start cluster only once --using parallel
    # if(parallelize){
    #     cl <- parallel::makeCluster(crs, type = "FORK")
    #     parallel::setDefaultCluster(cl)
    #     cli::cli_alert_info("Parallelization: {crs} cores")
    # }else{
    #     cl <- NA
    #     cli::cli_alert_info("Parallelization: No")
    # }
    ##
    if(modSelType != "cv" && modSelType != "stability")
        cli::cli_alert_warning(c("Mis-specified model selection strategy",
            "Using factor stability for model selection"))
    msg_suffix <- ifelse(modSelType == "cv", "cross-validation",
        "factor stability")
    cli::cli_alert_info("Model selection by {msg_suffix}")
    if(modSelType == "stability") cli::cli_alert_info("Bound: {bound}")
    set_parsimony <- NULL
    if(modSelType == "cv"){
        ## Specify set_parsimony when cv
        set_parsimony <- rep(TRUE, total_itr)
    }

    return(list(cl = cl, o_dir = o_dir, seqs_pos = seqs_pos,
                set_parsimony = set_parsimony))
}
## =============================================================================


process_innerChunk <- function(test_itr, innerChunksColl, config, lenOC,
                            seqs_ohe_mat, set_parsimony, outerChunkIdx,
                            bpparam){
    # globFactors <- vector("list", length(innerChunksColl))
    # globClustAssignments <- vector("list", length(innerChunksColl))
    # nClustEachIC <- rep(0, length(innerChunksColl))

    nmfResultEachIC <- lapply(seq_along(innerChunksColl), function(x){
        innerChunkIdx <- x
        cli::cli_h3(c("Inner chunk {innerChunkIdx} of ",
            "{length(innerChunksColl)} ",
            "[Size: {length(innerChunksColl[[innerChunkIdx]])}]"))
        ##
        ## Setting up sequences for the current chunk
        this_seqsMat <-
            seqs_ohe_mat[, innerChunksColl[[innerChunkIdx]]]
        ##
        chnksz <- config$chunkSize
        cvStep <- ifelse(test_itr == 1 || lenOC > 0.9*chnksz, 10, 5)
        thisNMFResult <- .handle_chunk_w_NMF2(innerChunkIdx,
            innerChunksColl, this_seqsMat,
            cgfglinear = TRUE, coarse_step = cvStep,
            askParsimony = set_parsimony[test_itr],
            config = config, test_itr = test_itr,
            oChunkIdx = outerChunkIdx, bpparam = bpparam)
        thisNMFResult
    })

    assertions <- lapply(nmfResultEachIC, .assert_seqArchR_NMFresult)

    globFactors <- lapply(nmfResultEachIC, function(x){x$forGlobFactors})

    globClustAssignments <- lapply(nmfResultEachIC, function(x){
                                    x$forGlobClustAssignments})

    nClustEachIC <- unlist(lapply(globClustAssignments, length))

    return(list(globFactors = globFactors,
                globClustAssignments = globClustAssignments,
                nClustEachIC = nClustEachIC))
}
## =============================================================================
