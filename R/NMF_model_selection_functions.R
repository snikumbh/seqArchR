# @title
# Get the reconstruction accuracy, \eqn{Q^{2}}, for a given choice of parameter
# values
#
# @param x The parameter values. Expected columns "k_vals", "alpha" and "fold".
# @param seed_val The seed to be used.
# @param verbose Default to \code{0} which will not print any messages, or can
# be set to \code{1} which will print messages.
#
# @return \eqn{Q^{2}}, a real value
# @importFrom MASS ginv
.get_q2_using_py <- function(x) {
    ##
    this_k <- as.numeric(x["k_vals"])
    this_alpha <- as.numeric(x["alpha"])
    this_seed <- as.numeric(x["seed_val"])
    #
    test_fold <- as.numeric(x["fold"])
    train_rows <-
        cvfolds$cvf_rows$subsets[cvfolds$cvf_rows$which != test_fold]
    train_cols <-
        cvfolds$cvf_cols$subsets[cvfolds$cvf_cols$which != test_fold]
    test_rows <-
        cvfolds$cvf_rows$subsets[cvfolds$cvf_rows$which == test_fold]
    test_cols <-
        cvfolds$cvf_cols$subsets[cvfolds$cvf_cols$which == test_fold]
    ##
    ## Split data matrix X -- separate the training and test parts
    ##        _      _
    ##   X = |  A  B  |
    ##       |_ C  D _|
    ##
    ## Reconstruct A, by performing NMF on D
    ## More details in Owen and Perry, Annals of Statistics, 2009
    ##
    ##
    submatrixD <- X[train_rows, train_cols]
    submatrixA <- X[test_rows, test_cols]
    submatrixB <- X[test_rows, train_cols]
    submatrixC <- X[train_rows, test_cols]

    ## NMF on submatrixD
    ## 1. Setup params
    ## 2. NMF call, using python/scikit-learn NMF
    nmf_submatrixD_result <- .perform_single_NMF_run(X = submatrixD,
                                kVal = as.integer(this_k),
                                alphaVal = this_alpha, seedVal = this_seed)
    D_W <- nmf_submatrixD_result$featuresMatrix
    D_H <- nmf_submatrixD_result$samplesMatrix
    ##
    reconstructed_submatrixA <- as.matrix(submatrixB) %*%
        MASS::ginv(D_H) %*% MASS::ginv(D_W) %*% as.matrix(submatrixC)
    ##
    q2 <- .compute_q2(as.matrix(submatrixA), reconstructed_submatrixA)
    return(q2)
}
## =============================================================================

.get_q2_using_py_serial <- function(x, X, cvfolds){

    ##
    this_k <- as.numeric(x["k_vals"])
    this_alpha <- as.numeric(x["alpha"])
    this_seed <- as.numeric(x["seed_val"])
    #
    test_fold <- as.numeric(x["fold"])
    train_rows <-
        cvfolds$cvf_rows$subsets[cvfolds$cvf_rows$which != test_fold]
    train_cols <-
        cvfolds$cvf_cols$subsets[cvfolds$cvf_cols$which != test_fold]
    test_rows <-
        cvfolds$cvf_rows$subsets[cvfolds$cvf_rows$which == test_fold]
    test_cols <-
        cvfolds$cvf_cols$subsets[cvfolds$cvf_cols$which == test_fold]
    ##
    ## Split data matrix X -- separate the training and test parts
    ##        _      _
    ##   X = |  A  B  |
    ##       |_ C  D _|
    ##
    ## Reconstruct A, by performing NMF on D
    ## More details in Owen and Perry, Annals of Statistics, 2009
    ##
    ##
    submatrixD <- X[train_rows, train_cols]
    submatrixA <- X[test_rows, test_cols]
    submatrixB <- X[test_rows, train_cols]
    submatrixC <- X[train_rows, test_cols]

    ## NMF on submatrixD
    ## 1. Setup params
    ## 2. NMF call, using python/scikit-learn NMF
    nmf_submatrixD_result <- .perform_single_NMF_run(X = submatrixD,
                                                    kVal = as.integer(this_k),
                                                    alphaVal = this_alpha,
                                                    seedVal = this_seed)
    D_W <- nmf_submatrixD_result$featuresMatrix
    D_H <- nmf_submatrixD_result$samplesMatrix
    ##
    reconstructed_submatrixA <- as.matrix(submatrixB) %*%
        MASS::ginv(D_H) %*% MASS::ginv(D_W) %*% as.matrix(submatrixC)
    ##
    q2 <- .compute_q2(as.matrix(submatrixA), reconstructed_submatrixA)
    return(q2)

}


.perform_single_NMF_run <- function(X, kVal, alphaVal, seedVal) {

    reticulate::source_python(system.file("python", "perform_nmf.py",
                                        package = "seqArchR",
                                        mustWork = TRUE)
    )
    ##
    nmf_result <- perform_nmf_func(X,
                            nPatterns = as.integer(kVal),
                            nIter = as.integer(2000),
                            givenAlpha = alphaVal,
                            givenL1_ratio = 1,
                            seed_val = as.integer(seedVal)
                            )
    D_W <- as.matrix(get_features_matrix(nmf_result))
    D_H <- as.matrix(get_samples_matrix(nmf_result))
    return(list(featuresMatrix = D_W, samplesMatrix = D_H))
}
## =============================================================================


# @title Compute \eqn{Q^2} Value
#
# @description Computes the reconstruction accuracy, \eqn{Q^2}, given the
# original matrix and the recontructed matrix to check against.
#
# @param A The original matrix
# @param recA The reconstructed matrix
#
# @return Q^2, the reconstruction accuracy
.compute_q2 <- function(A, recA) {
    ## input: original matrix, reconstructed portion
    ## return: computed q2 val
    if (!is.matrix(A)) {
        stop("Original matrix is not of type matrix")
    }
    if (!is.matrix(recA)) {
        stop("Reconstructed matrix is not of type matrix")
    }
    if (sum(dim(A)) == 2 && is.na(A)) {
        stop("Empty: Original matrix")
    }
    if (sum(dim(recA)) == 2 && is.na(recA)) {
        stop("Empty: Reconstructed matrix")
    }
    second_term_num <- sum((A - recA) ^ 2)
    second_term_den <- sum(A ^ 2) # TO-DO: make sure this is not zero
    if (second_term_den != 0.0) {
        q2 <- 1 - (second_term_num / second_term_den)
    } else {
        q2 <- 0.0
    }
    return(q2)
}
## =============================================================================

.grid_srch_par_df <- function(this_K, aBase, aPow, kFolds, nIter){
    dummy_seed <- 1234
    ## need this only as place filler for the seed_val column in the
    ## grid_search_params df
    grid_search_params <- expand.grid(list(
        k_vals = this_K,
        alpha = aBase ^ aPow,
        fold = seq_len(kFolds), iteration = seq_len(nIter),
        seed_val = dummy_seed
    ))
    ##
    seed_val_list <- get_n_seeds(n = kFolds*nIter)
    ##
    grid_search_params[, "seed_val"] <- seed_val_list
    return(grid_search_params)
}


get_n_seeds <- function(n){
    return(sample.int(.Machine$integer.max, size = n, replace = FALSE))
}

# should this be part of main help?
#  For 'result_aggl' and 'result_dist', it has been observed that
# agglomeration by 'ward.D' clustering and Euclidean distance work well
# together while complete linkage works well with 'correlation' as distance.
# The former works well for architectures that resemble those observed in
# Drosophila, and the later works well for architectures resembling those in
# mice.
#
performSearchForK <- function(X, cvfolds, startVal, endVal, step = 1,
                            prev_best_K = -1, best_K = 0, prev_df = NULL,
                            param_ranges, kFolds, nRuns,
                            parallelDo = FALSE, set_verbose = 1){
    vrbs <- ifelse(set_verbose == 1, TRUE, FALSE)
    dbg <- ifelse(set_verbose == 2, TRUE, FALSE)
    kValues <- seq(startVal, endVal, by = step)
    .msg_pstr("Checking K = [", paste(kValues, collapse = ","), "]", flg=vrbs)
    for (this_K in kValues){
        if(prev_best_K != best_K && !is.na(this_K) && this_K != 0){
            ##
            grid_search_params <- .grid_srch_par_df(this_K,
                #aBase=param_ranges$alphaBase, aPow=param_ranges$alphaPow,
                aBase = 0, aPow = 1,
                kFolds = kFolds, nIter = nRuns)
            ##
            if(parallelDo){
                q2_vals <- unlist(parallel::clusterApplyLB(cl = NULL,
                        seq_len(nrow(grid_search_params)), function(i) {
                            .get_q2_using_py( grid_search_params[i,] )
                        }))
            }else{
                q2_vals <- unlist(
                    lapply(seq_len(nrow(grid_search_params)),
                        function(i) {
                            .get_q2_using_py_serial( grid_search_params[i,],
                                                X = X, cvfolds = cvfolds)
                        }))
            }
            ##
            grid_search_results <- as.data.frame(grid_search_params[,
                    c("k_vals", "alpha", "fold", "iteration")],
                    col.names = c("k_vals", "alpha", "fold", "iteration"))
            grid_search_results$q2_vals <- q2_vals
            ##
            if (is.null(prev_df)) {
                best_K <- .get_best_K(grid_search_results)
                prev_df <- grid_search_results
            } else {
                new_df <- rbind(grid_search_results, prev_df)
                prev_best_K <- best_K
                best_K <- .get_best_K(new_df)
                prev_df <- new_df
            }
            .msg_pstr("Prev best K:", prev_best_K,
                        "Best K:", best_K,
                        "This K:", this_K, flg=dbg)
            .msg_pstr("Curr nrows:", nrow(grid_search_results),
                        "Total nrows:", nrow(prev_df), flg=dbg)
        }
    }
    returnObject <- list(best_K = best_K, prev_best_K = prev_best_K,
                        return_df = prev_df)
    returnObject
}


# @title Perform Model Selection via Cross-Validation
#
# @description The function performs model selection via cross-validation for
# choosing the optimal values of parameters for NMF. These are K, the number of
#  factors in NMF, and \eqn{\alpha}, the sparsity coefficient.
#
# @param X The given data matrix
# @param param_ranges An object holding the range of values for parameters
# \code{k}, \code{alphaBase}, and \code{alphaPow}.
# @param kFolds Numeric The number of cross-validation folds.
# @param parallelDo Set to \code{1} if you want to parallelize, else \code{0}.
# @param nRuns Number of runs for NMF.
# @param nCores If \code{parallel} is set to \code{1}
# @param seed_val The seed to be set.
# @param set_verbose Default to \code{0} which will not print any messages, or
# can be set to \code{1} which will print messages. Value passed to the
# \code{get_q2_val} function
#
# @return A data.frame/tibble of grid_search_results
#
# @importFrom methods is
# @importFrom parallel makeCluster stopCluster detectCores clusterEvalQ
# @importFrom parallel clusterExport getDefaultCluster
.cv_model_select_pyNMF2 <- function(X,
                                    param_ranges,
                                    kFolds = 5,
                                    parallelDo = FALSE,
                                    nCores = NA,
                                    nRuns = 20,
                                    returnBestK = TRUE,
                                    cgfglinear = TRUE,
                                    coarse_step = 10,
                                    askParsimony = FALSE,
                                    # monolinear = FALSE,
                                    debugFlag = FALSE,
                                    verboseFlag = TRUE) {
    dbg <- debugFlag
    vrbs <- verboseFlag
    if (!is.matrix(X) && !is(X, "dgCMatrix")) {
        stop("X not of type matrix/dgCMatrix")
    }
    ##
    if (kFolds < 3) {
        if (kFolds < 0) {
            stop("Number of cross-validation folds cannot be negative")
        } else {
            stop("Set at least 3 cross-validation folds")
        }
    } else if (kFolds > ncol(X)) {
        stop("CV folds should be less than or equal to #sequences.
                Standard values: 3, 5, 10.")
    }
    ## Check names in param_ranges list, the function relies on it below
    if (length(setdiff(names(param_ranges), c("alphaPow", "alphaBase",
                                    "k_vals"))) > 0) {
        stop("Expected elements in param ranges: alphaBase, alphaPow, k_vals")
    }
    ## Get cross-validation folds
    cvfolds <- .generate_folds(dim(X), kFolds)
    ##
    # if (parallelDo) {
        #######################
        #### New strategy
        if(returnBestK) {
            if(parallelDo){
                cl <- .setup_par_cluster(vlist=
                    c(".get_q2_using_py", ".compute_q2", "X", "cvfolds"))
            }
            ##
            set_verbose <- ifelse(debugFlag, 2, ifelse(verboseFlag, 1, 0))
            if(cgfglinear){
                # .msg_pstr("Coarse-fine grained binary search", flg=vrbs)
                prev_df <- NULL
                #go_fine <- FALSE
                ## when either lo or hi values are best, so we need to perform
                ## a fine-grained search
                #eureka <- FALSE
                ## to note when the mid value happens to be
                ## the best
                #coarse_step <- 10
                mi <- seq(coarse_step, max(param_ranges$k_vals), by=coarse_step)
                lo <- mi-1
                hi <- mi+1
                stopifnot(length(lo) == length(mi) && length(lo) == length(hi))
                ## For loop for coarse-grained search
                prev_best_K <- -1
                best_K <- 0
                for (kCGIdx in seq_along(lo)){
                    searchReturnCoarse <- performSearchForK(X = X,
                        cvfolds = cvfolds,
                        startVal = lo[kCGIdx],
                        endVal = hi[kCGIdx],
                        step = 1,
                        prev_best_K = prev_best_K,
                        best_K = best_K,
                        prev_df = prev_df,
                        param_ranges,
                        parallelDo = parallelDo,
                        kFolds = kFolds, nRuns = nRuns,
                        set_verbose = set_verbose)
                    best_K <- searchReturnCoarse$best_K
                    prev_best_K <- searchReturnCoarse$prev_best_K
                    prev_df <- searchReturnCoarse$return_df

                    if(best_K == mi[kCGIdx]){
                        ## Work done, leave loop?
                        .msg_pstr("mi value is best", flg=dbg)
                        if(askParsimony){
                            ## When mi value is chosen as best, in order for
                            ## the 1-SE rule, it would be a good idea compute
                            ## q2 values for at least 5 consecutive values of
                            ## K less than mi.
                            ## This means that, in thise case, we would set
                            ## fgIL to (mi-5)
                            ## Then, the 1-SE rule could lead to choosing a
                            ## smaller value.
                            ##
                            #go_fine <- TRUE
                            ## Scenario when the mi value is best, first check
                            ## the 1-SE rule,
                            ## What is the value chosen by this rule?
                            ## Does it already include values from the previous
                            ## triplet, if any?
                            ## If yes, we may need computations for a few
                            ## consecutive values below that value
                            # idx_best <- as.numeric(which.max(
                            # unlist(coarse_prev_df["q2_vals"])))
                            # threshold <- coarse_prev_df[idx_best, "q2"] -
                            #   coarse_prev_df[idx_best, "SE"]
                            #
                            ## TODO
                            fgIL <- max(mi[kCGIdx]-5, 1)
                            ## shield against setting 0
                            fgOL <- max(lo[kCGIdx]-1, 1)
                            .msg_range(fgIL, fgOL, vrbs)
                        }else{
                            ## Added to handle case when 1-SE rule is not
                            ## applied
                            ##
                            fgIL <- max(mi[kCGIdx]-1, 1)
                            ## shield against setting 0
                            fgOL <- max(mi[kCGIdx]-1, 1)
                            .msg_range(fgIL, fgOL, vrbs)
                        }
                        break
                        ##
                    } else if(best_K == lo[kCGIdx]){
                        if (best_K == min(lo)){ ## fgIL is 1
                            .msg_pstr("min(lo) is best", flg=dbg)
                            fgIL <- 1
                            fgOL <- min(lo)-1
                            .msg_range(fgIL, fgOL, vrbs)
                            break
                        } else{
                            ## go fine over interval (hi[kCGIdx]+1 , lo[kCGIdx])
                            .msg_pstr("lo value is best", flg=dbg)
                            fgOL <- lo[kCGIdx]-1    # fine-grained search outer
                            fgIL <- hi[kCGIdx-1]+1  # and inner limit
                            .msg_range(fgIL, fgOL, vrbs)
                            break
                        }
                        ##
                    } else if(kCGIdx > 1 && best_K == hi[kCGIdx-1]){
                        .msg_pstr("prev hi is best", flg=dbg)
                        fgIL <- hi[kCGIdx-1]+1
                        fgOL <- lo[kCGIdx]-1
                        .msg_range(fgIL, fgOL, vrbs)
                        break
                    }
                    else if(best_K == hi[kCGIdx]){
                        ## best_K is == hi, go to next coarse-grained iteration
                        .msg_pstr("Next interval of coarse-grained grid",
                                    flg=vrbs)

                    }
                }
                ## Fine-grained search
                searchReturnFine <- performSearchForK( X = X,
                    cvfolds = cvfolds,
                    startVal = fgIL,
                    endVal = fgOL,
                    step = 1,
                    prev_best_K = -1,
                    best_K = 0,
                    prev_df = NULL,
                    param_ranges, parallelDo = parallelDo,
                    kFolds = kFolds, nRuns = nRuns,
                    set_verbose = set_verbose)
                best_K <- searchReturnFine$best_K
                fine_prev_df <- searchReturnFine$return_df

                combined_df <- rbind(prev_df, fine_prev_df)
                best_K <- .get_best_K(combined_df, parsimony = askParsimony)
                ## Ensure chosen value is not the lower boundary
                # minKInDF <- min(as.numeric(unlist(combined_df["k_vals"])))
                kValsInDF <- as.numeric(unlist(combined_df["k_vals"]))
                if(best_K != 1 && !any((best_K-1) - kValsInDF == 0)){
                    .msg_pstr("Chosen best K is lowest in the range tested.",
                        "Making sure...", flg=vrbs)
                    attemptCount <- 1
                    makeSureK <- best_K
                    while(makeSureK != 1 &&
                            !any((makeSureK-1) - kValsInDF == 0)){
                        # Ensure the new value is not already computed
                        fgIL <- max(best_K - 1*attemptCount, 1)
                        fgOL <- fgIL
                        .msg_range(fgIL, fgOL, vrbs)
                        searchReturnFine <-
                            performSearchForK( X = X,
                                cvfolds = cvfolds,
                                startVal = fgIL,
                                endVal = fgOL,
                                step = 1,
                                prev_best_K = -1,
                                best_K = 0,
                                prev_df = combined_df,
                                param_ranges = param_ranges,
                                parallelDo = parallelDo,
                                kFolds = kFolds,
                                nRuns = nRuns,
                                set_verbose = set_verbose)
                        # temp_best_K <- searchReturnFine$best_K
                        combined_df <- searchReturnFine$return_df

                        # combined_df <- rbind(combined_df, fine_prev_df)
                        makeSureK <- .get_best_K(combined_df,
                                                    parsimony = askParsimony)

                        attemptCount <- attemptCount + 1
                        # minKInDF <- min(as.numeric(
                        #               unlist(combined_df["k_vals"])))
                        # message("MIN_K_IN_DF: ", minKInDF)
                        kValsInDF <- as.numeric(unlist(combined_df["k_vals"]))
                    }
                    best_K <- makeSureK
                    .msg_pstr("Made sure, best K is: ", best_K, flg=dbg)
                }

            }
        }
    # }
    if(best_K == max(param_ranges$k_vals)){
        cli::cli_alert_warning(c("Best K: {best_K} is already the maximum ",
                                "value for K specified. Try increasing ",
                                "the range"))
    }
    ##
    return(best_K)
}
## =============================================================================

.msg_range <- function(fgIL, fgOL, vrbs){
    .msg_pstr("Choosing interval:", fgIL, "-", fgOL, flg=vrbs)
}


.setup_par_cluster <- function(vlist){
    ## from perform_multiple_NMF_runs
    cl <- parallel::getDefaultCluster()
    parallel::clusterEvalQ(cl, suppressWarnings(require(MASS, quietly = TRUE)))
    parallel::clusterExport(cl = NULL, varlist = vlist,
        envir = parent.frame(n=1))
    # ## ^for pseudo-inverse using function `ginv` (CV-based model selection)

    return(cl)
}


check_par_conditions <- function(nCores){
    if (is.na(nCores)) {
        ## raise error or handle
        stop("'parallelize' is TRUE, but 'nCores' not specified")
    } else {
        if (nCores <= parallel::detectCores()) {
            ##
        } else {
            stop("Specified more than available cores. Stopping")
        }
    }
}

# return A list of two lists: one containing feature matrices, the other samples
# matrices
.perform_multiple_NMF_runs <- function(X, kVal, alphaVal, parallelDo = TRUE,
                                nCores = NA, nRuns = 100, bootstrap = TRUE) {
    ##
    new_ord <- lapply(seq_len(nRuns), function(x){seq_len(ncol(X))})
    if(bootstrap){
        new_ord <- lapply(seq_len(nRuns), function(x){
            sample(ncol(X), ncol(X), replace = FALSE)})
    }
    ##
    seed_val_list <- get_n_seeds(n = nRuns)
    ## In parallel
    if (parallelDo) {
        #TODO: ?
        check_par_conditions(nCores=nCores)
        ##
        cl <- .setup_par_cluster(vlist=c(".perform_single_NMF_run"))
        ##
        nmf_result_list <- parallel::clusterApplyLB(cl = cl, seq_len(nRuns),
                            function(i) {
                                .perform_single_NMF_run(
                                X = X[,new_ord[[i]]], kVal = kVal,
                                alphaVal = alphaVal, seedVal = seed_val_list[i])
                            })
        return(list(nmf_result_list = nmf_result_list, new_ord = new_ord))
        ##
    } else{
        ## In serial
        nmf_result_list <- lapply(seq_len(nRuns),
                            function(i) {
                                .perform_single_NMF_run(
                                X = X[,new_ord[[i]]], kVal = kVal,
                                alphaVal = alphaVal, seedVal = seed_val_list[i])
                            })
        return(list(nmf_result_list = nmf_result_list, new_ord = new_ord))
        ##
    }
}



# @title Generate Cross-Validation Data Splits
#
# @description This function generates the row and column indices for the
# cross-validation splits.
#
# @param Xdims Dimensions of matrix data matrix X.
# @param kFolds Number of cross-validation folds.
#
# @return A list of two elements: rowIDs and columnIDs for different cross-
# validation folds.
# @importFrom cvTools cvFolds
#
.generate_folds <- function(Xdims, kFolds) {
    ##
    ## Xdims gives the dimensions of the matrix X
    cvf_rows <- cvTools::cvFolds(Xdims[1], K = kFolds, type = "random")

    cvf_cols <- cvTools::cvFolds(Xdims[2], K = kFolds, type = "random")

    return(list(cvf_rows = cvf_rows, cvf_cols = cvf_cols))
}
## =============================================================================

# @title Get the best performing value of K (number of factors in NMF)
#
# @param x A data.frame/tibble, grid_search_results, as returned by
# \code{.cv_model_select_pyNMF2}
#
# @return A number The best performing value of K.
.get_best_K <- function(x, parsimony = FALSE) {
    # Assumes, max q2_val is best
    # Returns simply the best performing K value
    check_names_params(x)
    ##
    averages <-
        .get_q2_aggregates_chosen_var(x, chosen_var = x$k_vals, base::mean)
    ######
    ## Using the one-std error rule for selecting a parsimonious model
    sd_by_K <- .get_q2_aggregates_chosen_var(x, x$k_vals, stats::sd)
    se_by_K <- sd_by_K / sqrt(nrow(x)/nrow(sd_by_K))
    ##
    averages$SD <- sd_by_K$q2_vals
    averages$SE <- se_by_K$q2_vals
    ##
    ######
    idx_best <- as.numeric(which.max(unlist(averages["q2_vals"])))
    ##
    best_K <- as.numeric(averages[idx_best, "rel_var"])
    if(parsimony){
        ## Return the value of K using the one SE rule
        se_rule_threshold <- averages[idx_best, "q2_vals"] -
            averages[idx_best, "SE"]
        best_K_by_se_rule <- averages[
            which(unlist(averages["q2_vals"]) > se_rule_threshold), "rel_var"]
        best_K_by_se_rule_chosen <- min(best_K_by_se_rule)
        return(best_K_by_se_rule_chosen)
    }
    return(best_K)
}
## =============================================================================

check_names_params <- function(x){
    if (length(setdiff(
        names(x),
        c("k_vals", "alpha", "fold", "iteration",
            "q2_vals")
    )) > 0) {
        .msg_pstr(names(x), flg=FALSE)
        stop(
            "Check colnames in data.frame, expecting five element names: ",
            c("k_vals", "alpha", "fold", "iteration", "q2_vals")
        )
    }
}


# @title Aggregate \eqn{Q^2} Values
#
# @description Aggregate the \eqn{Q^2} values from the grid search results.
#
# @param x The return object from \code{\link{.cv_model_select_pyNMF2}}.
# @param chosen_var The variable to aggregate over.
# @param chosen_func The aggregate function to use (should be a function
# already existing wihtin R). Possible values are: \code{mean} and \code{sd}.
#
# @return The mean of $Q^2$ values per the chosen variable
# @importFrom stats aggregate
#
.get_q2_aggregates_chosen_var <- function(x, chosen_var, chosen_func) {
    ## Returns the mean of q2 values per the chosen variable
    averages <- stats::aggregate(x,
                                by = list(rel_var = chosen_var),
                                chosen_func,
                                simplify = TRUE)
    return(averages)
}
## =============================================================================


# @title Get Threshold Value for Selecting \eqn{\alpha}
#
# @description Get the threshold value for selection of \eqn{\alpha} by
# looking at cross-validation performance for K.
#
# @param model_selectK Cross-validation performance over K values.
#
# @return The \eqn{Q^2} threshold value.
# @importFrom stats sd
.get_q2_threshold_by_K <- function(model_selectK) {
    ##
    mean_by_K <- .get_q2_aggregates_chosen_var(model_selectK,
                                model_selectK$k_vals,
                                base::mean)
    sd_by_K <- .get_q2_aggregates_chosen_var(model_selectK,
                                model_selectK$k_vals,
                                stats::sd)
    se_by_K <- sd_by_K / sqrt(nrow(sd_by_K))
    ##
    best_K <- .get_best_K(model_selectK)
    idx_best_K <- which(sd_by_K$rel_var == best_K)
    ##
    q2_threshold <- list(mean = as.numeric(mean_by_K[idx_best_K, "q2_vals"]),
                        se = as.numeric(se_by_K[idx_best_K, "q2_vals"]))
    ##
    return(q2_threshold)
}
## =============================================================================


# @title Plot Cross-Validation Performance of K
#
# @description Plot showing performance of different values of K tested in
# cross-validation.
#
# @param averages The average of performance values for different combinations
# in grid search
#
# @return A ggplot object so you can simply call \code{print} or \code{save}
# on it later.
#
# @import ggplot2
.plot_cv_K <- function(averages) {
    # Using ggplot to plot
    if ("rel_var" %in% averages) {
        cat("Plotting Q2 vs. K: check colnames in the object
        (need 'rel_var' as 'k')")
        return(NULL)
    }
    if ("q2_vals" %in% averages) {
        cat("Plotting Q2 vs. K: check colnames in the object (need 'q2_vals')")
        return(NULL)
    }
    cat("Plotting Q2 as a function of K")
    p1 <- ggplot2::ggplot(averages, aes_string(x = "rel_var", y = "q2_vals")) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = averages$rel_var,
                            labels = averages$rel_var) +
        ggplot2::labs(
            title = "Reconstruction accuracy, Q\U00B2 = f(#Factors)",
            x = "#Factors (K)",
            y = paste0("Reconstruction accuracy (Q\U00B2)")
        )
    return(p1)
}
## =============================================================================
