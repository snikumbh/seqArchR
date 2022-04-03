## Function to check that the flags variable
## Expected to be:
## 1. not NULL
## 2. having the right fields and
## 3. all are LOGICAL values
.assert_seqArchR_flags <- function(flags) {
    ## Expected fields
    expNames <-  c("debugFlag", "verboseFlag", "plotVerboseFlag", "timeFlag")
    if (is.null(flags)) {
        stop("'flags' are NULL")
    }
    if (!is.list(flags)) {
        stop("flags variable is not a list")
    } else {
        matchNames <- names(flags) %in% expNames
        if (!is.null(names(flags)) && !all(matchNames)) {
            stop("Unexpected names or no elements in flags variable")
        } else {
            ## test all are set to logical values
            all_logical <- unlist(lapply(flags, is.logical))
            if (!all(all_logical)) {
                stop("Expected only LOGICAL values in flag variable, ",
                        "found otherwise")
            }
        }
    }
}
## =============================================================================


## Function to check properties of a features or a samples matrix from NMF
## This is combined function instead of two separate functions to check them
## separately.
## Expected to be:
## - not null
## - a matrix
## - ncol > 0 or nrow > 0
## - for featuresMat, if all matrix is 0, using a hopach can throw error,
##  check and safegaurd against it

.assert_seqArchR_featSampMatrix <- function(fs_mat, feat = TRUE){
    ##
    if (is.null(fs_mat)) {
        stop("NULL value found, instead of a matrix")
    }
    if (!is.matrix(fs_mat)) {
        stop("Expected a matrix, found otherwise")
    } else {
        if (any(is.na(fs_mat))) {
            stop("Matrix has NA values")
        }
        if (!all(dim(fs_mat) > 0))
            stop("Matrix has zero rows or columns")
        if (feat && all(fs_mat == 0)) {
            ## This will lead to an error if hopach is performed,
            ## hence throwing error
            stop("All zeroes as factors")
        }
    }
}
## =============================================================================



## Function to independently check validity of min_size param in config
## Expected to be:
## 1. not NULL
## 2. numeric and > 0
.assert_seqArchR_min_size_independent <- function(minSeqs_var) {
    .check_null_num(var = minSeqs_var, use_name = 'min_size')
    .check_gtZero(var = minSeqs_var, use_name = 'min_size')
}
## =============================================================================



## Function to check validity of minSeqs w.r.t. given sequences
## Expected to be:
## 1. independently valid (see f .assert_seqArchR_minSeqs_independent)
## 2. < #given sequences
##
.assert_seqArchR_min_size_in_tandem <- function(minSeqs_var, given_seqs_size) {
    .assert_seqArchR_min_size_independent(minSeqs_var)
    if (minSeqs_var > given_seqs_size) {
        stop("'min_size' is > number of input sequences.
            Typically, a small number")
    }
}
## =============================================================================



.check_null_num <- function(var, use_name){
    if (is.null(var)) {
        stop(use_name, " is NULL")
    }
    if (!is.numeric(var)) {
        stop(use_name, " should be numeric")
    }
}
## =============================================================================

.check_gtZero <- function(var, use_name){
    if (var < 0) {
        stop(use_name, " should be > 0")
    }
}
## =============================================================================




## Function to check validity of kFolds param in config
## Expected to be:
## 1. not NULL
## 2. numeric and > 0
##
.assert_seqArchR_kFolds_in_tandem <- function(kFolds_var,
                                            given_seqs_size = NULL) {
    ##
    .check_null_num(var = kFolds_var, use_name = 'kFolds')
    if (kFolds_var < 3) {
        stop("Set at least 3 cross-validation folds")
    }
    if (!is.null(given_seqs_size) && kFolds_var > given_seqs_size) {
        ## using leave-one-out cv will kill, so not suggested actually
        stop("CV folds should be less than or equal to #sequences. ",
            "Standard values: 3, 5, 10.")
    }

}
## =============================================================================


## Function to check validity of parallelization-related param in config.
## parallelize and nCoresUse
## parallelize is expected to be:
## 1. not NULL
## 2. LOGICAL
## If parallelize is TRUE, nCoresUse is expected to be:
## 1. not NULL
## 2. numeric and > 0
.assert_seqArchR_parallelization_params <- function(par_var, nCores_var) {
    if (is.null(par_var)) {
        stop("'parallelize' is NULL, expected LOGICAL")
    }
    if (!is.logical(par_var)) {
        stop("'parallelize' should be LOGICAL (TRUE/FALSE)")
    }
    ## nCoresUse, bother about it only when parallelize is TRUE
    if (par_var) {
        .check_null_num(var = nCores_var, use_name = 'nCoresUse')
        .check_gtZero(var = nCores_var, use_name = 'nCoresUse')
        if (nCores_var > parallel::detectCores()) {
            stop("Specified more than available cores. Available cores: ",
                parallel::detectCores())
        }
    }
}
## =============================================================================



## Function to check properties of collection of chunks for next iteration
## Expected to be
## 1. not NULL
## 2. A list
## 3. all elements of non-zero/positive length
##
.assert_seqArchR_OK_for_nextIteration <- function(nxtOuterChunksColl) {
    ## Check if it is NULL
    if (is.null(nxtOuterChunksColl)) {
        stop("Chunks for next iteration are NULL")
    }
    if (!is.list(nxtOuterChunksColl)) {
        stop("Chunks for next iteration expected as a list, found otherwise")
    }
    ## Check if any chunk in the collection is of zero-length?
    if (any(lapply(nxtOuterChunksColl, length) == 0)) {
        stop("Chunks for next iteration have a problem. ",
        "Index ", which(lapply(nxtOuterChunksColl, length) == 0),
                " of zero length")
    }
}
## =============================================================================


## Function to check properties of thresholdItr set by user
## Expected to be
## 1. not NULL
## 2. A numeric
## 3. > 0
##
.assert_seqArchR_thresholdIteration <- function(given_val) {
    .check_null_num(var = given_val, use_name = "Threshold iteration")
    .check_gtZero(var = given_val, use_name = "Threshold iteration")
}
## =============================================================================

## Function to check properties of paramRanges used for model selection w/ NMF
## Expected to be
## 1. not NULL
## 2. List with elements having a fixed set of names
## 3. acceptable range for alpha: any thing non-negative
## 4. acceptable range for K_vals: any thing positive (> 0)
##
.assert_seqArchR_model_selection_params <- function(param_ranges_var) {
    expNames <- c("alphaBase", "alphaPow", "k_vals")
    alphaBase <- param_ranges_var$alphaBase
    alphaPow <- param_ranges_var$alphaPow
    k_vals <- param_ranges_var$k_vals
    if (is.null(param_ranges_var)) {
        stop("'param_ranges' is NULL")
    }
    matchNames <- names(param_ranges_var) %in% expNames
    if (!all(matchNames)) {
        stop("Unexpected names of 'param_ranges' elements. Expected:", expNames)
    } else {
        if (is.numeric(k_vals) || is.numeric(alphaBase) ||
            is.numeric(alphaPow)) {
            ## all OK, check ranges now
            alphaVal <- alphaBase^alphaPow
            if (alphaVal < 0) {
                stop("Resulting alpha value is < 0. ",
                    "Check 'alpha_base' and 'alpha_pow'")
            }
            if (!all(k_vals > 0)) {
                stop("'k_vals' should be > 0")
            }
        } else {
            stop("Either of 'k_vals', 'alphaBase' or 'alphaPow' ",
                "is not a numeric value")
        }
    }
}
## =============================================================================



## Function to check validity of nRunsUse for NMF in config
## Expected to be:
## 1. not NULL
## 2. numeric and > 0
## 3. ?
.assert_seqArchR_nRuns <- function(nIter_var) {
    .check_null_num(var = nIter_var, use_name = 'n_runs')
    .check_gtZero(var = nIter_var, use_name = 'n_runs')
}
## =============================================================================




## Function to independently check validity of chunkSize
## Expected to be:
## 1. not NULL
## 2. numeric and >0
## When everything is satisfied, still whether it is valid depends on what is
## the total number of sequences
.assert_seqArchR_chunkSize_independent <- function(chunkSize_var) {
    .check_null_num(var = chunkSize_var, use_name = 'chunk_size')
    .check_gtZero(var = chunkSize_var, use_name = 'chunk_size')
}
## =============================================================================



## Function to check validity of chunkSize w.r.t. given sequences
## Expected to be:
## 1. independently valid (see f .assert_seqArchR_chunkSize_independent)
## 2. < #given sequences
##
.assert_seqArchR_chunkSize_in_tandem <- function(chunkSize_var,
                                                given_seqs_size) {
    .assert_seqArchR_chunkSize_independent(chunkSize_var)
    if (chunkSize_var > given_seqs_size) {
        stop("chunk_size should be <= number of input sequences")
    }
}
## =============================================================================



## Function to check properties of configuration variables
## Expected to be:
## 1. not NULL
## 1. A list with fixed set of names
## 2. Individual element assertions should pass
## 3.
.assert_seqArchR_config <- function(config_var, seqs_size = NA) {
    ##TODO: Check assertions for new args modSelType, tol and bound
    expNames <- c("modSelType", "bound", "kFolds",
                    "parallelize", "nCoresUse", "nRunsUse",
                    "paramRanges", "chunkSize", "minSeqs",
                    "flags", "checkpointing", "result_aggl", "result_dist")
    if (is.null(config_var)) {
        stop("'config' is NULL")
    }
    if (!is.list(config_var)) {
        stop("'config' is expected to be a list, found otherwise")
    } else {
        matchNames <- names(config_var) %in% expNames
        if (!all(matchNames)) {
            stop("Check names of elements in 'config'")
        }
        .assert_seqArchR_kFolds_in_tandem(config_var$kFolds, NULL)
        .assert_seqArchR_chunkSize_independent(config_var$chunkSize)
        if (!is.na(seqs_size)) {
            .assert_seqArchR_kFolds_in_tandem(config_var$kFolds, seqs_size)
            .assert_seqArchR_chunkSize_in_tandem(config_var$chunkSize,
                                                    seqs_size)
            .assert_seqArchR_min_size_in_tandem(config_var$minSeqs, seqs_size)
        }
        .assert_seqArchR_min_size_independent(config_var$minSeqs)
        .assert_seqArchR_model_selection_params(config_var$paramRanges)
        .assert_seqArchR_parallelization_params(config_var$parallelize,
                                            config_var$nCoresUse)
        .assert_seqArchR_nRuns(config_var$nRunsUse)
        .assert_seqArchR_flags(config_var$flags)
    }

}
## =============================================================================


.assert_seqArchR_collation <- function(set_ocollation, total_itr){
    if(is.null(set_ocollation)){
        stop("Please specify an outer chunk collation strategy. Found NULL")
    }
    if(length(set_ocollation) < total_itr){
        stop("Expecting length of set_ocollation to be same as total_itr. ",
                " Found ", total_itr, " and ", length(set_ocollation),
                call. = TRUE)
    }
}



## We need to write functions to sanity check the important variables
## throughout the seqArchR algorithm. These are:
## -- globFactors Matrix
## -- globFactorsClustering Hopach Object
## -- globClustAssignments List
## -- outerChunk Matrix?
## -- outerChunksColl List
## -- innerChunksColl List
## -- intClustFactors Matrix
## -- clustFactors Matrix
## -- seqsClustLabels Vector
## -- collatedClustAssignments List
## --
##
## =============================================================================


## Function to check validity of NMFresult from .handle_chunk_w_NMF function
## Expected to be:
## 1. not NULL
## 2. a list w elements 'forGlobFactors' (a matrix) & 'forGlobClustAssignments'
## (a list)
## 3.
.assert_seqArchR_NMFresult <- function(NMFresultObj) {
    if (is.null(NMFresultObj)) {
        stop("NMF result object is NULL")
    }
    expNames <- c("forGlobFactors", "forGlobClustAssignments")
    matchNames <- names(NMFresultObj) %in% expNames
    if (!all(matchNames)) {
        stop("Check names of elements in NMF result for inner chunk")
    }
    if (!is.list(NMFresultObj$forGlobFactors) &&
        !is.list(NMFresultObj$forGlobClustAssignments)) {
        stop("NMF result should hold information as list w/ elements per
                inner chunk")
    } else {
        .assert_seqArchR_featSampMatrix(NMFresultObj$forGlobFactors,
                                        feat = TRUE)
        if (!is.list(NMFresultObj$forGlobClustAssignments)) {
            stop("'Global cluster assignments' in NMF result should be a list")
        } else {
            ## Check number of clusters == #factors (ncol(featuresMat))
            if (length(NMFresultObj$forGlobClustAssignments) !=
                ncol(NMFresultObj$forGlobFactors)) {
                    stop("In NMF result, #clusters != #factors")
            }
        }
    }
}
## =============================================================================



## =============================================================================

## Function to check validity of seqsClustLabels
## Expected to be:
## 1. not NULL
## 2. non-empty vector
## 3. all entries in the vector of same lengths
.assert_seqArchR_seqsClustLabels <- function(seqsClustLabels) {
    if (is.null(seqsClustLabels)) {
        stop("Cluster labels for sequences NULL")
    }
    if (is.vector(seqsClustLabels) && length(seqsClustLabels) == 0) {
        stop("Expecting cluster labels for sequences as a non-empty vector")
    }
}

## Function to check validity of seqsClustLabels at the end of an iteration
## Expected to be:
## 1. general assertions passing
## 2. all entries in the vector of same lengths
.assert_seqArchR_seqsClustLabels_at_end <- function(seqsClustLabels) {
    splitChar <- "-"
    one <- 1
    .assert_seqArchR_seqsClustLabels(seqsClustLabels)
    all_lengths <- unlist(lapply(strsplit(seqsClustLabels,
                                            split = splitChar),
                                    length))
    check_length <- length(unique(all_lengths))
    if (check_length != one) {
        stop("Cluster labels for sequences have different lengths")
    }
}
## =============================================================================

## Function to check validity of seqArchRresult object
## Expected to be a nested list of various elements as follows:
## 1. seqsClustLabels and clustBasisVectors: list with as many elements as
## number of iterations
##   - each seqsClustLabels: a character vector of labels for all sequences per
##   iteration
##   - each clustBasisVectors: a list with elements 'nBasisVectors' and
##   'basisVectors'
##     - nBasisVectors: integer, number of basis vectors for that iteration
##     - basisVectors: double matrix, basisVectors, dimensions along rows,
##     different basis vectors along columns
## 2. clustSol: list storing the final clustering solution of seqArchR; this has
## three elements
##   - 'basisVectorsClust': a list storing clustered basis vectors from
##   last iteration of seqArchR (result from reorder_seqArchRresult)
##   - 'clusters': a list storing relevant sequence indices per cluster. Length
##   of 'clusters' and 'basisVectorsClust' should be same.
##   - 'seqsClustLabels': labels per sequence according to clustering in
##   clusters
## 3. rawSeqs: The DNAStringSet object storing all sequences provided to
## seqArchR. These sequences have names in the DNAStringSet object.
## 4. timeInfo: a list storing time required (in minutes) for each iteration
## of seqArchR
## 5. config: configuration set for processing the corresponding data with
## seqArchR
## 6. call: the function call itself.
##
##
.assert_seqArchRresult <- function(seqArchRresultObj) {
    ## seqArchRresult object always has same list elements.
    ## When time info is not requested, it is set to NA
    topLevel_elem_names <- c("seqsClustLabels", "clustBasisVectors",
        "clustSol", "rawSeqs", "timeInfo", "config", "call")

    if (is.null(seqArchRresultObj)) {
        stop("seqArchR result object is NULL")
    }

    if(!identical(names(seqArchRresultObj), topLevel_elem_names)){
        stop("seqArchR result object list elements' name do not match")
    }
    eq_itr_lengths <- length(seqArchRresultObj$seqsClustLabels)
    ##
    if(!all(length(seqArchRresultObj$clustBasisVectors) == eq_itr_lengths
            &&
            ifelse(seqArchRresultObj$config$flags$timeFlag,
                length(seqArchRresultObj$timeInfo) == eq_itr_lengths, TRUE)
            )){
        stop("seqArchR result object has erroneous iter info for labels,
                basis vectors and/or timeInfo")
    }
    ##
    eq_seqs_length <- length(seqArchRresultObj$rawSeqs)
    logi_seqs_length <- unlist(lapply(seqArchRresultObj$seqsClustLabels,
        function(x){
            length(x) == eq_seqs_length
        }))
    if(!all(logi_seqs_length)) stop("Sequence clust labels for iteration(s) of
                                inappropriate length")

    if(!(eq_seqs_length == length(seqArchRresultObj$clustSol$seqsClustLabels)))
        stop("Sequence labels for final clustering solution of inappropriate
                length")
    ##
    ## assert if dimensions of clustBasisVectors are fine.
    ## Sequences can be encoded as monomers, dimers, etc.
    ## this should be 4L, (4^2)L, (4^3)L, and so on
    factor_rows <- unlist(lapply(seqArchRresultObj$clustBasisVectors,
        function(x){
            tmp <- Biostrings::width(seqArchRresultObj$rawSeqs[1])
            log(nrow(x$basisVectors)/(tmp), base=4)
    }))
    if(!all(diff(factor_rows) == 0))
        stop("seqArchR result clustBasisVectors have unequal dimensions across
                iterations")

    factor_cols <- unlist(lapply(seqArchRresultObj$clustBasisVectors,
        function(x){
            x$nBasisVectors > 0 &&
            is(x$basisVectors, "matrix") &&
            x$nBasisVectors == ncol(x$basisVectors)
        }))
    if(!all(factor_cols))
        stop("nCols in basisVectors do not match nBasisVectors")

    .assert_seqArchR_seqsClustLabels(seqArchRresultObj$seqsClustLabels)

    ## assert attributes in clustSol
    ## - element names
    ## - lengths of clusters and basisVectorClust
    clustSol_exp_names <- names(seqArchRresultObj$clustSol)
    if(!identical(clustSol_exp_names, names(seqArchRresultObj$clustSol)))
        stop("Element names in clustSol of result object do not match")

    clustSol <- seqArchRresultObj$clustSol
    if(length(clustSol$basisVectorsClust) != length(clustSol$clusters))
        stop("Mismatch in number of clusters in clustSol$basisVectorsClust
                and clustSol$clusters")


}
## =============================================================================


## Function to check validity of lists in general
## Expected to be:
## 1. not NULL
## 2. a list
## 3. Not have any 0-length element
##
.assert_seqArchR_list_properties <- function(listVar) {
    returnMessage <- "FOO"
    if (is.null(listVar)) {
        returnMessage <- "NULL"
    }
    if (!is.list(listVar)) {
        returnMessage <- "Nlist"
    } else {
        if(length(listVar) == 0){
            returnMessage <- "0LengthList"
        }else{
            check_lengths <- lapply(listVar, length)
            if (any(check_lengths == 0)) {
                returnMessage <- "0LengthEntry"
            }
        }
    }
    return(returnMessage)
}
## =============================================================================



.assert_seqArchR_globClustAssignments <- function(given_var) {
    returnMessage <- .assert_seqArchR_list_properties(given_var)
    if (returnMessage == "NULL") {
        stop("Cluster assignments variable is NULL")
    }
    if (returnMessage == "Nlist") {
        stop("Cluster assignments variable is not a list")
    }
    if (returnMessage == "0LengthEntry") {
        stop("Cluster assignments variable has a 0-length entry")
    }
}
## =============================================================================


## Function to check consistency of nSeqs in clustLabels variable and in
## clustAssignments variable
## Expected to be:
## 1. Ensure list variable is OK and seqsClustLabels is OK
## 2. holding same number of sequences (the variable lengths)
##
.assert_seqArchR_consistent_nSeqs_w_clusters <- function(seqsClustLabels,
                                                    clustAssignments) {

    .assert_seqArchR_seqsClustLabels(seqsClustLabels)
    .assert_seqArchR_globClustAssignments(clustAssignments)
    nSeqs_in_labels <- length(seqsClustLabels)
    nSeqs_in_assignments <- length(unlist(clustAssignments))
    if (!nSeqs_in_labels == nSeqs_in_assignments) {
        stop("Number of sequences in seqsClustLabels and clustAssignments not",
                " equal")
    }
}
## =============================================================================
