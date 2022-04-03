## NMF model selection functions II
## Other approaches for model selection
##
## Note on regularization:
## We use dispersion for regularization together to Amari-type distance measure.
## The latter is a measure of stability of features, and the first of the
## clusters.
##

.stability_model_select_pyNMF2 <- function(X,
                            param_ranges,
                            nRuns = 100,
                            bootstrap = TRUE,
                            returnBestK = TRUE,
                            bound = 10^-6,
                            flags = list(debugFlag = FALSE,
                                verboseFlag = TRUE, plotVerboseFlag = FALSE,
                                timeFlag = FALSE),
                            bpparam
                            ){
    dbg <- flags$debugFlag
    vrbs <- flags$verboseFlag
    ##
    foo_scores <- expand.grid(list(kValue = param_ranges$k_vals,
                                    ScoreType = c("AmariTypeDistance"),
                                    nRuns = nRuns, Score = 100
                            ))
    ##
    .msg_pstr("Bound : ", bound, flg=dbg)
    ##
    bestK <- 1
    breakNow <- FALSE
    prev_amari <- NA
    for(kValue in param_ranges$k_vals){
        if(kValue > 1){
            .msg_pstr("Checking K = ", kValue, flg=vrbs)
        }

        ## Get NMF results
        resultList <- .perform_multiple_NMF_runs(X = X, kVal = kValue,
            alphaVal = 0, nRuns = nRuns,
            bootstrap = bootstrap, bpparam = bpparam)

        featMatList <- .get_feat_or_samp_matList(resultList, feat = TRUE)

        this_amari <- .get_amari_from_featMatList(featMatList)
        if(is.na(this_amari)){
            cli::cli_alert_warning("NA in Amari-type distance computation")
        }
        .msg_pstr("AmariTypeDist : ", this_amari, flg=dbg)

        foo_scores[kValue, "Score"] <- this_amari

        ##
        check_errant <- FALSE
        if(kValue %% 5 == 0){
            check_errant <- .check_no_mag_change_fail_condition(foo_scores)
            if(check_errant) breakNow <- TRUE
        }
        if(is.na(this_amari)){
            cli::cli_alert_info("errant TRUE")
            check_errant <- TRUE
            breakNow <- TRUE
        }
        ##
        if(!is.na(this_amari) && this_amari > bound) breakNow <- TRUE
        ##
        if(breakNow){
            ## greater than bound, choose and break loop
            ## Here, if kValue = 1, bestK would be assigned 0. Avoid this.
            if(kValue > 1) bestK <- kValue - 1
            if(check_errant) bestK <- 1
            break
        }
        ##
        prev_amari <- this_amari

    }
    if(returnBestK) return(bestK)
    return(foo_scores)
}


.check_no_mag_change_fail_condition <- function(foo_scores){
    score_pows <- abs(floor(log10(foo_scores$Score)))
    score_pows <- score_pows[which(!is.nan(score_pows))]
    if(all(score_pows >= 17)){
        return(TRUE)
    }
    return(FALSE)
}

# .mag_change <- function(A, B){
#     return(abs(floor(log10(A)) - floor(log10(B))))
# }


.get_feat_or_samp_matList <- function(resultList, bootstrap = TRUE,
                                        feat = TRUE){
    if(feat){
        featMatList <- lapply(resultList$nmf_result_list, get_features_matrix)
        return(featMatList)
    }
    sampMatList <- lapply(resultList$nmf_result_list, get_samples_matrix)

    if(bootstrap){
        sampMatListNew <-
            lapply(seq_len(length(sampMatList)), function(x){
                thisMat <- as.matrix(sampMatList[[x]])
                thisNR <- nrow(thisMat)
                thisNC <- ncol(thisMat)
                origOrdX <- matrix(rep(-100,thisNR*thisNC),
                    nrow = thisNR, ncol =  thisNC)
                origOrdX[,resultList$new_ord[[x]]] <- thisMat
                origOrdX
            })
        sampMatList <- sampMatListNew
    }
    return(sampMatList)
}


.get_amari_from_featMatList <- function(featMatList){
    return(computeAmariDistances(featMatList))
}

amariDistance <- function(matA, matB) {
    K <- dim(matA)[2]
    corrMat <- stats::cor(matA, matB)
    return(1 - (sum(apply(corrMat, 1, max)) +
            sum(apply(corrMat, 2, max))) / (2 * K))
}


computeAmariDistances <- function(matrices){
    B <- length(matrices)
    distances.list <- unlist(lapply(seq_len(B - 1), function(b) {
        distances <- lapply(seq(from=b + 1, to=B, by=1), function(b.hat) {
            amariDistance(matrices[[b]], matrices[[b.hat]])
        })
    })
    )
    return(mean(distances.list))
}
