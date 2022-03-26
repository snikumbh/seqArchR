context("seqArchR main functionality")
library(seqArchR)

tssSeqs_sinuc <- readRDS(system.file("extdata", "tssSinuc.rds",
                                     package = "seqArchR",
                                     mustWork = TRUE))
tssSeqsRaw <- readRDS(system.file("extdata", "tssSeqsRaw.rds",
                                  package = "seqArchR",
                                  mustWork = TRUE))

## Reduce size to avoid timing out on Bioconductor
tssSeqs_sinuc <- tssSeqs_sinuc[,seq(100)]
tssSeqsRaw <- tssSeqsRaw[seq(100)]


test_that("seqArchR (stability) works when timeFlag is FALSE/checkpointing", {
    ## Make toy objects and data
    ## Commented out. We could just save the sinuc and dinuc profiles instead
    ## of having to read and recompute them every time.
    ##
    # fname <- system.file("extdata", "example_data.fa.gz",
    #                                   package = "seqArchR",
    #                                   mustWork = TRUE)
    #
    #
    #
    # tssSeqs_sinuc <- suppressMessages(seqArchR::prepare_data_from_FASTA(fname))
    # tssSeqsRaw <- suppressMessages(seqArchR::prepare_data_from_FASTA(fname,
    #                                 raw_seq = TRUE))



    nSeqs <- ncol(tssSeqs_sinuc)
    positions <- seq(1,100)

    useFlags <- list(debug = FALSE,
                     verbose = TRUE,
                     plot = FALSE,
                     time = FALSE)
    toyConfig <- seqArchR::set_config(chunk_size = 100,
                                        k_min = 1, k_max = 20,
                                        parallelize = FALSE,
                                        mod_sel_type = "stability",
                                        n_runs = 50,
                                        n_cores = NA, checkpointing = TRUE,
                                        flags = useFlags)
    set.seed(1234)
    seqArchRresult <- seqArchR::seqArchR(toyConfig,
                                seqs_raw = tssSeqsRaw,
                                seqs_ohe_mat = tssSeqs_sinuc,
                                total_itr = 2, set_ocollation = c(TRUE,FALSE),
                                o_dir = "temp_test")
    ##
    expect_equal_to_reference(seqArchRresult,
                              "seqArchRresult_stability_check_timeFalse.rds")

    temp <- list.dirs(".")
    temp <- unlist(lapply(temp, basename))
    temp_l <- lapply(temp, function(x){unlist(strsplit(x, split="_"))})
    recent <- suppressWarnings(unlist(lapply(temp_l, function(x){
        foo <- tail(x, 1)
        foo <- ifelse(length(foo) == 1, as.numeric(foo), NA)
        foo
        })))
    recent <- tail(recent, 1)
    recent <- ifelse(is.na(recent), "temp_test", paste0("temp_test_", recent))
    rds_fname <- file.path(recent, "seqArchRresult_checkpoint1.rds")
    chkpoint <- readRDS(rds_fname)
    unlink(recent, recursive = TRUE)
    expect_equal_to_reference(chkpoint, "seqArchRresult_checkpoint1.rds")

    ##
})


test_that("seqArchR (stability) works when plot==TRUE, o_dir is NULL", {
    ## Make toy objects and data
    # fname <- system.file("extdata", "example_data.fa.gz",
    #     package = "seqArchR",
    #     mustWork = TRUE)
    #
    #
    #
    # tssSeqs_sinuc <- suppressMessages(seqArchR::prepare_data_from_FASTA(fname))
    # tssSeqsRaw <- suppressMessages(seqArchR::prepare_data_from_FASTA(fname,
    #     raw_seq = TRUE))

    nSeqs <- ncol(tssSeqs_sinuc)
    positions <- seq(1,100)

    useFlags <- list(debug = FALSE,
        verbose = TRUE,
        plot = TRUE,
        time = FALSE)
    toyConfig <- seqArchR::set_config(chunk_size = 100,
        k_min = 2, k_max = 20,
        parallelize = FALSE,
        mod_sel_type = "stability",
        n_runs = 50,
        n_cores = NA,
        flags = useFlags)
    set.seed(1234)
    ##
    expect_error(seqArchR::seqArchR(toyConfig,
        seqs_raw = tssSeqsRaw,
        seqs_ohe_mat = tssSeqs_sinuc,
        total_itr = 1, o_dir = NULL),
        ##
        paste("'plot' flag is TRUE but 'o_dir' is not provided.",
            "Did you forget to set 'o_dir'?")
        )
    ##
})




test_that("seqArchR (cv) works when timeFlag is FALSE", {
    ## Make toy objects and data
    # fname <- system.file("extdata", "example_data.fa.gz",
    #                                 package = "seqArchR",
    #                                 mustWork = TRUE)
    #
    #
    #
    # tssSeqs_sinuc <-
    #     suppressMessages(seqArchR::prepare_data_from_FASTA(fname))
    # tssSeqsRaw <-
    #     suppressMessages(seqArchR::prepare_data_from_FASTA(fname,
    #                         raw_seq = TRUE))

    nSeqs <- ncol(tssSeqs_sinuc)
    positions <- seq(1,100)

    useFlags <- list(debug = TRUE,
                     verbose = TRUE,
                     plot = FALSE,
                     time = FALSE)
    toyConfig <- seqArchR::set_config(chunk_size = 100,
                                   k_min = 2, k_max = 20, parallelize = TRUE,
                                   mod_sel_type = "cv",
                                   n_runs = 10,
                                   n_cores = 2,
                                   flags = useFlags)
    ## Test cross-validation-based model selection.
    ## This needs to parallel as TRUE.
    # skip_on_travis()
    set.seed(1234)
    seqArchRresult <- #suppressMessages(
        seqArchR::seqArchR(toyConfig, seqs_raw = tssSeqsRaw,
                                seqs_ohe_mat = tssSeqs_sinuc,
                                total_itr = 1, set_ocollation = TRUE)
        # )
    expect_equal_to_reference(seqArchRresult, "seqArchRresult_cv_check_timeFalse.rds")
})



test_that("seqArchR (stability) works when debug & timeFlag is FALSE", {
    ## Make toy objects and data
    # fname <- system.file("extdata", "example_data.fa.gz",
    #                                   package = "seqArchR",
    #                                   mustWork = TRUE)
    #
    #
    #
    # tssSeqs_sinuc <- suppressMessages(seqArchR::prepare_data_from_FASTA(fname))
    # tssSeqsRaw <- suppressMessages(seqArchR::prepare_data_from_FASTA(fname,
    #                             raw_seq = TRUE))

    nSeqs <- ncol(tssSeqs_sinuc)
    positions <- seq(1,100)
    ## keeping debug as TRUE doesn't alter the seqArchR result object which we
    ## check, but it helps hit many more debug message lines in the code
    useFlags <- list(debug = TRUE,
                     verbose = TRUE,
                     plot = FALSE,
                     time = FALSE)
    toyConfig <- seqArchR::set_config(chunk_size = 100,
                                   k_min = 1, k_max = 20, parallelize = FALSE,
                                   mod_sel_type = "stability",
                                   n_runs = 50,
                                   n_cores = NA,
                                   flags = useFlags)
    set.seed(1234)
    seqArchRresult <- #suppressMessages(
        seqArchR::seqArchR(toyConfig, seqs_raw = tssSeqsRaw,
                    seqs_ohe_mat = tssSeqs_sinuc, total_itr = 1,
                    set_ocollation = TRUE)
        #)
    ##
    expect_equal_to_reference(seqArchRresult,
                              "seqArchRresult_stability_check_debugTrue.rds")
    ##
})


test_that("Handles negative threshold iteration", {
    useFlags <- list(debug = TRUE,
        verbose = TRUE,
        plot = FALSE,
        time = FALSE)
    # toyResult <- seqArchR(toyConfig, seqs_ohe_mat = tssSeqs, total_itr = -1)
    expect_error(.assert_seqArchR_thresholdIteration(-1),
                 "Expecting number of iterations to be numeric and > 0")
    toyConfig <- seqArchR::set_config(chunk_size = 100,
        k_min = 2, k_max = 20, parallelize = FALSE,
        mod_sel_type = "stability",
        n_runs = 50,
        n_cores = NA,
        flags = useFlags)
    expect_error(seqArchR(toyConfig, seqs_raw = tssSeqsRaw,
        seqs_ohe_mat = tssSeqs_sinuc, total_itr = -1, set_ocollation = FALSE),
                 "Expecting number of iterations to be numeric and > 0")
})

# test_that("Handles negative threshold iteration from seqArchR main", {
#     expect_error(seqArchR(toyConfig, seqs_raw = tssSeqsRaw,
#     seqs_ohe_mat = tssSeqs_sinuc, total_itr = -1))
#
#     expect_error(seqArchR(toyConfig, seqs_raw = tssSeqsRaw,
#     seqs_ohe_mat = tssSeqs_sinuc, total_itr = -1),
#                  "Expecting threshold iteration to be numeric and > 0")
# })


test_that("Config handles: negative chunk_size", {
    expect_error(set_config(chunk_size = -500,
                                k_min = 1, k_max = 8, parallelize = TRUE,
                                cv_folds = 3, n_runs = 50,
                                n_cores = 2),
                "'chunk_size' should be > 0")
})

# This test is now null and void because if the flags variable is NULL,
# we handle it by assigning the default flag values in set_config func itself
# test_that("Config handles: NULL flags", {
#     expect_error(set_config(chunk_size = 500,
#                                 k_min = 2, k_max = 8, parallelize = TRUE,
#                                 cv_folds = 3, n_runs = 50,
#                                 n_cores = 2,
#                                 flags = NULL),
#                  "'flags' are NULL")
# })

# test_that("Config handles: improper (names in) flags", {
#     expect_error(set_config(chunk_size = 500,
#                                 k_min = 2, k_max = 8, parallelize = TRUE,
#                                 cv_folds = 3, n_runs = 50,
#                                 n_cores = 2,
#                                 flags = list(debug = FALSE,
#                                             plot = FALSE,
#                                             verboseFl = TRUE,
#                                             time = TRUE)),
#                  "Unexpected names or no elements in flags variable")
# })

test_that("Config handles: improper (non-logical) flags", {
    expect_error(set_config(chunk_size = 500,
                                k_min = 2, k_max = 8, parallelize = TRUE,
                                cv_folds = 3, n_runs = 50,
                                n_cores = 2,
                                flags = list(debug = NULL,
                                            plot = FALSE,
                                            verbose = TRUE,
                                            time = TRUE)),
            "Expected only LOGICAL values in flag variable, found otherwise")
})


test_that("Config handles: negative n_runs", {
    expect_error(set_config(chunk_size = 500,
                                k_min = 2, k_max = 8, parallelize = TRUE,
                                cv_folds = 3, n_runs = -50,
                                n_cores = 2,
                                flags = list(debug = FALSE,
                                            plot = FALSE,
                                            verbose = TRUE,
                                            time = TRUE)),
                 "'n_runs' should be > 0")
})



test_that("Config handles: negative min_size", {
    expect_error(set_config(chunk_size = 500,
                                k_min = 2, k_max = 8, parallelize = TRUE,
                                cv_folds = 3, n_runs = 50,
                                n_cores = 2,
                                min_size = -20,
                                flags = list(debug = FALSE,
                                                plot = FALSE,
                                                verbose = TRUE,
                                                time = TRUE)),
                    "'min_size' should be > 0")
})


test_that("Config handles: negative alphaVal", {
    expect_error(set_config(chunk_size = 500,
                                k_min = 2, k_max = 8, parallelize = TRUE,
                                cv_folds = 3, n_runs = 50,
                                n_cores = 2,
                                min_size = 20,
                                alpha_base = -3,
                                alpha_pow = 3,
                                flags = list(debug = FALSE,
                                                plot = FALSE,
                                                verbose = TRUE,
                                                timeFlag = TRUE)),
                    paste("Resulting alpha value is < 0.",
                            "Check 'alpha_base' and 'alpha_pow'"))
})


test_that("Config handles: chunk_size < nSeqs", {
    expect_error(.assert_seqArchR_chunkSize_in_tandem(500,450),
                    "'chunk_size' should be <= number of input sequences")
})

test_that("Config handles: kFolds NULL", {
    expect_error(.assert_seqArchR_kFolds_in_tandem(NULL,450),
                    "'kFolds' is NULL")
})

test_that("Config handles: kFolds is not numeric", {
    expect_error(.assert_seqArchR_kFolds_in_tandem("5",450),
                 "'kFolds' should be numeric and > 0")
})

test_that("Config handles: kFolds <= nSeqs", {
    expect_error(.assert_seqArchR_kFolds_in_tandem(451,450),
                 paste0("CV folds should be less than or equal to #sequences. ",
                 "Standard values: 3, 5, 10."))
})



test_that("next iteration chunks list has a 0-length entry", {
    toyList <- vector("list", 5)
    toyList <- lapply(seq_along(toyList),
                        function(x){
                            if ( x != 3) {
                                toyList[[x]] <- rep(seq(1:25),5)
                            }
                        })

    expect_error(.assert_seqArchR_OK_for_nextIteration(toyList),
                    "Index 3 of zero length")
})

test_that("next iteration chunks list is NULL", {
    toyList <- NULL
    expect_error(.assert_seqArchR_OK_for_nextIteration(toyList),
                 "Chunks for next iteration are NULL")
})
