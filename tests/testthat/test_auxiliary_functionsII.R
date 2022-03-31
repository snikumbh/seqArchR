context("Auxiliary Functions II")

test_that("setup_clustFactors return value", {
    toyMat <- matrix(rep(runif(1),1000), ncol = 5)
    returnVal <- .setup_clustFactors_for_seqArchR_result(toyMat)

    expect_identical(names(returnVal), c("nBasisVectors", "basisVectors"))
    expect_equal(returnVal$nBasisVectors, 5)
    expect_equal(returnVal$basisVectors, toyMat)
})

test_that("decide_process_outer_chunk works fine", {

    expect_message(.decide_process_outer_chunk(25, 24, 5),
                    "Sorry, will not process this small a chunk: 24")
    expect_true(.decide_process_outer_chunk(25, 24, 5))
    expect_false(.decide_process_outer_chunk(25, 30, 5))
    expect_error(.decide_process_outer_chunk(15, 30, 4),
                    "'min_size' should be at least 4 times 'kFolds'")
    expect_error(.decide_process_outer_chunk(25, 0, 5),
                    "lengthOfOC > 0 is not TRUE")
})


test_that("Null dir path is flagged in allSeqs Logo", {
  # fname <- system.file("extdata", "example_data.fa",
  #   package = "seqArchR",
  #   mustWork = TRUE)
  # tssSeqsRaw <- suppressMessages(seqArchR::prepare_data_from_FASTA(fname,
  #     raw_seq = TRUE))

  tssSeqsRaw <- readRDS(system.file("extdata", "tssSeqsRaw.rds",
                                    package = "seqArchR",
                                    mustWork = TRUE))

  expect_error(plot_all_seqs_logo(NULL, NULL, dpath = NULL),
              "seqs_raw is NULL")
  expect_error(plot_all_seqs_logo(seqs_raw = tssSeqsRaw[1:50],
                                  NULL, dpath = NULL),
              "directory path/name is NULL")
})




test_that("get_dimers is working fine", {
  alphabet = c("A", "C", "G", "T")
  expectAns <- c("AA", "CA", "GA", "TA", "AC", "CC", "GC", "TC", "AG", "CG",
                 "GG", "TG", "AT", "CT", "GT", "TT")
  ans <- get_dimers_from_alphabet(alphabet)
  expect_equal(expectAns, ans)
  expect_error(get_dimers_from_alphabet(NULL))
})


test_that("get_trimers is working fine", {
  alphabet = c("A", "C", "G", "T")
  expectAns <- c("AAA", "CAA", "GAA", "TAA", "ACA", "CCA", "GCA", "TCA", "AGA",
  "CGA", "GGA", "TGA", "ATA", "CTA", "GTA", "TTA", "AAC", "CAC", "GAC", "TAC",
  "ACC", "CCC", "GCC", "TCC", "AGC", "CGC", "GGC", "TGC", "ATC", "CTC", "GTC",
  "TTC", "AAG", "CAG", "GAG", "TAG", "ACG", "CCG", "GCG", "TCG", "AGG", "CGG",
  "GGG", "TGG", "ATG", "CTG", "GTG", "TTG", "AAT", "CAT", "GAT", "TAT", "ACT",
  "CCT", "GCT", "TCT", "AGT", "CGT", "GGT", "TGT", "ATT", "CTT", "GTT", "TTT")
  ans <- get_trimers_from_alphabet(alphabet)
  expect_equal(expectAns, ans)
  expect_error(get_trimers_from_alphabet(NULL))
})


# test_that("get_hopach_cluster_medoidsIdx handles null hopach object", {
#     hopachObj <- NULL
#     fMat <- matrix(rep(runif(1),1000), ncol = 5)
#     expect_error(.get_hopach_cluster_medoidsIdx(hopachObj),
#                  "Hopach object is NULL")
#     ## Make hopach object to test
#
# })


test_that("NMF result from py script OK", {
    nmfList <- vector("list", 2)
    samplesMat <- matrix(rep(runif(1),1000), ncol = 5)
    featuresMat <- matrix(rep(runif(1),1000), nrow = 5)
    expect_error(get_features_matrix(nmfList), "0LengthEntry")
    nmfList[[1]] <- featuresMat
    expect_error(get_samples_matrix(nmfList), "0LengthEntry")
    nmfList[[2]] <- samplesMat
    expect_equal(get_samples_matrix(nmfList), samplesMat)
    expect_equal(get_features_matrix(nmfList), featuresMat)
})


test_that("getting dimers works", {
    expect_error(get_dimers_from_alphabet(NULL))
    expect_equal(get_dimers_from_alphabet(c(1,0)), c("11", "01", "10", "00"))
})




test_that("get_factors_from_factor_clustering handles null hopach object", {
    hopachObj <- NULL
    fMat <- matrix(rep(runif(1),1000), ncol = 5)
    expect_equal(.get_factors_from_factor_clustering2(hopachObj, fMat),
                 fMat)
})


test_that("get_factors_from_factor_clustering handles all-zero
          featuresMatrix/factors", {
    hopachObj <- list(clustering = list(k = 5, sizes = rep(5,5), order = 1:25))
    fMat <- matrix(rep(0,1000), ncol = 5)
    expect_error(.get_factors_from_factor_clustering2(hopachObj, fMat),
                 "All zeroes as factors")
})


test_that("get_factors_from_factor_clustering handles NA in
          featuresMatrix/factors", {
    hopachObj <- list(clustering = list(k = 5, sizes = rep(5,5), order = 1:25))
    fMat <- matrix(rep(0,1000), ncol = 5)
    fMat[2,2] <- NA
    expect_error(.get_factors_from_factor_clustering2(hopachObj, fMat),
                 "Matrix has NA values")
})



