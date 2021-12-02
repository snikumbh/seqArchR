context("prepare_data_from_FASTA")

test_that("File non-existance handled", {
  testSeqs_filenameNotOK <- system.file("extdata", "example_data_.fa")
  # This returns empty string when file not found

  expect_error(prepare_data_from_FASTA(testSeqs_filenameNotOK), "not found")
})

test_that("Preparing sinuc/dinuc/trinuc data matrix", {
    testSeqs_filename <- system.file("extdata", "example_data.fa",
                                     package = "seqArchR",
                                     mustWork = TRUE)
    # This returns empty string when file not found
    rawSeqs <- prepare_data_from_FASTA(testSeqs_filename, raw_seq = TRUE)
    use_rawSeqs <- base::substr(rawSeqs[1:2], start = 1, stop = 5)
    # 20 x 2 sparse Matrix of class "dgCMatrix"
    #
    # [1,] . .
    # [2,] . .
    # [3,] 1 .
    # [4,] . .
    # [5,] . 1
    # [6,] . .
    # [7,] . 1
    # [8,] . .
    # [9,] . .
    # [10,] . .
    # [11,] 1 .
    # [12,] . .
    # [13,] . 1
    # [14,] . .
    # [15,] . .
    # [16,] . 1
    # [17,] 1 .
    # [18,] . .
    # [19,] 1 1
    # [20,] 1 .
    #
    example_sinuc <- matrix(rep(0, 40), nrow = 20, ncol = 2)
    example_sinuc[c(3,11,17,19,20),1] <- example_sinuc[c(5,7,13,16,19),2] <- 1
    example_sinuc <- Matrix::Matrix(example_sinuc, sparse = TRUE)
    rownames(example_sinuc) <-
          .get_feat_names(k=1, seqlen=Biostrings::width(use_rawSeqs[1]))

    expect_identical(get_one_hot_encoded_seqs(use_rawSeqs,
                                             sinuc_or_dinuc = "sinuc"),
                     example_sinuc)

    #
    # 80 x 2 sparse Matrix of class "dgCMatrix"
    #
    # [1,] . .
    # [2,] . .
    # [3,] . .
    # [4,] . .
    # [5,] . .
    # [6,] . .
    # [7,] . .
    # [8,] . .
    # [9,] . .
    # [10,] . .
    # [11,] . .
    # [12,] . .
    # [13,] . .
    # [14,] . .
    # [15,] . .
    # [16,] . .
    # [17,] 1 .
    # [18,] . .
    # [19,] . 1
    # [20,] . .
    # [21,] . .
    # [22,] . .
    # [23,] . .
    # [24,] . .
    # [25,] . .
    # [26,] . .
    # [27,] . .
    # [28,] . .
    # [29,] . .
    # [30,] . .
    # [31,] . .
    # [32,] . .
    # [33,] . .
    # [34,] . .
    # [35,] . .
    # [36,] . 1
    # [37,] . .
    # [38,] . .
    # [39,] . .
    # [40,] . .
    # [41,] . .
    # [42,] . .
    # [43,] . .
    # [44,] . .
    # [45,] . .
    # [46,] . .
    # [47,] . 1
    # [48,] . .
    # [49,] . .
    # [50,] . .
    # [51,] . .
    # [52,] . .
    # [53,] . .
    # [54,] . .
    # [55,] . .
    # [56,] . .
    # [57,] . .
    # [58,] . .
    # [59,] . .
    # [60,] . .
    # [61,] . .
    # [62,] . .
    # [63,] 1 .
    # [64,] . .
    # [65,] . .
    # [66,] . .
    # [67,] . .
    # [68,] . .
    # [69,] . .
    # [70,] . .
    # [71,] 1 .
    # [72,] . .
    # [73,] . 1
    # [74,] . .
    # [75,] . .
    # [76,] . .
    # [77,] . .
    # [78,] . .
    # [79,] 1 .
    # [80,] . .
    example_dinuc <- matrix(rep(0, 160), nrow = 80, ncol = 2)
    example_dinuc[c(17,63,71,79),1] <- example_dinuc[c(19,36,47,73),2] <- 1
    example_dinuc <- Matrix::Matrix(example_dinuc, sparse = TRUE)
    rownames(example_dinuc) <-
      .get_feat_names(k=2, seqlen=Biostrings::width(use_rawSeqs[1]))

    expect_identical(get_one_hot_encoded_seqs(use_rawSeqs,
                                             sinuc_or_dinuc = "dinuc"),
                     example_dinuc)

    example_trinuc <- matrix(rep(0, 640), nrow = 320, ncol = 2)
    example_trinuc[c(303,257,71),1] <- example_trinuc[c(287,196,73),2] <- 1
    example_trinuc <- Matrix::Matrix(example_trinuc, sparse = TRUE)
    rownames(example_trinuc) <-
      .get_feat_names(k=3, seqlen=Biostrings::width(use_rawSeqs[1]))

    expect_identical(get_one_hot_encoded_seqs(use_rawSeqs,
                                              sinuc_or_dinuc = "trinuc"),
                     example_trinuc)
})
