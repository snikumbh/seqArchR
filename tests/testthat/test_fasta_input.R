context("FASTA input")

test_that("All same length", {
  testSeqsNotOK <- c(
    "AACGTGACTAT",
    "ACCGATCGAT",
    "GGCATCATGC",
    "TCATCTAGAT"
  )
  testSeqsOK <- c(
    "AACGTGACTA",
    "ACCGATCGAT",
    "GGCATCATGC",
    "TCATCTAGAT"
  )
  expect_error(.assert_seq_attributes(testSeqsNotOK), "different")
  expect_message(.assert_seq_attributes(testSeqsOK), "OK")
})

test_that("None with length zero", {
  testSeqsNotOK1 <- c(
    "",
    "ACCGATCGAT",
    "GGCATCATGC",
    "TCATCTAGAT"
  )
  testSeqsNotOK2 <- c(
    "",
    "",
    "",
    ""
  )
  expect_error(.assert_seq_attributes(testSeqsNotOK1), "zero")
  expect_error(.assert_seq_attributes(testSeqsNotOK2), "zero")
})


test_that("Only characters from DNA alphabet", {
  testSeqsNotOK1 <- c(
    "AACGTGNNTA",
    "ACRGATCGAT",
    "GGCATCATGC",
    "TCATCTAGAT"
  )
  expect_warning(.assert_seq_attributes(testSeqsNotOK1), "Non DNA-alphabet character")
})
