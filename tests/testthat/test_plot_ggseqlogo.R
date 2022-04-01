context("Plotting seqlogo w/ ggseqlogo")
library(ggseqlogo)


test_that("Matrix has 4 rows", {
  testPwmMat <- matrix(rnorm(100), nrow = 2)
  testPositionLabels <- seq(25)
  expect_error(
    suppressWarnings(viz_pwm(testPwmMat, method = "bits",
                              pos_lab = testPositionLabels)),
    paste("Expecting a matrix with 4 rows corresponding to DNA chars",
      "'A', 'C', 'G', 'T'")
  )
})

test_that("Given object is matrix", {
  testPwmMat <- rnorm(200) # err
  testPositionLabels <- seq(25)
  expect_error(
    suppressWarnings(viz_pwm(testPwmMat, method = "bits",
                              pos_lab = testPositionLabels)),
    "Expecting a matrix with 4 rows"
  )
})

test_that("Handling empty matrix", {
  testPwmMat <- matrix()
  testPositionLabels <- seq(25)
  expect_error(
    suppressWarnings(viz_pwm(testPwmMat, method = "bits",
                              pos_lab = testPositionLabels)),
    "Empty"
  )
})

test_that("Position labels inadequate", {
  testPwmMat <- matrix(rnorm(100), nrow = 4)
  testPositionLabels <- seq(20)
  expect_error(
    suppressWarnings(viz_pwm(testPwmMat, method = "bits",
                              pos_lab = testPositionLabels)),
    "Inadequate"
  )
})

test_that("Position labels over-abundant", {
  testPwmMat <- matrix(rnorm(100), nrow = 4)
  testPositionLabels <- seq(50)
  expect_error(
    suppressWarnings(viz_pwm(testPwmMat, method = "bits",
                              pos_lab = testPositionLabels)),
    "Overabundant"
  )
})

test_that("ggseqlogo plotting works", {
  # setting seed enables proper comparison between ggplot objects since we use
  # rnorm
  set.seed(1234)
  # test variables
  # testPositionLabels <- seq(25)
  # testPwmMat <- matrix(rnorm(100), nrow = 4)
  # testPwmMat <- make_PWMs(testPwmMat)
  # p1 <- viz_pwm(testPwmMat, pos_lab = testPositionLabels)
  #
  res <- readRDS(system.file("extdata", "example_seqArchRresult.rds",
          package = "seqArchR", mustWork = TRUE))

  pwm <- seqArchR::make_PWMs(get_clBasVec_m(res,iter=1)[,1],
                         add_pseudo_counts = FALSE, sinuc = FALSE)

  p1 <- suppressWarnings(viz_pwm(pwm_mat = pwm, method = "bits",
                                  fixed_coord = TRUE))
  # test plot
  vdiffr::expect_doppelganger("ggseqlogo plot example", p1)
})
