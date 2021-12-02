context("Plotting features matrix as ggheatmap")
library(ggplot2)
library(reshape2)


test_that("Matrix has 4 rows", {
  testPwmMat <- matrix(rnorm(100), nrow = 2)
  testPositionLabels <- seq(25)
  expect_error(
    plot_ggheatmap(testPwmMat, pos_lab = testPositionLabels),
    paste0("Expecting a matrix with 4 rows corresponding to DNA chars ",
    "'A', 'C', 'G', 'T'")
  )
})

test_that("Given object is matrix", {
  testPwmMat <- rnorm(200) # err
  testPositionLabels <- seq(25)
  expect_error(
    plot_ggheatmap(testPwmMat, pos_lab = testPositionLabels),
    "Expecting a matrix with 4 rows"
  )
})

test_that("Handling empty matrix", {
  testPwmMat <- matrix()
  testPositionLabels <- seq(25)
  expect_error(
    plot_ggheatmap(testPwmMat, pos_lab = testPositionLabels),
    "Empty"
  )
})

test_that("Position labels inadequate", {
  testPwmMat <- matrix(rnorm(100), nrow = 4)
  testPositionLabels <- seq(20)
  expect_error(
    plot_ggheatmap(testPwmMat, pos_lab = testPositionLabels),
    "Inadequate"
  )
})

test_that("Position labels over-abundant", {
  testPwmMat <- matrix(rnorm(100), nrow = 4)
  testPositionLabels <- seq(50)
  expect_error(
    plot_ggheatmap(testPwmMat, pos_lab = testPositionLabels),
    "Overabundant"
  )
})

test_that("ggheatmap plotting works", {
  # setting seed enables proper comparison between ggplot objects since we use
  # rnorm
  set.seed(11223344)
  # test variables
  testPositionLabels <- seq(5)
  testPwmMat <- matrix(rnorm(20), nrow = 4)
  p1 <- plot_ggheatmap(testPwmMat, pos_lab = testPositionLabels)
  # test plot
  vdiffr::expect_doppelganger("ggheatmap plot example", p1)
})
