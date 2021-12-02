context("Test heatmap and seqlogo")
library(ggplot2)
library(ggseqlogo)

test_that("Given object is matrix", {
  testFeaturesMat <- rnorm(10000) # err
  testPositionLabels <- seq(25)
  expect_error(
    viz_bas_vec_heatmap_seqlogo(
      testFeaturesMat,
      pos_lab = testPositionLabels),
    "not of type matrix"
  )
})

test_that("Handling empty matrix", {
  testFeaturesMat <- matrix()
  testPositionLabels <- seq(25)
  expect_error(
    viz_bas_vec_heatmap_seqlogo(
      testFeaturesMat,
      pos_lab = testPositionLabels
    ),
    "Empty"
  )
})

test_that("Position labels inadequate", {
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200)
  testPositionLabels <- seq(20)
  expect_error(
    viz_bas_vec_heatmap_seqlogo(
      testFeaturesMat,
      pos_lab = testPositionLabels),
    "Inadequate"
  )
})

test_that("Position labels over-abundant", {
  testFeaturesMat <- matrix(rnorm(10000), nrow = 200)
  testPositionLabels <- seq(60)
  expect_error(
    viz_bas_vec_heatmap_seqlogo(
      testFeaturesMat,
      pos_lab = testPositionLabels),
    "Overabundant"
  )
})

## Commented after gridExtra changed to ggpubr
# test_that("Combined heatmap and seqlogo works", {
#   # setting seed enables proper comparison between ggplot objects since we use
#   # rnorm
#   set.seed(11223344)
#   # test variables
#   testPositionLabels <- seq(25)
#   testFeaturesMat <- matrix(rnorm(100), ncol = 1)
#   # test plot, include the function directly inside because it returns nothing.
#   vdiffr::expect_doppelganger(
#     "heatmap seqlogo plot ex",
#     viz_bas_vec_heatmap_seqlogo(
#       featuresMatrix = testFeaturesMat,
#       pos_lab = testPositionLabels)
#   )
# })
