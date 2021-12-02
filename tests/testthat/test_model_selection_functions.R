context("Model selection functions")


test_that("Handles < 3 cross-validation folds", {
  testDataX <- matrix(rnorm(10000), nrow = 200)
  testParamRanges <- list(alphaBase = c(), alphaPow = c(), k_vals = c())
  testKFolds <- 1
  expect_error(
    .cv_model_select_pyNMF2(testDataX,
      param_ranges = testParamRanges,
      kFolds = testKFolds
    ),
    "Set at least 3 cross-validation folds"
  )
})

test_that("Handles negative cross-validation folds", {
  testDataX <- matrix(rnorm(10000), nrow = 200)
  testParamRanges <- list(alphaBase = c(), alphaPow = c(), k_vals = c())
  testKFolds <- -4
  expect_error(
    .cv_model_select_pyNMF2(testDataX,
      param_ranges = testParamRanges,
      kFolds = testKFolds
    ),
    "cannot be negative"
  )
})

test_that("Handles CV folds > #sequences", {
  testDataX <- matrix(rnorm(10000), nrow = 200)
  testParamRanges <- list(alphaBase = c(), alphaPow = c(), k_vals = c())
  testKFolds <- ncol(testDataX) + 5
  expect_error(
    .cv_model_select_pyNMF2(testDataX,
      param_ranges = testParamRanges,
      kFolds = testKFolds
    ),
    "CV folds should be less than or equal to #sequences"
  )
})

test_that("Handles non-matrix type X", {
  testDataX <- rnorm(10000) # err
  testParamRanges <- list(alphaBase = c(), alphaPow = c(), k_vals = c())
  testKFolds <- 5
  expect_error(
    .cv_model_select_pyNMF2(testDataX,
      param_ranges = testParamRanges,
      kFolds = testKFolds
    ),
    "not of type matrix"
  )
})

test_that("Handles erroneous names in param_ranges list", {
  testDataX <- matrix(rnorm(10000), nrow = 200)
  testParamRanges <- list(alphabase = 0, alphaPow = 0, k_vals = 0)
  testKFolds <- 5
  expect_error(
    .cv_model_select_pyNMF2(testDataX,
      param_ranges = testParamRanges,
      kFolds = testKFolds
    ),
    "Expected elements in param ranges: alphaBase, alphaPow, k_vals"
  )
})

test_that("Returned names is 'generate_folds' return object are right", {
  testDataX <- matrix(rnorm(10000), nrow = 200)
  testNames <- list(cvf_rows = c(), cvf_cols = c())
  testKFolds <- 5
  givenAns <- .generate_folds(dim(testDataX), kFolds = testKFolds)
  expect_equal(names(givenAns), names(testNames))
})


test_that("Handles erroneous colnames in input data.frame (get_best_K)", {
  testDataX <- matrix(rnorm(10000), nrow = 200)
  testParamRanges <- list(alphabase = 0, alphaPow = 0, k_vals = 0)
  testKFolds <- 5
  testModelSelectResult <- data.frame(
    k_vals = seq(4), alpha = rnorm(4),
    fold = seq(4), q2vals = rnorm(4)
  )
  expect_error(
    .get_best_K(testModelSelectResult),
    "Check colnames in data.frame"
  )
})
