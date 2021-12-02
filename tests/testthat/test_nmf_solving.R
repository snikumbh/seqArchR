context("Solving NMF and reconstruction")


test_that("Q^2 computation: Original matrix empty", {
  matA <- matrix() # rnorm(20), nrow = 4)
  reMatA <- matrix(rnorm(20), nrow = 4)
  expect_error(.compute_q2(matA, reMatA), "Empty")
})

test_that("Q^2 computation: Reconstructed matrix empty", {
  matA <- matrix(rnorm(20), nrow = 4) # rnorm(20), nrow = 4)
  reMatA <- matrix()
  expect_error(.compute_q2(matA, reMatA), "Empty")
})

test_that("Q^2 computation: handling original non-matrix", {
  matA <- c(rnorm(20)) # matrix(rnorm(20), nrow = 4)#rnorm(20), nrow = 4)
  reMatA <- matrix(rnorm(20), nrow = 4)
  expect_error(.compute_q2(matA, reMatA), "not of type matrix")
})

test_that("Q^2 computation: handling reconstructed non-matrix", {
  matA <- matrix(rnorm(20), nrow = 4)
  reMatA <- c(rnorm(20)) # matrix()
  expect_error(.compute_q2(matA, reMatA), "not of type matrix")
})

# test_that("NMF solved correctly", {
#
#
#
# })
