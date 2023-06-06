# testthat::test_that("Valid arguments", {
#   expect_equal(LimmaWeight(LimmaFit)$Log2FC, LimmaOut$Log2FC, tolerance = 0.01)
# })

testthat::test_that("Valid arguments", {
  expect_equal(setNames(LimmaOut$Signed.Q.Value,row.names(LimmaOut)), weight, tolerance = 0.01)
})






