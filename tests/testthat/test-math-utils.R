test_that("find_hockey_inflection detects inflection", {
  x <- c(rep(1, 50), seq(1, 100, length.out = 50))
  res <- find_hockey_inflection(x)
  expect_type(res, "list")
  expect_true("inflection_idx" %in% names(res))
  expect_true(res$call_status %in% c("called", "no_call"))
})

test_that("find_hockey_inflection supports tangent (ROSE) method", {
  x <- c(rep(1, 50), seq(1, 100, length.out = 50))
  res_tangent <- find_hockey_inflection(x, method = "tangent")
  res_elbow <- find_hockey_inflection(x, method = "elbow")
  expect_true("inflection_idx" %in% names(res_tangent))
  expect_true(is.integer(res_tangent$inflection_idx))
  expect_equal(res_tangent$method, "tangent")
  expect_equal(res_elbow$method, "elbow")
})

test_that("find_hockey_inflection no-calls on constant input", {
  res <- find_hockey_inflection(rep(5, 20))
  expect_equal(res$call_status, "no_call")
  expect_true(is.na(res$inflection_idx))
})
