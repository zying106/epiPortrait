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

test_that("find_hockey_inflection rejects inverse left-loaded curves", {
  inverse <- c(1, seq(90, 100, length.out = 99))
  res <- find_hockey_inflection(inverse, method = "elbow")
  expect_equal(res$call_status, "no_call")
  expect_lt(res$quality_score, 1e-8)

  right_tail <- c(rep(1, 90), seq(2, 100, length.out = 10))
  expect_equal(find_hockey_inflection(right_tail, method = "elbow")$call_status,
               "called")
})

test_that("find_hockey_inflection validates min_quality", {
  expect_error(find_hockey_inflection(1:10, min_quality = NA_real_),
               "min_quality")
  expect_error(find_hockey_inflection(1:10, min_quality = 1.1),
               "min_quality")
})
