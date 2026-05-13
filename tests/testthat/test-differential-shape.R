context("Differential Shape Testing")

test_that("test_global_shape_shift validates inputs", {
  data(example_se)
  expect_error(test_global_shape_shift(example_se, group_var = "NotThere"),
               "Group variable")
})

test_that("test_global_shape_shift works with formula interface", {
  data(example_se)
  res <- test_global_shape_shift(example_se,
                                 formula = ~ Condition,
                                 workers = 1)
  expect_true(all(c("Shape_Shift_Score", "partial_eta_sq", "P.Value", "adj.P.Val") %in% colnames(res)))
  expect_equal(nrow(res), nrow(example_se))
})

test_that("test_global_shape_shift works with group_var (backward-compat)", {
  data(example_se)
  res <- test_global_shape_shift(example_se,
                                 group_var = "Condition",
                                 workers = 1)
  expect_true("partial_eta_sq" %in% colnames(res))
  expect_true(all(res$partial_eta_sq >= 0 & res$partial_eta_sq <= 1, na.rm = TRUE))
})

test_that("test_differential_feature supports all methods", {
  data(example_se)
  # limma
  res_l <- test_differential_feature(example_se, feature = "Intensity",
    group_var = "Condition", target_group = "Treatment", ref_group = "Control",
    method = "limma")
  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(res_l)))

  # t.test
  res_t <- test_differential_feature(example_se, feature = "Intensity",
    group_var = "Condition", target_group = "Treatment", ref_group = "Control",
    method = "t.test", workers = 1)
  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(res_t)))

  # wilcox
  res_w <- test_differential_feature(example_se, feature = "Intensity",
    group_var = "Condition", target_group = "Treatment", ref_group = "Control",
    method = "wilcox", workers = 1)
  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(res_w)))

  expect_error(test_differential_feature(example_se, feature = "NotThere",
    group_var = "Condition", target_group = "Treatment", ref_group = "Control"))
})
