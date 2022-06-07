test_that("shuffle_df shuffles an afq dataframe", {
  subjectID <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
  nodeID <- c(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3)
  dti_fa <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2)
  age <- c(22, 22, 22, 22, 33, 33, 33, 33, 44, 44, 44, 44)
  sex <- c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0)
  group <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)

  df_afq <- data.frame(subjectID,
                       nodeID,
                       dti_fa,
                       age,
                       sex,
                       group)

  set.seed(123)
  df_shuffle <- shuffle_df(df_afq, "dti_fa")

  df_shuffle$group <- as.numeric(df_shuffle$group)
  df_shuffle$nodeID <- as.numeric(df_shuffle$nodeID)

  expect_setequal(unique(df_shuffle$sex), c(0, 1))
  expect_setequal(df_afq$age, df_shuffle$age)
  expect_setequal(df_afq$group, df_shuffle$group)
  expect_false(isTRUE(all.equal(df_afq$age, df_shuffle$age)))
  expect_false(isTRUE(all.equal(df_afq$group, df_shuffle$group)))
  expect_false(isTRUE(all.equal(df_afq$sex, df_shuffle$sex)))
  expect_equal(df_afq$nodeID, df_shuffle$nodeID)
  expect_equal(df_afq$subjectID, df_shuffle$subjectID)
  expect_equal(df_afq$dti_fa, df_shuffle$dti_fa)
})
