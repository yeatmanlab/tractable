test_that("load.afq.sarica loads the sarica dataset", {
  df_sarica <- load.afq.sarica()

  expect_equal(nrow(df_sarica), 96000)
  expect_equal(ncol(df_sarica), 8)
  expect_equal(length(unique(df_sarica$subjectID)), 48)
})

test_that("load.afq.weston.havens loads the WH dataset", {
  df_hbn <- load.afq.weston.havens()

  expect_equal(nrow(df_hbn), 154000)
  expect_equal(ncol(df_hbn), 8)
  expect_equal(length(unique(df_hbn$subjectID)), 77)
})

test_that("load.afq.hbn loads the HBN dataset", {
  df_hbn <- load.afq.hbn(truncate = TRUE)

  expect_equal(nrow(df_hbn), 49)
  expect_equal(ncol(df_hbn), 7)
  expect_equal(length(unique(df_hbn$subjectID)), 1)
})
