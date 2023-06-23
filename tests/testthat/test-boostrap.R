test_that("bootstrap_df resamples an afq dataframe by subject", 
    {

    #checks if size is correct
   df_afq <- read.afq.sarica()
   set.seed(0)
   df_boot <- bootstrap_resample(df_afq,
              subject_id_col='subjectID')
   df_ref <- read.csv("sarica_bootstrap_test_df")

  expect_true(dplyr::all_equal(dim(df_boot)[1], dim(df_ref)[1]))
})

test_that("bootstrap_df resamples an afq dataframe by group", 
    {

    #checks if size is correct
   df_afq <- read.afq.sarica()
   set.seed(0)
   df_boot <- bootstrap_resample(df_afq,
              subject_id_col='subjectID')
   df_ref <- read.csv("sarica_bootstrap_test_df")

  expect_true(dplyr::all_equal(dim(df_boot)[1], dim(df_ref)[1]))
})