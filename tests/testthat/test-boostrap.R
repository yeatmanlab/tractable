test_that("bootstrap_tracts generates a df with correct dimensions", 
    {
    #checks if size is correct
   df_afq <- read.afq.sarica()
   set.seed(0)
   df_boot <- tractr::bootstrap(df_afq,
              subject_id_col='subjectID')
        
  expect_equal(dim(df_afq)[1], dim(df_boot)[1])
})

# test_that("bootstrap_df resamples an afq dataframe by group", 
#     {
#    df_afq <- read.afq.sarica()
#    set.seed(0)
#    df_boot <- tractr::bootstrap(df_afq,
#               subject_id_col='subjectID',
#                 grouping_id_col='age')

#     length(df_afq)
        
#   expect_equal(dim(df_afq)[1], dim(df_boot)[1])

# })