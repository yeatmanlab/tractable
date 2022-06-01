shuffle_df <- function(input_df, dwi_metric) {
  # Spread the input dataframe to one row per participant
  col_names <- colnames(input_df)
  wide_df <- tidyr::pivot_wider(
    input_df,
    names_from = "nodeID",
    values_from = dwi_metric
  )

  # Then shuffle participants' age and group
  wide_df$group <- sample(wide_df$group, length(wide_df$group))
  wide_df$age <- sample(wide_df$age, length(wide_df$age))

  # Randomly assign the subjects' sex
  wide_df$sex <- sample(c(0, 1), length(wide_df$sex), replace=TRUE)

  # Gather back to long format (one row per node)
  output_df <- tidyr::pivot_longer(
    wide_df,
    -dplyr::any_of(c("subjectID", "age", "sex", "tractID", "group")),
    names_to = "nodeID",
    values_to = dwi_metric
  )

  output_df <- dplyr::select(output_df, dplyr::all_of(col_names))
  output_df$nodeID <- as.integer(output_df$nodeID)

  return(output_df)
}
