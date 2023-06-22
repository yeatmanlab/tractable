select_bundle <- function(df_afq = NULL,
                          tract,
                          dwi_metric,
                          participant_id,
                          covariates,
                          group_by,
                          ...) {


  if (is.null(df_afq)) {
    df_afq <- read.afq.files(
    index = participant_id,
    dwi_metrics = dwi_metric,
    ... = ...)
  }

  if (!is.null(group_by)) {
    df_afq[[group_by]] <- factor(df_afq[[group_by]])
  }

  df_afq[[participant_id]] <- factor(df_afq[[participant_id]])

  cols <- unique(c(participant_id,
                   "nodeID",
                   "tractID",
                   dwi_metric,
                   group_by,
                   covariates))
  df_afq <- dplyr::select(df_afq, dplyr::any_of(cols))

  if (tract != "all") {
    df_tract <- df_afq[which(df_afq$tractID == tract), ]
    df_tract <- stats::na.omit(df_tract)
    tract_names <- c(tract)
  } else {
    df_tract <- df_afq
    old_group_by <- group_by
    group_by <- paste0("tractID", group_by)
    df_tract[[group_by]] <- as.factor(paste(df_tract$tractID,
                                            df_tract[[old_group_by]]))
    df_tract <- stats::na.omit(df_tract)
    tract_names <- unique(df_tract$tractID)
  }
  return(list(df_tract = df_tract,
              tract_names = tract_names))
}