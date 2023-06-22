select_bundle <- function(df_afq = NULL,
                          tract,
                          dwi_metric,
                          participant.id,
                          covariates,
                          group.by
                          ){

  if (is.null(df_afq)) {
    df_afq <- read.afq.files(...,
    index = participant.id,
    dwi_metrics = c(dwi_metric))
  }

  if (!is.null(group.by)) {
    df_afq[[group.by]] <- factor(df_afq[[group.by]])
  }

  df_afq[[participant.id]] <- factor(df_afq[[participant.id]])

  cols <- unique(c(participant.id,
                   "nodeID",
                   "tractID",
                   dwi_metric,
                   group.by,
                   covariates))
  df_afq <- dplyr::select(df_afq, dplyr::any_of(cols))

  if (tract != "all") {
    df_tract <- df_afq[which(df_afq$tractID == tract), ]
    df_tract <- stats::na.omit(df_tract)
    tract_names <- c(tract)
  } else {
    df_tract <- df_afq
    old_group.by <- group.by
    group.by <- paste0("tractID", group.by)
    df_tract[[group.by]] <- as.factor(paste(df_tract$tractID,
                                            df_tract[[old_group.by]]))
    df_tract <- stats::na.omit(df_tract)
    tract_names <- unique(df_tract$tractID)
  }
  return(list(df_tract=df_tract, tract_names=tract_names))
}