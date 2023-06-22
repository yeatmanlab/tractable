#' Shuffle an AFQ dataframe
#'
#' This function shuffles participants' age, group, and sex,
#' thereby destroying correlations between the participants'
#' tract profiles and phenotypic data.
#'
#' @param input_df The input AFQ dataframe
#' @param dwi_metric The diffusion MRI metric (e.g. "FA", "MD")
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param participant.id The name of the column that encodes participant ID
#' @param shuffle_vars List of strings of column names that should be shuffled
#' @param sample_uniform Boolean flag. If TRUE, shuffling should sample
#'     uniformly from the unique values in the columns. If FALSE, shuffling
#'     will shuffle without replacement.
#'
#' @return A shuffled AFQ dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' df_shuffle <- shuffle_df(df_afq, "dti_fa")
#' }
shuffle_df <- function(input_df, dwi_metric, group.by = "group", participant.id = "subjectID", shuffle_vars = NULL, sample_uniform = FALSE) {
  # Spread the input dataframe to one row per participant
  col_names <- colnames(input_df)
  wide_df <- tidyr::pivot_wider(
    input_df,
    names_from = "nodeID",
    values_from = dwi_metric
  )

  if (is.null(shuffle_vars)) {
    shuffle_vars <- col_names[-which(col_names %in% c("nodeID", "tractID", dwi_metric, participant.id))]
  }

  # Then shuffle participants' shuffle_vars and the grouping variable
  for (svar in unique(c(shuffle_vars, group.by))) {
    if (sample_uniform) {
      # Sample uniformly from the unique values
      wide_df[[svar]] <- sample(unique(wide_df[[svar]]), length(wide_df[[svar]]), replace=TRUE)
    } else {
      # Shuffle the existing values
      wide_df[[svar]] <- sample(wide_df[[svar]], length(wide_df[[svar]]))
    }
  }

  # Gather back to long format (one row per node)
  output_df <- tidyr::pivot_longer(
    wide_df,
    -dplyr::any_of(c(participant.id, shuffle_vars, "tractID", group.by)),
    names_to = "nodeID",
    values_to = dwi_metric
  )

  output_df <- dplyr::select(output_df, dplyr::all_of(col_names))
  output_df$nodeID <- as.integer(output_df$nodeID)

  return(output_df)
}

#' Perform repeated sampling tests on an AFQ dataframe.
#'
#' When permutation_test == FALSE (the default), this function bootstrap
#' samples from an AFQ dataframe and returns pairwise differences at each node
#' for each bootstrap sample. These results can then be used to construct
#' bootstrap confidence intervals for the node-wise differences.
#'
#' When permutation_test == TRUE, this function simulates the null distribution
#' using permutation testing. That is, it shuffles to destroy any relationships
#' between covariates and dwi_metrics. It then computes node-wise differences
#' for each shuffled sample.
#'
#' @param df_tract AFQ Dataframe of node metric values for single tract
#' @param n_samples Number of sample tests to perform
#' @param dwi_metric The diffusion metric to model (e.g. "FA", "MD")
#' @param tract AFQ tract name
#' @param group.by The grouping variable used to group nodeID smoothing terms
#' @param participant.id The name of the column that encodes participant ID
#' @param sample_uniform Boolean flag. If TRUE, shuffling should sample
#'     uniformly from the unique values in the columns. If FALSE, shuffling
#'     will shuffle without replacement.
#' @param covariates List of strings of GAM covariates, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param k Dimension of the basis used to represent the node smoothing term
#' @param family Distribution to use for the gam. Must be 'gamma' or 'beta'
#' @param formula Optional explicit formula to use for the GAM. If provided,
#'     this will override the dynamically generated formula build from the
#'     target, covariate, and k inputs. Default = NULL.
#' @param factor_a First group factor, string
#' @param factor_b Second group factor, string
#' @param permute Boolean flag. If TRUE, perform a permutation test.
#'     Otherwise, do a bootstrap simulation.
#'
#' @return Dataframe with bootstrap or permutation test coefficients
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' df_tract <- df_afq[which(df_afq$tractID == tract), ]
#' bootstrap_coefs <- sampling_test(df_afq,
#'                                  dwi_metric = "dti_fa",
#'                                  covariates = list("group", "sex"),
#'                                  family = "gamma",
#'                                  k = 40,
#'                                  n_samples = 1000)
#' }
sampling_test <- function(df_tract,
                          n_samples,
                          dwi_metric,
                          tract,
                          group.by = "group",
                          participant.id = "subjectID",
                          sample_uniform = FALSE,
                          covariates = NULL,
                          smooth_terms = NULL,
                          k = NULL,
                          family = NULL,
                          formula = NULL,
                          factor_a = NULL,
                          factor_b = NULL,
                          permute = FALSE) {
  if (group.by %in% covariates) {
    coefs <- vector(mode = "list", length = n_samples)
    pvalues <- vector(mode = "list", length = n_samples)
  }
  node_diffs <- data.frame(nodeID = numeric(0),
                           est = numeric(0),
                           permIdx = numeric(0))

  pb <- progress::progress_bar$new(total = n_samples)
  for (idx in 1:n_samples) {
    pb$tick()

    if (permute) {
      df_shuffle <- shuffle_df(input_df = df_tract,
                               dwi_metric = dwi_metric,
                               group.by = group.by,
                               participant.id = participant.id,
                               shuffle_vars = covariates,
                               sample_uniform = sample_uniform)
    } else {
      df_shuffle <- bootstrap_df(input_df = df_tract,
                                 dwi_metric = dwi_metric,
                                 group.by = group.by,
                                 participant.id = participant.id)
    }

    gam_shuffle <- fit_gam(df_tract = df_shuffle,
                           target = dwi_metric,
                           covariates = covariates,
                           smooth_terms = smooth_terms,
                           group.by = group.by,
                           participant.id = participant.id,
                           formula = formula,
                           k = k,
                           family = family)
    ff <- summary(gam_shuffle)
    pvalues[[idx]] <- ff$p.table[,"Pr(>|t|)"][[paste0(group.by, factor_b)]]
    if (group.by %in% covariates) {
      coef_name <- grep(paste0("^", group.by),
                        names(gam_shuffle$coefficients),
                        value = TRUE)

      coefs[[idx]] <- gam_shuffle$coefficients[[coef_name]]
    }

    if (!is.null(factor_a) & !is.null(factor_b)) {
      df_pair <- spline_diff(gam_model = gam_shuffle,
                             tract = tract,
                             group.by = group.by,
                             factor_a,
                             factor_b,
                             save_output = FALSE,
                             out_dir = NULL)

      df_pair$permIdx <- idx
      df_pair <- dplyr::select(df_pair, c("nodeID", "est", "permIdx"))
      node_diffs <- rbind(node_diffs, df_pair)
    }
  }

  df_sampling_test <- tidyr::pivot_wider(
    node_diffs,
    names_from = "nodeID",
    names_prefix = "node_",
    values_from = "est",
    id_cols = "permIdx"
  )

  if (group.by %in% covariates) {
    df_sampling_test$group_coefs <- unlist(coefs)
  }

  df_sampling_test$pvalue <- unlist(pvalues)
  return(df_sampling_test)
}

#' Bootstrap data by family
#' Takes a long dataframe with Family_ID, and resamples families
#' with replacement. Returns resampled dataframe with column specifying 
#' how many times that family was resampled. 
#' @param df input dataframe, can be long or wide
#' @param resample_num how many resamples do you want
#' @param subject_id_col column with subject ids
#' @param grouping_id_col column to group by - if null, groups by subject

bootstrap <- function(tract_df, resample_num=NULL, subject_id_col="subject", grouping_id_col=NULL) { 
  
    if (is.null(grouping_id_col)) { 
        if (is.null(resample_num ))  { 
            resample_num <- length(unique(tract_df[[subject_id_col]])) 
        } 
     print("Grouping by Subject") 
     nested_df <- tract_df %>% 
            nest(data = everything(), .by=subject_id_col) %>%
            dplyr::slice_sample(n = resample_num, replace=TRUE) 
     counter <- sapply(unique(as.character(tract_df[[subject_id_col]])), function(x) 0)                        
    }  else { 
        if (is.null(resample_num)) {  
            resample_num <- length(unique(tract_df[[grouping_id_col]])) 
            } 
        print(paste("Grouping by", grouping_id_col, sep=' '))
        nested_df <- tract_df %>% 
            nest(data = everything(), .by=grouping_id_col) %>%
            dplyr::slice_sample(n = resample_num, replace=TRUE) 

            counter <- sapply(unique(as.character(tract_df[[grouping_id_col]])), 
                              function(x) 0) } 

    for (ii in 1:length(nested_df$data)) { 
        if (!is.null(grouping_id_col)) { 
            group_id <- as.character(unique(nested_df$data[[ii]][[grouping_id_col]]))
            count <- counter[[group_id]]
            nested_df$data[[ii]][[grouping_id_col]] <- paste(
                group_id, count, sep='_')
           nested_df$data[[ii]][[subject_id_col]] <-paste(
               nested_df$data[[ii]][[subject_id_col]], count, sep='_')               
            }
         else { 
            subject_id <- unique(nested_df$data[[ii]][[subject_id_col]]) 
            count <- counter[[subject_id]] 
            nested_df$data[[ii]][[subject_id_col]] <- paste(subject_id, count, sep='_') 
            counter[[subject_id]] =count + 1    
            }     
        } 
    boot_df <- unnest(nested_df, cols = c(data), names_repair="minimal")
    return(boot_df)     
}



cv_split <- function(df_tract, k=5, group.by = NULL) {

  df_tract <- df_tract %>% pivot_wider(names_from = c('tractID', 'nodeID'), 
                            values_from = all_of(sel_metrics), 
                            names_sep = '/')
    
    group_fold <- group_vfold_cv(profiles_wide, group=group.by, v=k) 

    return(group_fold)
    } 

