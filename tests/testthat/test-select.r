test_that("select_bundle runs as expected", {
 sarica <- read.afq.sarica()
 sarica$group <- factor(sarica$class)
 sarica$subjectID <- unclass(factor(sarica$subjectID))

 gam_fit <- expect_no_error(
    select_bundle(df_afq = sarica,
                  tract = "Right Corticospinal",
                  dwi_metric = "fa",
                  covariates = c("age", "group"),
                  participant_id = "subjectID",
                  group_by = "group")
                )
    })
