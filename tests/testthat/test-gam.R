

test_that("gams run as expected", {
 sarica <- read.afq.sarica()
 sarica$group <- factor(sarica$class)
 sarica$subjectID <- unclass(factor(sarica$subjectID))
 tractr_single_bundle(df_afq = sarica,
                      out_dir = ".",
                      tract = "Right Corticospinal",
                      participant_id = "subjectID",
                      group_by = "group",
                      covariates = c("age","group"),
                      dwi_metric = "fa",
                      comp_list = c("ALS", "CTRL"))})