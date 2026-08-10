options(stringsAsFactors = FALSE, width = 180)
d <- readRDS(file.path("02_processed_data", "pacc_pdac_analysis_cohort.rds"))
table_dir <- file.path("04_results", "tables")

fmt_n_pct <- function(n, denom) sprintf("%s (%.1f%%)", format(n, big.mark = ",", scientific = FALSE), 100 * n / denom)

rows <- list(); j <- 0
add_row <- function(variable, level, pdac, pacc, smd = "", p = "") {
  j <<- j + 1
  rows[[j]] <<- data.frame(Variable = variable, Level = level, PDAC = pdac, pACC = pacc,
                            SMD = smd, P_value = p, check.names = FALSE)
}

z0 <- d[d$histology == "PDAC", ]; z1 <- d[d$histology == "pACC", ]
n0 <- nrow(z0); n1 <- nrow(z1)
add_row("Patients", "N", format(n0, big.mark = ","), format(n1, big.mark = ","))

age_smd <- (mean(z1$age) - mean(z0$age)) / sqrt((var(z1$age) + var(z0$age)) / 2)
age_p <- wilcox.test(z1$age, z0$age)$p.value
add_row("Age at diagnosis, years", "Median (IQR)",
        sprintf("%.0f (%.0f–%.0f)", median(z0$age), quantile(z0$age, .25), quantile(z0$age, .75)),
        sprintf("%.0f (%.0f–%.0f)", median(z1$age), quantile(z1$age, .25), quantile(z1$age, .75)),
        sprintf("%.3f", age_smd), format.pval(age_p, digits = 3, eps = 0.001))

categorical_block <- function(variable, x, levels_to_show = NULL) {
  x <- droplevels(factor(x))
  if (!is.null(levels_to_show)) x <- factor(x, levels = levels_to_show)
  tab <- table(d$histology, x, useNA = "no")
  p <- suppressWarnings(chisq.test(tab)$p.value)
  props <- prop.table(tab, 1)
  denom <- sqrt((props[1, ] * (1 - props[1, ]) + props[2, ] * (1 - props[2, ])) / 2)
  level_smd <- ifelse(denom > 0, (props[2, ] - props[1, ]) / denom, 0)
  block_smd <- max(abs(level_smd), na.rm = TRUE)
  first <- TRUE
  for (lv in colnames(tab)) {
    add_row(variable, lv, fmt_n_pct(tab[1, lv], sum(tab[1, ])), fmt_n_pct(tab[2, lv], sum(tab[2, ])),
            if (first) sprintf("%.3f", block_smd) else "",
            if (first) format.pval(p, digits = 3, eps = 0.001) else "")
    first <- FALSE
  }
}

categorical_block("Sex", d$sex)
categorical_block("Race", d$race3, c("White", "Black", "Other/unknown"))
categorical_block("Ethnicity", d$origin, c("Non-Hispanic", "Hispanic"))
categorical_block("Primary site", d$site, c("Head", "Body/tail", "Other/unspecified"))
categorical_block("Summary stage", d$stage, c("Localized", "Regional", "Distant", "Unknown"))
categorical_block("Tumor grade", d$grade, c("1", "2", "3", "4", "Unknown"))
categorical_block("Cancer-directed surgery", d$surgery, c("No/unknown", "Yes"))
categorical_block("Chemotherapy", d$chemotherapy, c("No/unknown", "Yes"))
categorical_block("Radiotherapy", d$radiation, c("No/unknown", "Yes"))
categorical_block("Diagnosis era", d$era, c("2004-2009", "2010-2015", "2016-2017", "2018-2023"))

table1 <- do.call(rbind, rows)
write.csv(table1, file.path(table_dir, "Table1_patient_characteristics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Compact survival table for the manuscript body.
km <- read.csv(file.path(table_dir, "04_km_survival_estimates.csv"), check.names = FALSE)
table2 <- km[km$endpoint == "OS" & km$stage %in% c("All", "Localized", "Regional", "Distant") & km$time %in% c(1, 3, 5), ]
table2$survival_95CI <- sprintf("%.1f%% (%.1f–%.1f)", 100 * table2$survival, 100 * table2$lower, 100 * table2$upper)
table2 <- table2[, c("histology", "stage", "time", "n", "events", "survival_95CI", "median_survival")]
names(table2) <- c("Histology", "Stage", "Year", "N", "Deaths", "OS_95CI", "Median_OS_years")
write.csv(table2, file.path(table_dir, "Table2_stage_specific_OS.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Conditional survival table with the prespecified 3-year horizon.
cs <- read.csv(file.path(table_dir, "05_conditional_survival.csv"), check.names = FALSE)
table3 <- cs[cs$endpoint == "OS" & cs$horizon == 3 & cs$stage %in% c("All", "Localized", "Regional", "Distant"), ]
table3$conditional_OS_95CI <- sprintf("%.1f%% (%.1f–%.1f)", 100 * table3$estimate, 100 * table3$lower, 100 * table3$upper)
table3 <- table3[, c("histology", "stage", "landmark", "n_risk", "conditional_OS_95CI")]
names(table3) <- c("Histology", "Stage", "Landmark_year", "N_at_risk", "Three_year_conditional_OS_95CI")
write.csv(table3, file.path(table_dir, "Table3_three_year_conditional_OS.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Overall five-year conditional OS, retained as a compact supplementary table
# because stage-specific late risk sets are sparse for pACC.
table_s4 <- cs[cs$endpoint == "OS" & cs$horizon == 5 & cs$stage == "All", ]
table_s4$conditional_OS_95CI <- sprintf("%.1f%% (%.1f–%.1f)", 100 * table_s4$estimate, 100 * table_s4$lower, 100 * table_s4$upper)
table_s4 <- table_s4[, c("histology", "landmark", "n_risk", "events", "conditional_OS_95CI")]
names(table_s4) <- c("Histology", "Landmark_years", "N_at_landmark", "Subsequent_deaths", "Five_year_conditional_OS_95CI")
write.csv(table_s4, file.path(table_dir, "Supplementary_Table_S4_five_year_conditional_OS.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Three-year conditional probability of remaining free from cancer-specific
# death, estimated as 1 minus the Aalen-Johansen cancer-death CIF.
ccr <- read.csv(file.path(table_dir, "20_conditional_competing_risk_survival.csv"), check.names = FALSE)
table_s5 <- ccr[ccr$horizon == 3 & ccr$stage %in% c("All", "Localized", "Regional", "Distant"), ]
table_s5$conditional_cancer_death_free_95CI <- sprintf(
  "%.1f%% (%.1f–%.1f)", 100 * table_s5$estimate, 100 * table_s5$lower95, 100 * table_s5$upper95
)
table_s5 <- table_s5[, c(
  "histology", "stage", "landmark", "n_risk", "cancer_deaths", "competing_deaths",
  "conditional_cancer_death_free_95CI"
)]
names(table_s5) <- c(
  "Histology", "Stage", "Landmark_years", "N_at_landmark", "Cancer_deaths_next_3y",
  "Competing_deaths_next_3y", "Three_year_conditional_cancer_death_free_95CI"
)
write.csv(table_s5, file.path(table_dir, "Supplementary_Table_S5_conditional_competing_risk.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Supplementary competing-risk table.
fg <- read.csv(file.path(table_dir, "15_fine_gray_models.csv"), check.names = FALSE)
ci <- read.csv(file.path(table_dir, "16_cancer_death_cif.csv"), check.names = FALSE)
ci <- ci[ci$time == 5 & ci$stage %in% fg$stage, c("stage", "histology", "cancer_death_cif")]
ci_wide <- reshape(ci, idvar = "stage", timevar = "histology", direction = "wide")
fg <- merge(fg, ci_wide, by = "stage", all.x = TRUE, sort = FALSE)
fg <- fg[match(c("All", "Localized", "Regional", "Distant"), fg$stage), ]
fg$subdistribution_HR_95CI <- sprintf("%.3f (%.3f-%.3f)", fg$sHR, fg$lower95, fg$upper95)
fg$P_value <- format.pval(fg$p, digits = 3, eps = 0.001)
fg$PDAC_five_year_CIF <- sprintf("%.1f%%", 100 * fg$cancer_death_cif.PDAC)
fg$pACC_five_year_CIF <- sprintf("%.1f%%", 100 * fg$cancer_death_cif.pACC)
table_s2 <- fg[, c("stage", "PDAC_five_year_CIF", "pACC_five_year_CIF", "subdistribution_HR_95CI", "P_value")]
names(table_s2) <- c("Stage", "PDAC_5y_cancer_death_CIF", "pACC_5y_cancer_death_CIF", "pACC_vs_PDAC_sHR_95CI", "P_value")
write.csv(table_s2, file.path(table_dir, "Supplementary_Table_S2_competing_risk.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Supplementary pACC temporal-trend table.
trend <- read.csv(file.path(table_dir, "17_pacc_survival_trend.csv"), check.names = FALSE)
trend$HR_95CI <- sprintf("%.3f (%.3f-%.3f)", trend$HR, trend$lower95, trend$upper95)
trend$P_value <- ifelse(trend$p < 0.001, "<0.001", sprintf("%.3f", trend$p))
trend$term <- c("Per one-year increase in diagnosis year", "2010-2015 vs 2004-2009", "2016-2023 vs 2004-2009")
table_s3 <- trend[, c("analysis", "term", "n", "events", "HR_95CI", "P_value")]
names(table_s3) <- c("Analysis", "Contrast", "N", "Deaths", "Adjusted_HR_95CI", "P_value")
write.csv(table_s3, file.path(table_dir, "Supplementary_Table_S3_temporal_trend.csv"), row.names = FALSE, fileEncoding = "UTF-8")

message("Manuscript tables written.")
