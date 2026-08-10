options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(survival))

processed_file <- file.path("02_processed_data", "pacc_pdac_analysis_cohort.rds")
table_dir <- file.path("04_results", "tables")
d <- readRDS(processed_file)

tidy_histology <- function(fit, analysis, endpoint, stage = "All") {
  sm <- summary(fit)
  i <- match("histologypACC", rownames(sm$coefficients))
  if (is.na(i)) return(NULL)
  data.frame(
    analysis = analysis, endpoint = endpoint, stage = stage,
    HR = sm$conf.int[i, "exp(coef)"],
    lower95 = sm$conf.int[i, "lower .95"],
    upper95 = sm$conf.int[i, "upper .95"],
    p = sm$coefficients[i, "Pr(>|z|)"], n = fit$n, events = fit$nevent
  )
}

fit_model <- function(z, endpoint = c("OS", "CSS"), covariates, analysis, stage = "All") {
  endpoint <- match.arg(endpoint)
  if (endpoint == "CSS") z <- z[z$css_known, , drop = FALSE]
  z$event_tmp <- if (endpoint == "OS") z$os_event else z$css_event
  f <- as.formula(paste("Surv(os_time, event_tmp) ~ histology +", covariates, "+ cluster(patient_id)"))
  fit <- coxph(f, data = z, ties = "efron", na.action = na.omit, x = TRUE)
  tidy_histology(fit, analysis, endpoint, stage)
}

base_cov <- "age + sex + race3 + origin + stage + site + era"
treatment_cov <- paste(base_cov, "+ grade + surgery + chemotherapy + radiation")
sens <- list(); j <- 0

for (ep in c("OS", "CSS")) {
  j <- j + 1; sens[[j]] <- fit_model(d, ep, base_cov, "Primary baseline-adjusted")
  j <- j + 1; sens[[j]] <- fit_model(d, ep, treatment_cov, "Additionally adjusted for grade/treatment")
  j <- j + 1; sens[[j]] <- fit_model(d[d$single_primary, , drop = FALSE], ep, base_cov, "Single-primary tumors")
  j <- j + 1; sens[[j]] <- fit_model(d[d$stage != "Unknown", , drop = FALSE], ep, base_cov, "Known-stage cohort")
  j <- j + 1; sens[[j]] <- fit_model(d[d$year <= 2018, , drop = FALSE], ep, base_cov, "Diagnosis through 2018")
}
sens <- do.call(rbind, sens)
write.csv(sens, file.path(table_dir, "09_sensitivity_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Proportional-hazards diagnostics for the primary model. The cluster term is
# omitted only for the residual test; point estimates are reported from robust models elsewhere.
ph_rows <- list(); j <- 0
for (ep in c("OS", "CSS")) {
  z <- if (ep == "CSS") d[d$css_known, , drop = FALSE] else d
  z$event_tmp <- if (ep == "OS") z$os_event else z$css_event
  fit <- coxph(as.formula(paste("Surv(os_time, event_tmp) ~ histology +", base_cov)), data = z, ties = "efron", x = TRUE)
  ph <- cox.zph(fit, transform = "km")$table
  for (term in c("histology", "GLOBAL")) {
    if (term %in% rownames(ph)) {
      j <- j + 1
      ph_rows[[j]] <- data.frame(endpoint = ep, term = term, chisq = ph[term, "chisq"], df = ph[term, "df"], p = ph[term, "p"])
    }
  }
}
ph_out <- do.call(rbind, ph_rows)
write.csv(ph_out, file.path(table_dir, "10_proportional_hazards_tests.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Interval-specific adjusted HRs describe how the histology association changes
# among patients alive at successive time intervals.
intervals <- data.frame(start = c(0, 1, 3, 5), stop = c(1, 3, 5, 10), label = c("0-1 y", "1-3 y", "3-5 y", "5-10 y"))
tv_rows <- list(); j <- 0
for (ep in c("OS", "CSS")) for (st in c("All", "Localized", "Regional", "Distant")) for (ii in seq_len(nrow(intervals))) {
  start <- intervals$start[ii]; stop <- intervals$stop[ii]
  z <- d[(st == "All" | d$stage == st) & d$os_time > start, , drop = FALSE]
  if (ep == "CSS") z <- z[z$css_known, , drop = FALSE]
  z$interval_time <- pmin(z$os_time, stop) - start
  z$interval_event <- if (ep == "OS") as.integer(z$os_event == 1 & z$os_time <= stop) else as.integer(z$css_event == 1 & z$os_time <= stop)
  rhs <- if (st == "All") base_cov else "age + sex + race3 + origin + site + era"
  fit <- coxph(as.formula(paste("Surv(interval_time, interval_event) ~ histology +", rhs, "+ cluster(patient_id)")),
               data = z, ties = "efron", na.action = na.omit)
  out <- tidy_histology(fit, "Interval-specific adjusted", ep, st)
  out$interval <- intervals$label[ii]
  out$n_at_interval_start <- nrow(z)
  j <- j + 1; tv_rows[[j]] <- out
}
tv <- do.call(rbind, tv_rows)
tv <- tv[, c("analysis", "endpoint", "stage", "interval", "n_at_interval_start", "HR", "lower95", "upper95", "p", "n", "events")]
write.csv(tv, file.path(table_dir, "11_interval_specific_hr.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Wald test for histology-by-stage interaction.
interaction_rows <- list()
for (ep in c("OS", "CSS")) {
  z <- d[d$stage != "Unknown" & (ep == "OS" | d$css_known), , drop = FALSE]
  z <- droplevels(z)
  z$event_tmp <- if (ep == "OS") z$os_event else z$css_event
  fit <- coxph(Surv(os_time, event_tmp) ~ histology * stage + age + sex + race3 + origin + site + era + cluster(patient_id),
               data = z, ties = "efron", x = TRUE)
  b <- coef(fit); V <- vcov(fit)
  ids <- grep("^histologypACC:stage", names(b))
  wald <- as.numeric(t(b[ids]) %*% qr.solve(V[ids, ids, drop = FALSE], b[ids]))
  interaction_rows[[ep]] <- data.frame(endpoint = ep, chisq = wald, df = length(ids), p = pchisq(wald, length(ids), lower.tail = FALSE))
}
interaction_out <- do.call(rbind, interaction_rows)
write.csv(interaction_out, file.path(table_dir, "12_histology_stage_interaction.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Conditional-survival contrast (percentage-point difference) from the primary estimates.
cs <- read.csv(file.path(table_dir, "05_conditional_survival.csv"), check.names = FALSE)
pd <- cs[cs$histology == "PDAC", c("endpoint", "stage", "landmark", "horizon", "estimate", "n_risk")]
pa <- cs[cs$histology == "pACC", c("endpoint", "stage", "landmark", "horizon", "estimate", "n_risk")]
names(pd)[5:6] <- c("PDAC_estimate", "PDAC_n_risk")
names(pa)[5:6] <- c("pACC_estimate", "pACC_n_risk")
contrast <- merge(pa, pd, by = c("endpoint", "stage", "landmark", "horizon"), all = TRUE)
contrast$absolute_difference_pctpt <- 100 * (contrast$pACC_estimate - contrast$PDAC_estimate)
write.csv(contrast, file.path(table_dir, "13_conditional_survival_contrast.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Log-rank comparisons by stage.
lr_rows <- list(); j <- 0
for (ep in c("OS", "CSS")) for (st in c("All", "Localized", "Regional", "Distant")) {
  z <- d[st == "All" | d$stage == st, , drop = FALSE]
  if (ep == "CSS") z <- z[z$css_known, , drop = FALSE]
  event <- if (ep == "OS") z$os_event else z$css_event
  lr <- survdiff(Surv(z$os_time, event) ~ z$histology, rho = 0)
  j <- j + 1
  lr_rows[[j]] <- data.frame(endpoint = ep, stage = st, chisq = unname(lr$chisq), df = length(lr$n) - 1,
                             p = pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE))
}
lr_out <- do.call(rbind, lr_rows)
write.csv(lr_out, file.path(table_dir, "14_logrank_tests.csv"), row.names = FALSE, fileEncoding = "UTF-8")

message("Extended analyses complete.")
