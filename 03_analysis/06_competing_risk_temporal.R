options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(survival))

processed_file <- file.path("02_processed_data", "pacc_pdac_analysis_cohort.rds")
table_dir <- file.path("04_results", "tables")
d <- readRDS(processed_file)

# Competing-risk status: 0=censored/alive, 1=index-cancer death,
# 2=death from another known cause. Deaths with missing/unknown cause are
# excluded, consistent with the conventional CSS analysis.
z <- d[d$css_known, , drop = FALSE]
z$cr_status <- factor(
  ifelse(z$os_event == 0, "censor", ifelse(z$css_event == 1, "cancer", "other")),
  levels = c("censor", "cancer", "other")
)

tidy_fg <- function(data, stage_label = "All") {
  rhs <- if (stage_label == "All") {
    "histology + age + sex + race3 + origin + stage + site + era + patient_id"
  } else {
    "histology + age + sex + race3 + origin + site + era + patient_id"
  }
  fg <- finegray(
    as.formula(paste("Surv(os_time, cr_status) ~", rhs)),
    data = data,
    etype = "cancer"
  )
  cox_rhs <- if (stage_label == "All") {
    "histology + age + sex + race3 + origin + stage + site + era + cluster(patient_id)"
  } else {
    "histology + age + sex + race3 + origin + site + era + cluster(patient_id)"
  }
  fit <- coxph(
    as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", cox_rhs)),
    data = fg,
    weights = fgwt,
    ties = "efron"
  )
  sm <- summary(fit)
  i <- match("histologypACC", rownames(sm$coefficients))
  data.frame(
    analysis = "Fine-Gray competing-risk model",
    endpoint = "Cancer-specific death",
    stage = stage_label,
    sHR = sm$conf.int[i, "exp(coef)"],
    lower95 = sm$conf.int[i, "lower .95"],
    upper95 = sm$conf.int[i, "upper .95"],
    p = sm$coefficients[i, "Pr(>|z|)"],
    n = nrow(data),
    cancer_deaths = sum(data$cr_status == "cancer"),
    competing_deaths = sum(data$cr_status == "other")
  )
}

fg_rows <- list(tidy_fg(z, "All"))
for (st in c("Localized", "Regional", "Distant")) {
  fg_rows[[length(fg_rows) + 1]] <- tidy_fg(droplevels(z[z$stage == st, , drop = FALSE]), st)
}
fg_out <- do.call(rbind, fg_rows)
write.csv(fg_out, file.path(table_dir, "15_fine_gray_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Aalen-Johansen cumulative incidence estimates for cancer-specific death.
cif_rows <- list(); k <- 0
for (st in c("All", "Localized", "Regional", "Distant")) {
  zz <- if (st == "All") z else droplevels(z[z$stage == st, , drop = FALSE])
  fit <- survfit(Surv(os_time, cr_status) ~ histology, data = zz)
  ss <- summary(fit, times = c(1, 3, 5, 10), extend = TRUE)
  for (i in seq_along(ss$time)) {
    k <- k + 1
    cif_rows[[k]] <- data.frame(
      stage = st,
      histology = sub("^histology=", "", ss$strata[i]),
      time = ss$time[i],
      n_risk = ss$n.risk[i],
      cancer_death_cif = ss$pstate[i, "cancer"],
      lower95 = ss$lower[i, "cancer"],
      upper95 = ss$upper[i, "cancer"]
    )
  }
}
cif_out <- do.call(rbind, cif_rows)
write.csv(cif_out, file.path(table_dir, "16_cancer_death_cif.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Conditional competing-risk estimates. At each landmark, the analysis is
# restricted to patients still alive, time is reset to zero, and the
# Aalen-Johansen estimator is used for cancer-specific death with noncancer
# death retained as a competing event. The reported probability is 1-CIF.
conditional_cif <- function(data, landmark, horizon) {
  zz <- data[data$os_time >= landmark, , drop = FALSE]
  if (nrow(zz) < 2) {
    return(data.frame(n_risk = nrow(zz), cancer_deaths = NA_integer_, competing_deaths = NA_integer_,
                      estimate = NA_real_, lower95 = NA_real_, upper95 = NA_real_))
  }
  zz$residual_time <- zz$os_time - landmark
  fit <- survfit(Surv(residual_time, cr_status) ~ 1, data = zz)
  ss <- summary(fit, times = horizon, extend = TRUE)
  cif <- unname(ss$pstate[1, "cancer"])
  cif_lower <- unname(ss$lower[1, "cancer"])
  cif_upper <- unname(ss$upper[1, "cancer"])
  data.frame(
    n_risk = nrow(zz),
    cancer_deaths = sum(zz$cr_status == "cancer" & zz$residual_time <= horizon),
    competing_deaths = sum(zz$cr_status == "other" & zz$residual_time <= horizon),
    estimate = 1 - cif,
    lower95 = 1 - cif_upper,
    upper95 = 1 - cif_lower
  )
}

ccs_rows <- list(); k <- 0
for (st in c("All", "Localized", "Regional", "Distant")) {
  for (h in levels(z$histology)) {
    zz <- z[z$histology == h & (st == "All" | z$stage == st), , drop = FALSE]
    for (lm in c(0, 1, 2, 3, 5)) for (hz in c(3, 5)) {
      out <- conditional_cif(zz, lm, hz)
      out$histology <- h
      out$stage <- st
      out$landmark <- lm
      out$horizon <- hz
      k <- k + 1
      ccs_rows[[k]] <- out
    }
  }
}
conditional_cif_out <- do.call(rbind, ccs_rows)
conditional_cif_out <- conditional_cif_out[, c(
  "histology", "stage", "landmark", "horizon", "n_risk", "cancer_deaths",
  "competing_deaths", "estimate", "lower95", "upper95"
)]
write.csv(
  conditional_cif_out,
  file.path(table_dir, "20_conditional_competing_risk_survival.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Contemporary pACC survival trend. Diagnosis year is evaluated both as a
# continuous covariate and in clinically interpretable eras. The model is
# descriptive and adjusts for measured case mix; it does not estimate incidence.
p <- droplevels(d[d$histology == "pACC" & d$stage != "Unknown", , drop = FALSE])
p$year_centered <- p$year - 2004
p$trend_era <- factor(
  ifelse(p$year <= 2009, "2004-2009", ifelse(p$year <= 2015, "2010-2015", "2016-2023")),
  levels = c("2004-2009", "2010-2015", "2016-2023")
)
fit_year <- coxph(
  Surv(os_time, os_event) ~ year_centered + age + sex + race3 + origin + stage + site,
  data = p,
  ties = "efron"
)
fit_era <- coxph(
  Surv(os_time, os_event) ~ trend_era + age + sex + race3 + origin + stage + site,
  data = p,
  ties = "efron"
)

tidy_term <- function(fit, terms, analysis) {
  sm <- summary(fit)
  do.call(rbind, lapply(terms, function(term) {
    i <- match(term, rownames(sm$coefficients))
    data.frame(
      analysis = analysis,
      term = term,
      HR = sm$conf.int[i, "exp(coef)"],
      lower95 = sm$conf.int[i, "lower .95"],
      upper95 = sm$conf.int[i, "upper .95"],
      p = sm$coefficients[i, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent
    )
  }))
}

trend_out <- rbind(
  tidy_term(fit_year, "year_centered", "Adjusted linear trend per calendar year"),
  tidy_term(fit_era, c("trend_era2010-2015", "trend_era2016-2023"), "Adjusted diagnosis-era comparison")
)
write.csv(trend_out, file.path(table_dir, "17_pacc_survival_trend.csv"), row.names = FALSE, fileEncoding = "UTF-8")

era_desc <- do.call(rbind, lapply(levels(p$trend_era), function(e) {
  q <- p[p$trend_era == e, , drop = FALSE]
  sf <- survfit(Surv(os_time, os_event) ~ 1, data = q, conf.type = "log-log")
  s <- summary(sf, times = c(1, 3, 5), extend = TRUE)
  data.frame(
    era = e,
    n = nrow(q),
    deaths = sum(q$os_event),
    time = s$time,
    os = s$surv,
    lower95 = s$lower,
    upper95 = s$upper
  )
}))
write.csv(era_desc, file.path(table_dir, "18_pacc_os_by_era.csv"), row.names = FALSE, fileEncoding = "UTF-8")

message("Competing-risk and temporal analyses complete.")
