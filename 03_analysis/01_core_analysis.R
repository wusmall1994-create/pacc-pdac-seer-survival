options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(survival))

input_file <- file.path("01_data_extraction", "pacc_pdac_seer17_2004_2023_raw.txt")
processed_file <- file.path("02_processed_data", "pacc_pdac_analysis_cohort.rds")
table_dir <- file.path("04_results", "tables")
dir.create(dirname(processed_file), recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading raw SEER export...")
raw <- read.delim(
  input_file, sep = "\t", quote = "\"", na.strings = "NA",
  check.names = FALSE, colClasses = "character"
)

v <- function(name) raw[[name]]
flow <- data.frame(step = character(), n = integer())
add_flow <- function(step, keep) {
  flow <<- rbind(flow, data.frame(step = step, n = sum(keep, na.rm = TRUE)))
  invisible(keep)
}

keep <- rep(TRUE, nrow(raw))
add_flow("Raw pancreas records exported", keep)
keep <- keep & v("Histologic Type ICD-O-3") %in% c("8140", "8500", "8550")
add_flow("Target main histologies (8140, 8500, 8550)", keep)
target_histology <- ifelse(v("Histologic Type ICD-O-3") == "8550", "pACC", "PDAC")
confirmation_ok <- v("Diagnostic Confirmation") %in% c(
  "Positive histology",
  "Pos hist AND immunophenotyping AND/OR pos genetic studies"
)
confirmation_flow <- do.call(rbind, lapply(c("PDAC", "pACC"), function(h) {
  eligible <- keep & target_histology == h
  data.frame(
    histology = h,
    before_positive_histology_filter = sum(eligible),
    after_positive_histology_filter = sum(eligible & confirmation_ok),
    excluded_by_positive_histology_filter = sum(eligible & !confirmation_ok)
  )
}))
keep <- keep & v("Diagnostic Confirmation") %in% c(
  "Positive histology",
  "Pos hist AND immunophenotyping AND/OR pos genetic studies"
)
add_flow("Positive histologic confirmation", keep)
keep <- keep & !v("Type of Reporting Source") %in% c("Autopsy only", "Death certificate only")
add_flow("Exclude autopsy/DCO reporting source", keep)
keep <- keep & !is.na(v("Survival Days")) & !is.na(v("Vital status recode (study cutoff used)"))
add_flow("Known survival time and vital status", keep)

d <- raw[keep, , drop = FALSE]
rm(raw)

num <- function(x) suppressWarnings(as.numeric(x))
d$patient_id <- d[["Patient ID"]]
d$year <- num(d[["Year of diagnosis"]])
d$age <- num(sub("[^0-9].*$", "", d[["Age recode with single ages and 85+"]]))
d$age_topcoded <- grepl("85\\+", d[["Age recode with single ages and 85+"]])
d$sex <- factor(d[["Sex"]])
d$race <- factor(d[["Race recode (W, B, AI, API)"]])
d$race3 <- factor(ifelse(d[["Race recode (W, B, AI, API)"]] == "White", "White",
                  ifelse(d[["Race recode (W, B, AI, API)"]] == "Black", "Black", "Other/unknown")),
                  levels = c("White", "Black", "Other/unknown"))
d$origin <- factor(ifelse(grepl("Non-Spanish", d[["Origin recode NHIA (Hispanic, Non-Hisp)"]]), "Non-Hispanic",
                          ifelse(grepl("Spanish|Hispanic", d[["Origin recode NHIA (Hispanic, Non-Hisp)"]]), "Hispanic", "Unknown")))
d$hist_code <- d[["Histologic Type ICD-O-3"]]
d$histology <- factor(ifelse(d$hist_code == "8550", "pACC", "PDAC"), levels = c("PDAC", "pACC"))
d$pdac_subtype <- factor(ifelse(d$hist_code == "8140", "8140 adenocarcinoma NOS",
                         ifelse(d$hist_code == "8500", "8500 duct carcinoma", "8550 pACC")))

# This is the SEER system-supplied cross-era Summary Stage variable, not an
# author-derived merge of edition-specific AJCC stage groups. Expanded regional
# categories are collapsed below to the conventional localized/regional/distant
# framework used throughout the analysis.
stage_raw <- d[["Combined Summary Stage with Expanded Regional Codes (2004+)"]]
d$stage <- factor(ifelse(grepl("^Localized", stage_raw), "Localized",
                  ifelse(grepl("^Regional", stage_raw), "Regional",
                  ifelse(grepl("^Distant", stage_raw), "Distant", "Unknown"))),
                  levels = c("Localized", "Regional", "Distant", "Unknown"))

site_raw <- d[["Primary Site - labeled"]]
d$site <- factor(ifelse(grepl("^C25.0", site_raw), "Head",
                 ifelse(grepl("^C25.[12]", site_raw), "Body/tail", "Other/unspecified")))
# Era cut points align with major registry staging/data-collection transitions:
# AJCC 6 (2004-2009), AJCC 7 (2010-2015), the 2016-2017 transition interval,
# and EOD 2018/AJCC 8 (2018 onward). Era is an adjustment variable; the primary
# stage variable remains the system-supplied Combined Summary Stage above.
d$era <- factor(ifelse(d$year <= 2009, "2004-2009",
                ifelse(d$year <= 2015, "2010-2015",
                ifelse(d$year <= 2017, "2016-2017", "2018-2023"))))

grade_pre <- d[["Grade Recode (thru 2017)"]]
grade_post <- d[["Derived Summary Grade 2018 (2018+)"]]
grade_from <- function(x) {
  ifelse(grepl("Grade I$|category \\(1\\)", x), "1",
  ifelse(grepl("Grade II$|category \\(2\\)", x), "2",
  ifelse(grepl("Grade III$|category \\(3\\)", x), "3",
  ifelse(grepl("Grade IV$|category \\(4\\)", x), "4", "Unknown"))))
}
d$grade <- factor(ifelse(d$year <= 2017, grade_from(grade_pre), grade_from(grade_post)),
                  levels = c("1", "2", "3", "4", "Unknown"))

surg_old <- d[["RX Summ--Surg Prim Site (1998-2022)"]]
surg_new <- d[["RX Summ--Surg Prim Site 2023 (2023+)"]]
surg_code <- ifelse(d$year <= 2022, surg_old, surg_new)
d$surgery <- factor(ifelse(is.na(surg_code) | surg_code %in% c("Blank(s)", "00", "A000", "99", "A990"), "No/unknown", "Yes"),
                    levels = c("No/unknown", "Yes"))
d$chemotherapy <- factor(ifelse(d[["Chemotherapy recode (yes, no/unk)"]] == "Yes", "Yes", "No/unknown"),
                         levels = c("No/unknown", "Yes"))
rad_raw <- d[["Radiation recode"]]
d$radiation <- factor(ifelse(rad_raw %in% c("None/Unknown", "Refused (1988+)",
                                            "Recommended, unknown if administered"), "No/unknown", "Yes"),
                      levels = c("No/unknown", "Yes"))

d$os_days_raw <- num(d[["Survival Days"]])
d$os_time <- pmax(d$os_days_raw, 0.5) / 365.25
d$os_event <- as.integer(d[["Vital status recode (study cutoff used)"]] == "Dead")
css_raw <- d[["SEER cause-specific death classification"]]
d$css_event <- as.integer(css_raw == "Dead (attributable to this cancer dx)")
d$css_known <- css_raw != "Dead (missing/unknown COD)"
d$single_primary <- d[["Sequence number"]] == "One primary only" &
                    d[["Total number of in situ/malignant tumors for patient"]] == "01"

saveRDS(d, processed_file, compress = "xz")
write.csv(flow, file.path(table_dir, "01_cohort_flow.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(confirmation_flow, file.path(table_dir, "01b_histologic_confirmation_flow.csv"), row.names = FALSE, fileEncoding = "UTF-8")

count_table <- as.data.frame(addmargins(table(d$histology, d$stage)))
names(count_table) <- c("histology", "stage", "n")
write.csv(count_table, file.path(table_dir, "02_stage_counts.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summarize_group <- function(z) {
  c(
    n = nrow(z),
    age_median = median(z$age, na.rm = TRUE),
    age_q1 = unname(quantile(z$age, 0.25, na.rm = TRUE)),
    age_q3 = unname(quantile(z$age, 0.75, na.rm = TRUE)),
    female_pct = mean(z$sex == "Female", na.rm = TRUE) * 100,
    distant_pct = mean(z$stage == "Distant", na.rm = TRUE) * 100,
    surgery_pct = mean(z$surgery == "Yes", na.rm = TRUE) * 100,
    chemotherapy_pct = mean(z$chemotherapy == "Yes", na.rm = TRUE) * 100,
    os_events = sum(z$os_event),
    median_followup_reverse_km = {
      fit <- survfit(Surv(z$os_time, 1 - z$os_event) ~ 1)
      unname(summary(fit)$table["median"])
    }
  )
}
desc <- do.call(rbind, lapply(split(d, d$histology), summarize_group))
desc <- data.frame(histology = rownames(desc), desc, row.names = NULL, check.names = FALSE)
write.csv(desc, file.path(table_dir, "03_basic_description.csv"), row.names = FALSE, fileEncoding = "UTF-8")

km_estimate <- function(z, endpoint = c("OS", "CSS"), times = c(1, 3, 5, 10)) {
  endpoint <- match.arg(endpoint)
  if (endpoint == "CSS") z <- z[z$css_known, , drop = FALSE]
  event <- if (endpoint == "OS") z$os_event else z$css_event
  fit <- survfit(Surv(z$os_time, event) ~ 1, conf.type = "log-log")
  s <- summary(fit, times = times, extend = TRUE)
  med <- unname(summary(fit)$table["median"])
  data.frame(endpoint = endpoint, time = times, n = nrow(z), events = sum(event),
             survival = s$surv, lower = s$lower, upper = s$upper, median_survival = med)
}

km_rows <- list(); idx <- 0
for (ep in c("OS", "CSS")) {
  for (h in levels(d$histology)) {
    for (st in c("All", levels(d$stage))) {
      z <- d[d$histology == h & (st == "All" | d$stage == st), , drop = FALSE]
      if (nrow(z) == 0) next
      out <- km_estimate(z, ep)
      out$histology <- h; out$stage <- st
      idx <- idx + 1; km_rows[[idx]] <- out
    }
  }
}
km <- do.call(rbind, km_rows)
km <- km[, c("endpoint", "histology", "stage", "time", "n", "events", "survival", "lower", "upper", "median_survival")]
write.csv(km, file.path(table_dir, "04_km_survival_estimates.csv"), row.names = FALSE, fileEncoding = "UTF-8")

conditional_estimate <- function(z, endpoint, landmark, horizon) {
  if (endpoint == "CSS") z <- z[z$css_known, , drop = FALSE]
  event <- if (endpoint == "OS") z$os_event else z$css_event
  at_landmark <- z$os_time >= landmark
  z <- z[at_landmark, , drop = FALSE]
  event <- event[at_landmark]
  if (nrow(z) < 2) return(data.frame(n_risk = nrow(z), events = sum(event), estimate = NA, lower = NA, upper = NA))
  residual_time <- z$os_time - landmark
  fit <- survfit(Surv(residual_time, event) ~ 1, conf.type = "log-log")
  s <- summary(fit, times = horizon, extend = TRUE)
  data.frame(n_risk = nrow(z), events = sum(event), estimate = s$surv, lower = s$lower, upper = s$upper)
}

cs_rows <- list(); idx <- 0
for (ep in c("OS", "CSS")) for (h in levels(d$histology)) for (st in c("All", levels(d$stage))) {
  z <- d[d$histology == h & (st == "All" | d$stage == st), , drop = FALSE]
  for (lm in c(0, 1, 2, 3, 5)) for (hz in c(3, 5)) {
    out <- conditional_estimate(z, ep, lm, hz)
    out$endpoint <- ep; out$histology <- h; out$stage <- st; out$landmark <- lm; out$horizon <- hz
    idx <- idx + 1; cs_rows[[idx]] <- out
  }
}
cs <- do.call(rbind, cs_rows)
cs <- cs[, c("endpoint", "histology", "stage", "landmark", "horizon", "n_risk", "events", "estimate", "lower", "upper")]
write.csv(cs, file.path(table_dir, "05_conditional_survival.csv"), row.names = FALSE, fileEncoding = "UTF-8")

piecewise_hazard <- function(z, max_year = 10) {
  do.call(rbind, lapply(0:(max_year - 1), function(start) {
    stop <- start + 1
    exposure <- pmax(0, pmin(z$os_time, stop) - start)
    deaths <- sum(z$os_event == 1 & z$os_time > start & z$os_time <= stop)
    py <- sum(exposure)
    data.frame(interval_start = start, interval_end = stop, n_risk = sum(z$os_time > start),
               deaths = deaths, person_years = py, hazard_per_100py = ifelse(py > 0, 100 * deaths / py, NA))
  }))
}
haz_rows <- list(); idx <- 0
for (h in levels(d$histology)) for (st in c("All", levels(d$stage))) {
  z <- d[d$histology == h & (st == "All" | d$stage == st), , drop = FALSE]
  out <- piecewise_hazard(z)
  out$histology <- h; out$stage <- st
  idx <- idx + 1; haz_rows[[idx]] <- out
}
haz <- do.call(rbind, haz_rows)
haz <- haz[, c("histology", "stage", "interval_start", "interval_end", "n_risk", "deaths", "person_years", "hazard_per_100py")]
write.csv(haz, file.path(table_dir, "06_piecewise_os_hazard.csv"), row.names = FALSE, fileEncoding = "UTF-8")

tidy_cox <- function(fit, model, endpoint, stage = "All", comparator = "8140+8500") {
  co <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  target <- grep("^histologypACC$", rownames(co))
  if (!length(target)) return(NULL)
  data.frame(model = model, endpoint = endpoint, stage = stage, comparator = comparator,
             term = rownames(co)[target], HR = ci[target, "exp(coef)"],
             lower95 = ci[target, "lower .95"], upper95 = ci[target, "upper .95"],
             p = co[target, "Pr(>|z|)"], n = fit$n, events = fit$nevent)
}

fit_cox <- function(z, endpoint = c("OS", "CSS"), adjusted = TRUE, stage_label = "All", comparator = "8140+8500") {
  endpoint <- match.arg(endpoint)
  if (endpoint == "CSS") z <- z[z$css_known, , drop = FALSE]
  z$event_tmp <- if (endpoint == "OS") z$os_event else z$css_event
  rhs <- if (adjusted) "histology + age + sex + race3 + origin + stage + site + era" else "histology"
  if (stage_label != "All") rhs <- if (adjusted) "histology + age + sex + race3 + origin + site + era" else "histology"
  f <- as.formula(paste("Surv(os_time, event_tmp) ~", rhs, "+ cluster(patient_id)"))
  fit <- coxph(f, data = z, ties = "efron", na.action = na.omit)
  tidy_cox(fit, ifelse(adjusted, "Adjusted baseline factors", "Unadjusted"), endpoint, stage_label, comparator)
}

cox_rows <- list(); idx <- 0
for (ep in c("OS", "CSS")) for (adj in c(FALSE, TRUE)) for (st in c("All", "Localized", "Regional", "Distant")) {
  z <- d[st == "All" | d$stage == st, , drop = FALSE]
  out <- fit_cox(z, ep, adj, st)
  if (!is.null(out)) { idx <- idx + 1; cox_rows[[idx]] <- out }
}

for (code in c("8140", "8500")) {
  z <- d[d$hist_code %in% c("8550", code), , drop = FALSE]
  z$histology <- factor(ifelse(z$hist_code == "8550", "pACC", "PDAC"), levels = c("PDAC", "pACC"))
  for (ep in c("OS", "CSS")) {
    out <- fit_cox(z, ep, TRUE, "All", paste0(code, " only"))
    if (!is.null(out)) { idx <- idx + 1; cox_rows[[idx]] <- out }
  }
}
cox <- do.call(rbind, cox_rows)
write.csv(cox, file.path(table_dir, "07_cox_histology_effects.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Automated value gate: sufficient rare-tumor sample, usable conditional-risk sets,
# clinically meaningful adjusted effect, and directionally stable comparator definitions.
pacc_n <- sum(d$histology == "pACC")
pacc_cs3 <- subset(cs, endpoint == "OS" & histology == "pACC" & stage %in% c("Localized", "Regional", "Distant") & landmark == 3 & horizon == 3)
adequate_stage_landmarks <- sum(pacc_cs3$n_risk >= 20)
main_hr <- subset(cox, model == "Adjusted baseline factors" & endpoint == "OS" & stage == "All" & comparator == "8140+8500")$HR[1]
sens_hr <- subset(cox, model == "Adjusted baseline factors" & endpoint == "OS" & stage == "All" & comparator %in% c("8140 only", "8500 only"))$HR
clinically_meaningful <- is.finite(main_hr) && (main_hr <= 0.80 || main_hr >= 1.25)
direction_stable <- length(sens_hr) == 2 && all(sens_hr < 1) == (main_hr < 1)
value_gate <- data.frame(
  criterion = c("Final pACC sample >=300", "At least two stages have >=20 pACC survivors at 3-year landmark",
                "Adjusted OS HR is clinically meaningful (<=0.80 or >=1.25)", "Direction stable for 8140-only and 8500-only comparators"),
  observed = c(as.character(pacc_n), as.character(adequate_stage_landmarks), sprintf("%.3f", main_hr), paste(sprintf("%.3f", sens_hr), collapse = "; ")),
  pass = c(pacc_n >= 300, adequate_stage_landmarks >= 2, clinically_meaningful, direction_stable)
)
write.csv(value_gate, file.path(table_dir, "08_value_gate.csv"), row.names = FALSE, fileEncoding = "UTF-8")

qc_counts <- data.frame(
  item = c(
    "Zero-day survival observations assigned 0.5 day",
    "Records excluded from cancer-specific analyses because cause of death was missing/unknown",
    "PDAC records excluded from cancer-specific analyses because cause of death was missing/unknown",
    "pACC records excluded from cancer-specific analyses because cause of death was missing/unknown"
  ),
  n = c(
    sum(d$os_days_raw == 0, na.rm = TRUE),
    sum(!d$css_known),
    sum(d$histology == "PDAC" & !d$css_known),
    sum(d$histology == "pACC" & !d$css_known)
  )
)
write.csv(qc_counts, file.path(table_dir, "19_endpoint_qc_counts.csv"), row.names = FALSE, fileEncoding = "UTF-8")

capture.output(sessionInfo(), file = file.path(table_dir, "99_session_info.txt"))
message("Core analysis complete. Final cohort n=", nrow(d), "; pACC n=", pacc_n)
