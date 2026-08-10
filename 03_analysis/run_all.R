options(warn = 1)

scripts <- c(
  "00_label_audit.R",
  "01_core_analysis.R",
  "02_extended_analysis.R",
  "06_competing_risk_temporal.R",
  "04_manuscript_tables.R",
  "03_figures.R"
)

for (script in scripts) {
  message("\n===== Running ", script, " =====")
  source(file.path("03_analysis", script), local = new.env(parent = globalenv()))
}

message("\nFull analysis pipeline completed.")
