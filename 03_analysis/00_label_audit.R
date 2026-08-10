options(stringsAsFactors = FALSE, width = 180)

input_file <- file.path("01_data_extraction", "pacc_pdac_seer17_2004_2023_raw.txt")
output_file <- file.path("04_results", "tables", "00_label_audit.txt")

d <- read.delim(
  input_file,
  sep = "\t",
  quote = "\"",
  na.strings = "NA",
  check.names = FALSE,
  colClasses = "character"
)

keys <- c(
  "Diagnostic Confirmation",
  "Survival months flag",
  "SEER cause-specific death classification",
  "SEER other cause of death classification",
  "Vital status recode (study cutoff used)",
  "Combined Summary Stage with Expanded Regional Codes (2004+)",
  "Sequence number",
  "Record number recode",
  "Total number of in situ/malignant tumors for patient",
  "Type of Reporting Source",
  "Histologic Type ICD-O-3",
  "Behavior code ICD-O-3",
  "Primary Site - labeled",
  "Grade Recode (thru 2017)",
  "Derived Summary Grade 2018 (2018+)",
  "RX Summ--Surg Prim Site (1998-2022)",
  "RX Summ--Surg Prim Site 2023 (2023+)",
  "Chemotherapy recode (yes, no/unk)",
  "Radiation recode"
)

sink(output_file)
cat("Dimensions:", nrow(d), "rows x", ncol(d), "columns\n")
for (k in keys) {
  cat("\n### ", k, "\n", sep = "")
  print(sort(table(d[[k]], useNA = "ifany"), decreasing = TRUE))
}
cat("\n### Histology by diagnostic confirmation\n")
print(addmargins(table(d[["Histologic Type ICD-O-3"]], d[["Diagnostic Confirmation"]], useNA = "ifany")))
cat("\n### Histology by summary stage\n")
print(addmargins(table(d[["Histologic Type ICD-O-3"]], d[["Combined Summary Stage with Expanded Regional Codes (2004+)"]], useNA = "ifany")))
cat("\n### Histology by zero-day survival\n")
print(addmargins(table(d[["Histologic Type ICD-O-3"]], d[["Survival Days"]] == "0", useNA = "ifany")))
cat("\n### Patient ID multiplicity\n")
print(table(table(d[["Patient ID"]])))
sink()

cat("Wrote", output_file, "\n")
