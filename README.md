# Conditional survival of pACC versus PDAC by SEER Summary Stage

This repository contains the R analysis code for a population-based comparison of pancreatic acinar cell carcinoma (pACC) and pancreatic ductal adenocarcinoma (PDAC) using SEER 17 registry data from 2004 through 2023.

The analysis estimates overall survival, conventional cancer-specific survival, competing-risk cumulative incidence, Fine-Gray subdistribution hazards, stage-specific associations, interval-specific hazards, and landmark conditional survival.

## Revision analysis package

The revision scripts and complete usage notes are in [03_analysis/revision](03_analysis/revision/README.md). They add follow-up-restricted and 8500-only conditional OS, included/excluded comparisons, PDAC calendar trends, era-specific conditional OS, pACC landmark risk-set composition and independent raw-export validation. The extraction used SEER*Stat 9.0.43.0 with follow-up through December 31, 2023.

After running the original pipeline, run:

```text
Rscript --vanilla 03_analysis/run_revision.R
python 03_analysis/revision/validate_independent.py
python 03_analysis/revision/validate_remaining.py
```

Outputs are written under the ignored `04_results/revision/` directory. The original pipeline remains available for reproduction; the new conditional-OS grid and calendar models supersede the corresponding historical outputs. Overall pACC conditional survival improved in common diagnosis cohorts, but the data do not support a general pACC-PDAC convergence conclusion.

## Data access and restrictions

Patient-level SEER data are not included in this repository and must not be redistributed. Researchers must obtain access independently from the [SEER Program](https://seer.cancer.gov/data-software/) and accept the applicable SEER Research Data Use Agreement.

After exporting the required case-listing variables in SEER*Stat, save the tab-delimited file locally as:

```text
01_data_extraction/pacc_pdac_seer17_2004_2023_raw.txt
```

This path is excluded by `.gitignore`. Do not commit the raw export, SEER dictionary, SEER*Stat session, processed patient-level data, or files containing patient identifiers.

The required variables and their analytic roles are listed in `01_data_extraction/variable_contract.csv`.

## Cohort definitions

- Pancreatic primary sites: ICD-O-3 C25.0-C25.9
- pACC: ICD-O-3 histology 8550/3
- PDAC comparator: ICD-O-3 histology 8140/3 or 8500/3
- Diagnosis years: 2004-2023
- Age: 18 years or older
- Positive histologic confirmation required
- Autopsy-only and death-certificate-only records excluded

Mixed acinar-neuroendocrine carcinoma (8154), acinar cystadenocarcinoma (8551), and mixed acinar-ductal carcinoma (8552) are not part of the primary pACC definition.

## Software

The verified analysis environment used:

- R 4.6.1
- survival 3.8-6
- ggplot2 4.0.3
- patchwork 1.3.2

Check package availability with:

```bash
Rscript 03_analysis/check_dependencies.R
```

## Reproducing the analysis

Run all commands from the repository root:

```bash
Rscript 03_analysis/run_all.R
```

The pipeline runs in this order:

1. `00_label_audit.R`: verifies raw SEER labels and category frequencies.
2. `01_core_analysis.R`: constructs the cohort and performs the primary survival analyses.
3. `02_extended_analysis.R`: performs sensitivity, interaction, proportional-hazards, and interval-specific analyses.
4. `06_competing_risk_temporal.R`: performs Fine-Gray, cumulative-incidence, conditional competing-risk, and temporal analyses.
5. `04_manuscript_tables.R`: creates aggregate tables.
6. `03_figures.R`: creates publication figures and aggregate figure source files.

Generated patient-level intermediate data are written to `02_processed_data/`. Aggregate tables and figures are written to `04_results/`. Both directories are ignored to prevent accidental publication from a local analysis run; aggregate outputs intended for publication should be reviewed separately before deposit in an archival repository.

## Repository contents

```text
01_data_extraction/
  README.md
  variable_contract.csv
03_analysis/
  00_label_audit.R
  01_core_analysis.R
  02_extended_analysis.R
  03_figures.R
  04_manuscript_tables.R
  06_competing_risk_temporal.R
  check_dependencies.R
  install_figure_packages.R
  run_all.R
REPOSITORY_AUDIT.md
```

## Scope of inference

The analyses estimate prognostic associations and time-updated survival. Recorded treatment variables are used descriptively and in exploratory prognostic adjustment, not to estimate causal treatment effects.

## Citation and licence

A formal citation and software licence will be added after the author list, repository record, and journal submission details are finalized. Until a licence is added, the code remains copyrighted and no reuse licence is granted.
