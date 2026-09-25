# Revision analyses

Run from the repository root after the original pipeline:

```text
Rscript --vanilla 03_analysis/run_all.R
Rscript --vanilla 03_analysis/run_revision.R
python 03_analysis/revision/validate_independent.py
python 03_analysis/revision/validate_remaining.py
```

The revision runner accepts an alternative project root and output directory. The two analysis R scripts and Python validators also accept PROJECT_ROOT OUTPUT_DIRECTORY. Defaults use the current directory and `04_results/revision/conditional_OS` or `04_results/revision/selection_and_calendar`. The figure script accepts its output directory as its first argument. R requires survival; Python validators use the standard library. Figure exports use R graphics, including TIFF support.

## Inputs and outputs

Required private inputs are the raw export in `01_data_extraction` and the RDS created by the original pipeline in `02_processed_data`. The first analysis also audits `04_results/tables/05_conditional_survival.csv` produced by the original pipeline. Saved dictionary/session files are optional provenance inputs; they are not needed to calculate results and must not be committed. Input hashes record only files actually present, using project-relative paths.

- `analysis_revision.R`: pooled and 8500-only conditional OS at landmarks 0, 1, 2, 3, 5 years, prediction horizons 3 and 5 years, overall and by SEER Summary Stage. Full, landmark-matched and common diagnosis-window analyses include risk sets, horizon-end support, deaths, censoring and log-log 95% confidence intervals. No extrapolation in the revised grid.
- `analysis_remaining.R`: included/excluded comparisons at both the confirmation-filter step and otherwise-eligible population; flow reconciliation; calendar Cox models for pACC, pooled PDAC and 8500-only PDAC with patient-clustered robust variance; era-specific conditional OS with adequate potential follow-up.
- `validate_independent.py`: independently reconstructs all new conditional OS estimates directly from raw TXT without R objects or survival libraries.
- `validate_remaining.py`: independently checks selection counts, means-based SMDs, medians, and all era-specific conditional OS estimates and risk sets. It does not replicate Cox fitting.
- `riskset_composition.R`: pACC landmark risk-set composition by age at diagnosis, diagnosis year and SEER Summary Stage for the full and common diagnosis cohorts.
- `flow_diagram.R`: audited sequential flow diagram, PDF/SVG/PNG/600-dpi TIFF. Counts are fixed to this cohort and must be updated if eligibility changes.

## Definitions and interpretation

SEER*Stat 9.0.43.0; SEER 17 Nov 2025 submission; administrative cutoff 2023-12-31. Registry-coded pACC with histologic confirmation does not imply central pathology review. Stage means SEER Summary Stage, not AJCC TNM stage.

Matched windows end in 2023 minus landmark minus horizon. Common windows end in 2015 for the three-year horizon and 2013 for the five-year horizon. Cohorts are selected by diagnosis year, never by actually surviving eight or ten years. Numbers at landmarks are records alive and still observed.

Overall pACC conditional OS improved across landmarks in common diagnosis cohorts. Between-histology differences depended on comparator definition, cohort and horizon; there is no general convergence or equivalence conclusion. Calendar models describe associations, not treatment effects. Most confirmation-excluded records had positive cytology without positive histology; the bias direction cannot be assigned from the mixed observed characteristics.

## Changes from the original analysis

The original pipeline is retained for reproducibility. The revision grid supersedes its historical conditional OS grid: five unknown-stage historical cells had extrapolated beyond observed support; original published overall/known-stage OS cells did not. Updated calendar models supersede the original pACC temporal variance calculation: its year HR is unchanged, with robust P=0.075 replacing P=0.069. For revised Table S4, use `deaths_within_horizon` from the new five-year grid, not all deaths over subsequent follow-up.

All 600 conditional-OS rows (544 supported) and 192 selection/era checks were independently validated on the study export. No patient-level or aggregate study outputs are deposited in this code repository; aggregate results accompany the revision as Supplementary Data.
