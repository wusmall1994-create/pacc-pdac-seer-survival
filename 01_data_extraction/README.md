# Local SEER data directory

Place the authorized SEER*Stat case-listing export in this directory using the filename:

```text
pacc_pdac_seer17_2004_2023_raw.txt
```

The raw export is patient-level third-party restricted data. It is intentionally absent from this repository and excluded by `.gitignore`.

Do not commit:

- SEER case-listing exports (`.txt` or `.gz`)
- export dictionaries (`.dic`)
- SEER*Stat sessions or matrices (`.sl` or `.slm`)
- files containing `Patient ID` values

Use `variable_contract.csv` to configure and verify the SEER*Stat export.
