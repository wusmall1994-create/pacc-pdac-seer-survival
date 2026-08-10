# Repository release audit

## Included

- R scripts required to reconstruct the analytic cohort and reproduce all statistical analyses, tables, and figures
- A variable-level SEER extraction contract without patient values
- Reproduction and data-access instructions
- Dependency versions used for the verified analysis

## Excluded

- Patient-level SEER exports and compressed copies
- SEER export dictionaries, sessions, and matrices
- Processed patient-level R objects
- Patient identifiers or record-level extracts
- Manuscript drafts, Word documents, references, review notes, and rendered manuscript pages
- Local virtual environments, caches, logs, temporary files, and machine-specific paths
- Personal names, email addresses, phone numbers, postal addresses, institutional directory paths, and account credentials

## Automated checks before publication

The release directory is scanned for:

- common email, phone, credential, and absolute-path patterns
- author names and local user or institution identifiers
- manuscript and submission file extensions
- restricted SEER data extensions
- files exceeding normal source-code size expectations

The `.gitignore` file provides a second barrier against accidental addition of restricted data during local reproduction.

## Remaining author decisions

- Choose a software licence.
- Add a formal repository citation after the author list is finalized.
- Archive a tagged release in a DOI-issuing repository if required by the target journal.
