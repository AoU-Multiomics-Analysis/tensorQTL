# Pre-shard Variant QC Implementation Plan

> **For agentic workers:** Use the existing checkout. Implement independent file groups in parallel. Verify all changes before completion.

**Goal:** Apply MAC >= 10 and missing fraction <= 0.05 before cis shard planning.

**Architecture:** Stream PGEN batches for the phenotype samples before imputation. Use the mapping hard-call/dosage mode, write retained positions and IDs, and plan and extract only retained variants.

**Tech Stack:** Python, NumPy, pgenlib 0.93.0, PLINK2 2.00a5.12, WDL 1.0.

**Spec:** `docs/superpowers/specs/2026-11-09-pre-shard-variant-qc-design.md`

## Global constraints

- Defaults: `min_mac=10`, `max_missing_fraction=0.05`; both configurable.
- Use analysis samples and observed diploid-coded calls/dosages before imputation.
- No local Docker build or cloud submission. Preserve Terra File types and script-only image triggers.

## Review focus

- Samples absent from the phenotype BED cannot affect QC.
- Dosage availability and dosage MAC can differ from hard-call values.
- Equality at each threshold passes; missingness uses all analysis samples.
- Empty windows are reported before scatter; an empty analysis fails clearly.
- QC and extraction must retain exactly the same variants at both window boundaries.

## Task 1: QC and shard preparation

- [x] Add failing tests for CLI limits; add real PGEN boundary/sample-subset cases.
- [x] Add bounded PGEN QC, retained list/index, summary, and empty-window handling in `cis_permutations/scripts/sharded_cis.py`.
- [x] Run focused tests and verify retained QC counts. PLINK extraction is covered by the configured CI smoke.

## Task 2: Terra inputs and documentation

- [x] Add WDL `Int min_mac=10`, `Float max_missing_fraction=0.05`, and preparation `dosages` input.
- [x] Add `File variant_qc` output, named arguments, cloud-localization tests, and input example/documentation.
- [x] Run WDL validation and rendered command tests.

## Task 3: CI smoke and final verification

- [x] Extend the real PLINK/tensorQTL smoke fixture for subset MAC, missingness, dosage mode, and preparation skips.
- [x] Run the 54-test cis suite with real PGEN dependencies and review the diff. A fresh reviewer found no material defects.
- [x] Document local checks and distinguish configured CI smoke from untested Terra execution.

## Verification notes

The 54 cis tests, 11 nominal tests, and 28 trans tests pass. WDL validation passes with
existing runtime and optional-input warnings. The full PLINK/tensorQTL/R
smoke remains for GitHub Actions. No local Docker image was built. No
workflow was submitted to Terra. Changes remain in the existing checkout.
