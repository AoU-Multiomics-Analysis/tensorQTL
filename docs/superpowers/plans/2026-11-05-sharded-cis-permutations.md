# Sharded cis permutations implementation plan

Use native execution in the current checkout. The user approved the design and
delegated the remaining implementation choices. No further planning approval
is needed. Design: `../specs/2026-11-05-sharded-cis-permutations-design.md`.

- [x] Add failing preparation tests in `cis_permutations/tests/test_sharding.py`.
   Cover exact cis boundaries, chromosome aliases, window unions, group
   integrity, mixed chromosomes, size limits, and invalid inputs.
- [x] Add `scripts/sharded_cis.py` with typed named CLI arguments for preparation,
   mapping, and merge. Use PLINK2 extraction without loading source genotypes.
   Call tensorQTL's API for raw permutations and global Storey qvalues.
- [x] Add the WDL 1.0 descriptor and command tests. Exercise `File` localization
   for preparation, mapping, and merge. Reject workflow-scope file writes.
- [x] Add a pinned micromamba image and root GitHub Actions workflow. Build only
   in Actions; filter image builds to changes in the included script. Add real
   PLINK2, tensorQTL, and R qvalue smoke tests, including sharded versus whole
   input comparison. Validate with miniwdl and Cromwell's parser.
- [x] Add a Terra input example, documentation, and Dockstore entry. Run existing
   and new tests, inspect the diff, and request an independent review. Report
   the local evidence and the untested Terra/GPU/Actions limits.

## Verification record

- Baseline: 18 existing cis permutation tests passed before implementation.
- Preparation and WDL tests failed before the corresponding code was added.
- Independent design check found no material unanswered questions.
- Independent code review found two important defects: an IID-only PLINK keep
  file lost source family IDs, and pandas could convert text IDs to NA values
  during merge. Each received a failing regression test and a passing fix.
- Local suites: 38 cis permutation, 11 nominal cis, and 28 trans tests passed.
- miniwdl accepts the new descriptor. Expected lint warnings concern the
  Terra machine runtime key and the optional machine input with a default.
- actionlint accepts the two changed root Actions workflows. The example
  Terra input file passes WDL input type checks. JSON/YAML parsing and
  `git diff --check` pass.
- A Linux micromamba dry-run resolves the pinned environment. R 4.3.3 and
  qvalue 2.34.0 are used with rpy2 3.5.11; the initial R 4.4 combination was
  not compatible with that rpy2 build.
- Real PLINK2/tensorQTL/R and container smoke tests are configured in Actions,
  including hard calls, fractional dosages, grouped tests, and empty shards.
  They have not been run during this task.
- No local Docker build, image publication, GitHub Actions run, GPU run, or
  Terra submission occurred. Local Cromwell parser validation was unavailable
  because no Java runtime is installed. It is configured in Actions.
- Changes remain in the existing checkout; no commit, push, or merge was made.
