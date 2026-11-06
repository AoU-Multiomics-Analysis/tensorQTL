# Sharded cis permutation workflow

## Intent and approved design

Reduce genotype transfers to tensorQTL workers. Assign each phenotype to one
shard. Extract all variants in its full cis window before the scatter.
The user approved coordinate-based shards and left the remaining details to
the implementer. Shards may contain more than one chromosome. Chromosome
boundaries do not limit a shard, but each cis window stays on its chromosome.

## Constraints

- Add a WDL 1.0 entry point; preserve the existing entry points.
- Target Terra managed Cromwell. Keep all input paths as `File` values until
  command rendering. Do not create files at workflow scope.
- Use typed inputs and named CLI arguments. Create file lists in task scope.
- Log validation, preparation, mapping, and merge stages.
- Use a micromamba image with pinned packages from conda-forge and bioconda.
- Do not build Docker locally or submit Terra jobs. Add GitHub Actions checks.
- Rebuild the image automatically only when its included script changes.

## Approaches

1. Coordinate-based packing with phenotype and variant limits (selected).
   Nearby windows share genotype data. Distant gaps do not enter the extracted
   data. Two limits control task size; an indivisible group can exceed a limit.
2. Fixed phenotype row counts. This is simple, but gives less control over
   genotype size and can put distant windows in the same shard.
3. Whole chromosomes. This limits the number of downloads, but large
   chromosomes can produce large workers. The user removed this requirement.

## Components and interfaces

Add `cis_permutations/workflows/tensorqtl_cis_permutations_sharded.wdl` and one
image-owned Python CLI, `cis_permutations/scripts/sharded_cis.py`.

The workflow accepts one PLINK2 triplet (`plink_pgen`, `plink_pvar`,
`plink_psam`), a plain or gzip BED phenotype file, covariates, an optional
two-column groups file, `prefix`, `cis_window`, `max_phenotypes_per_shard`,
`max_variants_per_shard`, permutation count, MAF threshold, seed, FDR, optional
qvalue lambda, image, and separate CPU/GPU resource inputs. One source can
contain a chromosome or the complete genome. The caller can run this workflow
separately for independent data sets, but genome-wide FDR needs one merge.

1. **Prepare, CPU:** Localize the source once. Validate files, sample IDs,
   coordinates, and group membership. Read variant positions without loading
   the genotype matrix. Normalize unambiguous standard chromosome aliases.
   Sort phenotypes by chromosome and coordinate. Pack complete phenotype
   groups using phenotype count and the count of variants in the window union.
   Extract ranges with real PLINK2; preserve dosages and sample order. Select
   phenotype samples in PLINK to reduce transfer further. Emit one phenotype
   BED, groups file, and PLINK2 triplet for each shard. Empty groups files are
   valid for ungrouped runs. Output typed arrays and a TSV audit manifest.
2. **Map, GPU:** Localize only that shard's files plus covariates. Call
   `cis.map_cis` directly. Keep the same seed in each shard, so permutations
   match an unsharded call. Do not calculate per-shard q-values. Preserve sample
   IDs as text. Report skipped phenotypes or groups. Return full precision
   permutation summaries and logs.
3. **Merge, CPU:** Read explicit `Array[File]` outputs through a task-local
   newline list. Reject duplicate tests. Run tensorQTL's `calculate_qvalues`
   once on all valid tests. Use Storey's R qvalue implementation, not BH.
   Fail clearly if qvalue cannot estimate FDR; allow an explicit lambda.

## Correctness rules

For BED coordinates `(start, end)`, use the inclusive PLINK interval
`[max(1, start + 1 - cis_window), end + cis_window]`. This matches tensorQTL's
conversion from zero-based BED start to one-based positions. Merge overlapping
intervals within each chromosome. Retain both cis boundaries. A variant may
occur in several shards when their windows overlap; each phenotype occurs in
one shard only. The mapper still applies the original per-phenotype window.

Keep group members together. Before mapping, require the same chromosome and
cis variant range for all members, as tensorQTL requires. Groups do not need
identical phenotype coordinates if their selected variant range is identical.
Reject ambiguous chromosome aliases, duplicate phenotype/sample IDs, malformed
coordinates, nonfinite values, and missing samples. Report phenotypes with no
cis variants or constant values. Fail if no valid tests remain.

A single large cis window or group is indivisible. Record when it exceeds a
requested shard limit. Limits are targets, not permission to remove variants.

## Outputs and verification

Return the combined cis table with global q-values, raw shard tables, mapping
logs, skip reports, and the preparation manifest. The manifest records counts,
cis window, interval counts, actual file sizes, and oversized shards.

Test interval boundaries, overlapping windows, interval phenotypes, chromosome
aliases, mixed chromosome shards, group integrity, bad inputs, no-variant
phenotypes, safe quoting, and cloud-to-local rewriting in every new task.
Statically reject workflow-scope file-writing functions. In GitHub Actions,
run real PLINK2 and CPU tensorQTL and compare sharded results with one unsharded
call using the same seed. Check global qvalues with the real R package. Build
and smoke-test the container in Actions. No complete Terra run is claimed.

## Builder decisions

Initial defaults: cis window 1,000,000 bp, 250 phenotypes per shard, 100,000
variants per shard, 10,000 permutations, MAF 0, seed 1, FDR 0.05, L4 GPU workers.
Use the existing checkout. Leave existing unsharded behavior unchanged.
PLINK1 input and automated cloud submission are outside this change.
