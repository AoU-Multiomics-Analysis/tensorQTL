# Cis permutations with phenotype and genotype shards

Use `workflows/tensorqtl_cis_permutations_sharded.wdl` on Terra managed
Cromwell. The descriptor uses WDL 1.0. The existing workflows are unchanged.
The complete workflow has not been tested on Terra.

## Task sequence

1. One CPU task localizes the source PLINK2 files. It sorts phenotypes by
   chromosome and position. It puts nearby phenotypes into the same shard.
2. The task extracts the union of the full cis windows for each shard. It
   retains both window boundaries. It excludes gaps between windows and
   samples that are absent from the phenotype file.
3. GPU tasks run cis permutations on the smaller PLINK and phenotype files.
   Each phenotype occurs in one shard. A shard can contain several chromosomes.
4. One CPU task combines the raw permutation results and calculates global
   Storey q-values with tensorQTL and R `qvalue`.

The preparation task downloads the source once. PLINK can read the localized
source several times during extraction. GPU tasks download only their shard.
Overlapping windows can put the same variant in several shards. This is
required to retain complete windows. The manifest records actual genotype
bytes per shard; it does not promise a fixed reduction in total transfer.

## Statistical method

For each phenotype, tensorQTL permutes sample labels and records the strongest
cis association in each permutation. It fits a beta distribution to these
results and returns `pval_beta`. This procedure accounts for correlated
variants in the cis window. It does not use an explicit haplotype count.

After the merge, tensorQTL sends all phenotype-level `pval_beta` values to the
R `qvalue` package. That package estimates the proportion of true null tests
and calculates Storey q-values. Shard-level q-values are not calculated. For
grouped phenotypes, there is one test per group, and the same global step
uses the group-level p-values. See the tensorQTL
[permutation implementation](https://github.com/broadinstitute/tensorqtl/blob/v1.0.10/tensorqtl/cis.py#L577)
and [q-value implementation](https://github.com/broadinstitute/tensorqtl/blob/v1.0.10/tensorqtl/post.py#L16).

Every shard uses the same seed and sample order. This retains the permutation
scheme used in an unsharded call. Small numeric differences can occur between
CPU and GPU execution or between matrix sizes.

## Inputs

Start with `inputs.sharded.example.json`. Replace all bucket paths and the
image tag. The example is a Terra workflow input file. Tasks use named CLI
arguments; they do not use JSON argument wrappers.

| Input | Default | Requirement |
|---|---|---|
| `plink_pgen`, `plink_pvar`, `plink_psam` | Required | One PLINK2 triplet, for one chromosome or the complete genome |
| `phenotype_bed` | Required | Plain or gzip BED; columns `#chr`, `start`, `end`, `phenotype_id`, then samples |
| `covariates` | Required | TSV; covariates as rows, samples as columns |
| `prefix` | Required | Output filename prefix |
| `docker_image` | Required | Published image digest or full commit tag |
| `phenotype_groups` | Absent | Optional two-column TSV without a header; one phenotype ID and group ID per row |
| `cis_window` | `1000000` | Nonnegative distance in bp on each side of the phenotype |
| `max_phenotypes_per_shard` | `250` | Positive target; groups remain complete |
| `max_variants_per_shard` | `100000` | Positive target for the union of cis variants |
| `permutations` | `10000` | Number of permutations, at least 10 |
| `maf_threshold` | `0.0` | In-sample MAF filter; applied during mapping |
| `seed` | `1` | Same seed in every shard |
| `dosages` | `false` | Read PGEN dosage values when true; otherwise read hard calls |
| `fdr` | `0.05` | Global q-value threshold |
| `qvalue_lambda` | Absent | Optional R qvalue tuning parameter |

For BED interval `[start, end)`, extraction uses the one-based inclusive
range `[max(1, start + 1 - cis_window), end + cis_window]`. This matches
tensorQTL's BED conversion. For a TSS at BED `[99, 100)` and a window of
10 bp, retain variant positions 90 through 110. Variant positions 89 and 111
are outside that window. Interval phenotypes use both their start and end.

Standard aliases such as `1` and `chr1` are matched. The input PVAR must not
contain both names for the same chromosome. Nonstandard contig names must
match exactly. PVAR positions must be sorted within each chromosome. Use
uncompressed PVAR files. Split multiallelic variants before mapping. IDs for
selected variants, phenotypes, and samples must be unique and contain no
whitespace. Sample IDs, including leading zeros, are read as text.

All phenotype samples must occur in both the PSAM and covariates. The mapper
aligns covariates by sample ID. Missing or nonfinite phenotype and covariate
values cause an error. Constant phenotypes and phenotypes without cis
variants are listed in `preparation_skipped`.

When a groups file is supplied, it must cover exactly the phenotype IDs in
the BED file. Members of a group stay in one shard. They must have the same
cis variant range on the same chromosome, as required by tensorQTL. A single
group or window can exceed a shard target. The task records that condition
in the manifest; it does not remove variants to meet the target.

## Resources

The preparation defaults are 4 CPUs, 32 GB memory, and 500 GB disk. Set
`prepare_disk_gb` for the source files, all extracted shards, and temporary
files. Overlapping windows increase the required disk space. PLINK receives
half the requested preparation memory; Python also needs memory for variant
positions and phenotype values.

GPU workers default to 16 CPUs, 64 GB memory, 100 GB disk, one NVIDIA L4, and
`g2-standard-16` in `us-central1-c`. Use `map_cpu`, `map_memory_gb`,
`map_disk_gb`, `machine_type`, `gpu_type`, `gpu_count`, and `zones` to change
these settings. The machine and GPU settings must be compatible. Set
`machine_type` to null or an empty string for attached-GPU machine selection.
A worker fails if a GPU is requested but CUDA is unavailable. tensorQTL uses
one GPU per worker; use `gpu_count=1` for GPU runs.

Merge defaults are 8 GB memory and 20 GB disk. Tasks allow two preemptible
retries. Use `preemptible=0` when retries are not wanted.

## Outputs

| Output | Content |
|---|---|
| `cis_qtl` | Combined cis summary, with global `qval` and tensorQTL significance thresholds when available |
| `shard_results` | Raw permutation summaries without q-values |
| `manifest` | Per-shard chromosomes, phenotype/test/variant counts, window size, bytes, and oversized status |
| `preparation_skipped` | Constant phenotypes and phenotypes without cis variants |
| `mapping_skipped` | Phenotypes or group members without valid variants after MAF and monomorphic filters |
| `preparation_log`, `shard_logs`, `merge_log` | Stage logs |

An empty mapping shard returns a header-only result and a skip report. Merge
requires at least two valid tests. It rejects duplicate tests and invalid
beta p-values. R qvalue can fail when the p-value set does not support its
default null-proportion estimate. The task fails with a clear error. An
explicit `qvalue_lambda=0` uses a conservative null proportion of one. No
automatic substitution with BH or another FDR method occurs.

## Image and checks

`envs/Dockerfile.sharded` uses micromamba with pinned Python, PLINK2,
tensorQTL 1.0.10, CUDA-enabled PyTorch 2.5.1, R, `qvalue`, and rpy2. The image
includes `scripts/sharded_cis.py`. The root image Actions workflow builds
on changes to that script and provides a manual trigger. A Dockerfile-only
change requires the manual trigger. Successful runs on `main` publish an
image tagged with the full commit SHA. Set `docker_image` to that tested tag
or its digest. The image has not been built or published during this task.

The root validation Actions workflow checks all three WDL descriptors with
miniwdl and Cromwell's parser. It runs unit and command tests, then a real
PLINK2/tensorQTL/R CPU smoke test. The smoke test compares raw permutations
and global q-values from sharded and unsharded inputs. It checks cis
boundaries, extracted genotypes, sample IDs, groups, an empty shard, and
mixed chromosome shards. It also checks fractional dosage retention and
the dosage mapping path. The image workflow runs that smoke test in its
candidate container before publication. These Actions checks have been
configured; a local check does not prove that they passed on GitHub.

Preparation returns five explicit `Array[File]` outputs from `glob()`:
PGEN, PVAR, PSAM, phenotype BED, and phenotype groups files. Each file type
uses the same fixed-width shard basenames, such as `shard_000000`.
WDL 1.0 specifies Bash filename order for `glob()`. The indexed scatter
selects the same position from each array to keep each shard's files together.
The task-local `shards.json` is audit metadata only. It does not declare
shard files for transfer.

Files remain typed `File` values through task inputs. Merge receives an
explicit `Array[File]`; its command writes the localized paths to a local
newline list. A static AST check rejects file writes at workflow
scope. Command tests simulate cloud-to-local path rewriting for every new
task, but do not prove Terra localization or GPU execution.

## Terra output transfer error

An earlier WDL used `read_json("shards.json")` to return shard files in
struct fields. That output declaration can cause this error when a mapping
task downloads its inputs:

```text
ERROR: (gcloud.storage.cp) The following URLs matched no objects or files:
```

Cromwell plans output transfers before the preparation command runs. At
that time, it cannot read `shards.json` to find the shard file paths.
Declaring its fields as `File` does not add those files to the transfer plan.
The explicit `glob()` outputs let Cromwell plan transfers for all shard
files. The mapping input download uses Cromwell's own `gcloud storage cp`
command. See the [Cromwell Batch output collection code](https://github.com/broadinstitute/cromwell/blob/92/supportedBackends/google/batch/src/main/scala/cromwell/backend/google/batch/actors/GcpBatchAsyncBackendJobExecutionActor.scala#L844-L868).

This correction changes only the WDL. The original image and preparation
script can be used. No image rebuild is required. Refresh the workflow in
Dockstore and Terra to use the corrected WDL. Run preparation again. Do not
reuse the cached preparation result from the faulty WDL.

The complete corrected workflow has not been tested on Terra. Local checks
and GitHub Actions checks do not prove Terra output transfer behavior.

Local checks:

```bash
miniwdl check workflows/tensorqtl_cis_permutations_sharded.wdl
python -m unittest discover -s tests
```

Do not claim Terra validation until the complete workflow runs on Terra.
No cloud job was submitted during this task.

## Terra merge file-list error

The merge task can fail with this error when it uses an earlier WDL:

```text
script: line 43: gs://.../call-merge_cis/write_lines_HASH.tmp: No such file or directory
```

The earlier command used `sub(write_lines(results), ...)` for the result
list path. `write_lines()` returns a `File`, but `sub()` returns a `String`.
The command then used the resulting cloud URI as a local filename. Bash
cannot open that URI. See the [Cromwell `sub()` implementation](https://github.com/broadinstitute/cromwell/blob/92/wdl/transforms/new-base/src/main/scala/wdl/transforms/base/linking/expression/values/EngineFunctionEvaluators.scala#L804-L820).

The corrected task keeps `results` as `Array[File]` until command rendering.
The command writes those localized paths to `shard_results.list` with a
quoted heredoc. It rejects an empty array, unresolved cloud URIs, and files
that it cannot read before it starts the merge script.

Refresh the workflow in Dockstore and Terra, then submit it again with the
corrected WDL. Use the same Docker image; this change does not require an
image rebuild. Keep call caching enabled for this merge error. Terra can
reuse successful preparation and mapping calls when matching cache entries
and their output files are available. Cache reuse is not guaranteed. The
earlier output-transfer error requires the separate preparation steps above.

Command tests simulate cloud-to-local path changes. The existing CPU smoke
test checks the merge with PLINK2, tensorQTL, and R. These checks do not prove
Terra localization. The complete corrected workflow has not been tested on
Terra. No local Docker image was built, and no cloud job was submitted for
this correction.
