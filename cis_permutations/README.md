# tensorQTL cis Permutations

This repository contains WDL workflows for running **cis-QTL mapping with permutations** using [tensorQTL](https://github.com/broadinstitute/tensorqtl). Two workflows are provided:

| Workflow | WDL file | Use case |
|---|---|---|
| `tensorqtl_cis_permutations_workflow` | `workflows/tensorqtl_cis_permutations.wdl` | Standard cis-eQTL mapping (expression QTLs) |
| `tensorqtl_cis_permutations_workflow` | `workflows/tensorqtl_cis_permutations_sQTL.wdl` | cis-sQTL mapping (splicing QTLs) – requires a phenotype groups file |

The workflows are registered on [Dockstore](https://dockstore.org/) and intended to be executed on a platform that supports WDL (e.g. Terra / Cromwell). Each task runs tensorQTL in `cis` permutation mode and produces a genome-wide cis-QTL summary table.

---

## Workflows

### `tensorqtl_cis_permutations.wdl` – eQTL workflow

Runs tensorQTL cis permutation testing for expression QTLs. The `phenotype_groups` file is optional and may be omitted when all phenotypes are independent.

### `tensorqtl_cis_permutations_sQTL.wdl` – sQTL workflow

Runs tensorQTL cis permutation testing for splicing QTLs. A `phenotype_groups` file is **required** because sQTL analyses group multiple intron clusters (or splice junctions) by gene.

---

## Inputs

### Required inputs

| Input | Type | Description |
|---|---|---|
| `plink_pgen` | `File` | PLINK2 `.pgen` file containing genotype data |
| `plink_pvar` | `File` | PLINK2 `.pvar` file containing variant information |
| `plink_psam` | `File` | PLINK2 `.psam` file containing sample information |
| `phenotype_bed` | `File` | bgzipped and tabix-indexed BED file of phenotype values (samples as columns, one phenotype per row) |
| `covariates` | `File` | Tab-separated covariates file (samples as columns, one covariate per row; first column is the covariate name) |
| `prefix` | `String` | Output file prefix (e.g. `my_cohort.chr1`) |
| `memory` | `Int` | Memory to allocate in GB |
| `disk_space` | `Int` | Disk space to allocate in GB |
| `num_threads` | `Int` | Number of CPU threads |
| `num_gpus` | `Int` | Number of GPUs; the eQTL default machine includes one L4 |
| `num_preempt` | `Int` | Number of preemptible retries |

### Optional inputs

| Input | Type | Default | Description |
|---|---|---|---|
| `phenotype_groups` | `File` | – | Two-column file grouping phenotypes by gene; **required for sQTL**, optional for eQTL (needed when multiple phenotypes share a gene) |
| `fdr` | `Float` | `0.05` | FDR threshold for significant associations |
| `qvalue_lambda` | `Float` | – | Lambda parameter for q-value calculation |
| `pval_threshold` | `Float` | – | p-value threshold to filter nominal results |
| `seed` | `Int` | – | Random seed for reproducibility |
| `flags` | `String` | – | Additional tensorQTL arguments; the eQTL workflow parses shell-style quotes without executing shell expressions |

---

## Outputs

| Output | Description |
|---|---|
| `cis_qtl` | `{prefix}.cis_qtl.txt.gz` – Compressed table of cis-QTL results (one row per phenotype, top variant and permutation-derived p-values) |
| `log` | `{prefix}.tensorQTL.cis.log` – tensorQTL run log |

---

## Data preparation

### Genotype data (PLINK2 format)

tensorQTL expects genotypes in **PLINK2 binary format** (`.pgen` / `.pvar` / `.psam`). To convert from VCF:

```bash
plink2 --vcf input.vcf.gz \
       --make-pgen \
       --out genotypes
```

Recommended pre-processing steps:
- Filter to samples present in your phenotype and covariate files.
- Remove variants with missing genotype rate > 5% and MAF < 0.01.
- Split multi-allelic sites and left-align indels.
- Optionally split by chromosome for parallelism.

### Phenotype file (BED format)

The phenotype file must be a **bgzipped, tabix-indexed BED file** where:
- The first four columns are `#chr`, `start`, `end`, `phenotype_id`.
- Subsequent columns are sample IDs matching those in the PLINK `.psam` file.
- Rows are phenotypes (genes for eQTL; splice junctions / intron clusters for sQTL).
- Values should be normalized (e.g. inverse-normal transformed).

Create and index the file:

```bash
bgzip phenotypes.bed
tabix -p bed phenotypes.bed.gz
```

### Covariates file

A tab-separated text file with:
- First column: covariate name.
- Remaining columns: one column per sample (in the same order as the phenotype BED file).
- Typical covariates include genotype PCs, expression/splicing PCs, and known technical factors.

Example:
```
ID    SAMPLE1  SAMPLE2  ...
PC1   0.021    -0.013   ...
PC2   0.004     0.031   ...
```

### Phenotype groups file (sQTL / optional eQTL)

A two-column, tab-separated file (no header) mapping each phenotype ID to its group (gene):

```
ENSG00000000003:1:100  ENSG00000000003
ENSG00000000003:2:200  ENSG00000000003
ENSG00000000005:1:300  ENSG00000000005
```

This file ensures that multiple phenotypes belonging to the same gene share a single permutation test.

---

## Scripts

The `scripts/` directory is reserved for helper scripts used to prepare or post-process data for these workflows. Currently no scripts are included; contributions are welcome.

---

## Runtime environment

Tasks run inside a Docker container:

```
gcr.io/broad-cga-francois-gtex/tensorqtl:latest
```

The eQTL workflow uses WDL 1.0 and targets Terra managed Cromwell. Its
runtime defaults to `g2-standard-16`, NVIDIA L4, and `us-central1-c`.
Existing input keys retain the prefix
`tensorqtl_cis_permutations_workflow.tensorqtl_cis_permutations.`.

| eQTL input | Type | Default |
|---|---|---|
| `machine_type` | `String?` | `g2-standard-16` |
| `gpu_type` | `String` | `nvidia-l4` |
| `zones` | `Array[String]` | `["us-central1-c"]` |

A predefined machine type fixes the CPU count and memory. The GPU type and
count must be compatible with that machine.
Use `num_threads=16`, `memory=64` (GB), and `num_gpus=1` with the default machine.
To select another predefined GPU machine, set `machine_type` to a compatible
type such as `a2-highgpu-1g` (one A100), and choose zones that support it.
Set `gpu_type="nvidia-tesla-a100"` and `num_gpus=1` for that example.

To request an attached GPU such as a T4, set `machine_type` to `null`
or an empty string. The runtime expression converts `null` to an empty
`String`. Google Batch then selects a machine from `num_threads`, `memory`,
`gpu_type`, and `num_gpus`. A direct `String?` runtime value causes
`Invalid 'predefinedMachineType': String value required but got String?`.
Add these entries to the existing Terra workflow inputs:

```json
{
  "tensorqtl_cis_permutations_workflow.tensorqtl_cis_permutations.machine_type": null,
  "tensorqtl_cis_permutations_workflow.tensorqtl_cis_permutations.gpu_type": "nvidia-tesla-t4",
  "tensorqtl_cis_permutations_workflow.tensorqtl_cis_permutations.num_gpus": 1,
  "tensorqtl_cis_permutations_workflow.tensorqtl_cis_permutations.zones": ["us-central1-c"]
}
```

The selected GPU must have sufficient regional quota and available capacity.
The global GPU quota can also apply. Changing the GPU type can avoid the L4
quota limit, but it does not guarantee that the job will start.
See [Cromwell runtime attributes](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#predefinedmachinetype-alpha)
and [Google Cloud GPU quotas](https://docs.cloud.google.com/compute/resource-usage#gpu_quota).

The eQTL task logs validation and mapping stages. Input files remain typed
`File` values until command rendering. The task checks localized paths and
links all three PLINK files to one basename. Optional arguments are passed
only when present. The output prefix must start with a letter or digit and
contain only letters, digits, dots, underscores, or hyphens.

The sQTL entry point uses the eQTL task and its runtime defaults.
Its required groups input is `tensorqtl_cis_permutations_workflow.phenotype_groups`.
Other task input keys keep the prefix `tensorqtl_cis_permutations_workflow.tensorqtl_cis_permutations.`.
The file requires WDL 1.0. Keep its imported eQTL descriptor in the same directory.

The copied Dockerfile is an empty source placeholder. Its nested Docker Actions file is inactive.

---

## Dockstore registration

Both WDL workflows are registered on Dockstore via `.dockstore.yml`:

| Name | Primary descriptor |
|---|---|
| `tensorQTL_permutation` | `workflows/tensorqtl_cis_permutations.wdl` |
| `tensorQTL_permutation_sQTL` | `workflows/tensorqtl_cis_permutations_sQTL.wdl` |


## Validation

The root GitHub Actions workflow validates both WDL files and runs command tests. Tests cover L4
defaults, machine/GPU/zone overrides, T4 selection without a predefined machine,
optional arguments, safe path quoting, and unresolved cloud URI rejection.
An AST check rejects workflow-scope file-writing functions. These checks do
not build or publish a Docker image.

Local checks:

```bash
miniwdl check workflows/tensorqtl_cis_permutations.wdl
miniwdl check workflows/tensorqtl_cis_permutations_sQTL.wdl
python -m unittest discover -s tests
```

Command tests use a tensorQTL test double. They do not validate the tensorQTL
engine, the production container, GPU execution, or Terra localization.
The complete workflow has not been tested on Terra.
