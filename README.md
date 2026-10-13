# tensor_qtl_cis

Nominal cis tensorQTL mapping for Terra managed Cromwell, using WDL 1.0.
The source is the supplied `tensorqtl_cis_nominal_v1-0_BETA_modified.1.wdl`.

## Runtime

The default runtime uses `g2-standard-16`, NVIDIA L4, `us-central1-c`,
and the existing `gcr.io/broad-cga-francois-gtex/tensorqtl:latest` image.
Use `num_gpus=1`, `num_threads=16`, and `memory=64` (GB) for this machine.
The GPU count defaults to one. Disk size and preemptible retries remain inputs.
The existing image uses a mutable tag; this change does not build a new image.

The machine type, GPU type, and zones are configurable:

| Input | Type | Default |
|---|---|---|
| `machine_type` | `String?` | `g2-standard-16` |
| `gpu_type` | `String` | `nvidia-l4` |
| `zones` | `Array[String]` | `["us-central1-c"]` |

A predefined machine type selects its own CPU, memory, and GPUs. Cromwell
ignores `num_threads`, `memory`, `gpu_type`, and `num_gpus` for that selection.
For example, select `a2-highgpu-1g` to request a machine with one A100 GPU.
Select zones that support that machine and have sufficient quota and capacity.

To request an attached GPU such as a T4, set `machine_type` to `null`.
Cromwell then selects a machine from `num_threads`, `memory`, `gpu_type`,
and `num_gpus`. Add these entries to the existing Terra workflow inputs:

```json
{
  "tensorqtl_cis_nominal_workflow.tensorqtl_cis_nominal.machine_type": null,
  "tensorqtl_cis_nominal_workflow.tensorqtl_cis_nominal.gpu_type": "nvidia-tesla-t4",
  "tensorqtl_cis_nominal_workflow.tensorqtl_cis_nominal.num_gpus": 1,
  "tensorqtl_cis_nominal_workflow.tensorqtl_cis_nominal.zones": ["us-central1-c"]
}
```

Changing the GPU type can avoid the L4 quota limit. The selected GPU must
have sufficient regional quota. The global GPU quota can also apply.
An available quota does not guarantee GPU capacity.
See [Cromwell runtime attributes](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#predefinedmachinetype-alpha)
and [Google Cloud GPU quotas](https://docs.cloud.google.com/compute/resource-usage#gpu_quota).

## Inputs and outputs

Input keys retain the prefix
`tensorqtl_cis_nominal_workflow.tensorqtl_cis_nominal.`:

- Required files: `plink_pgen`, `plink_pvar`, `plink_psam`, `phenotype_bed`, `covariates`.
- Optional files: `interaction`, `phenotype_groups`.
- Required settings: `prefix`, `memory`, `disk_space`, `num_threads`, `num_preempt`.
- Optional settings: `num_gpus` (default 1), `machine_type`, `gpu_type`, `zones`.

All input files remain typed Files until command rendering. The task checks
localized paths and links the PLINK files to a common basename, even when Terra
places them in different directories. Optional flags are passed only when their
files are supplied. The prefix must start with a letter or digit and contain
only letters, digits, dots, underscores, or hyphens.

The task and workflow expose `chr_parquet` and `log`. The supplied workflow's
nominal cis mode and analysis options are preserved.

## Validation

GitHub Actions checks WDL syntax, runtime defaults and overrides, the rendered command, optional
files, safe path quoting, and rejection of unresolved cloud URIs. An AST check
rejects workflow-scope file-writing functions. A CPU smoke test runs the actual
tensorQTL engine on generated PLINK2 data and checks a known association.
No local Docker build is required.

CPU and command tests do not validate GPU execution, the production image, or
Terra localization. The complete workflow has not been tested on Terra.
