# tensorQTL

Terra WDL pipelines for nominal cis, cis permutation, and trans QTL mapping.

<!-- workflow-badges:start -->
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/tensorQTL/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/tensorQTL/actions/workflows/update-readme-badges.yml)
[![Validate nominal cis WDL](https://github.com/AoU-Multiomics-Analysis/tensorQTL/actions/workflows/validate-cis-nominal.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/tensorQTL/actions/workflows/validate-cis-nominal.yml)
[![Validate cis permutations WDL](https://github.com/AoU-Multiomics-Analysis/tensorQTL/actions/workflows/validate-cis-permutations.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/tensorQTL/actions/workflows/validate-cis-permutations.yml)
[![Validate WDL and chromosome mapping](https://github.com/AoU-Multiomics-Analysis/tensorQTL/actions/workflows/validate-trans.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/tensorQTL/actions/workflows/validate-trans.yml)
<!-- workflow-badges:end -->

## Pipelines

| Pipeline | WDL file | Documentation |
|---|---|---|
| Nominal cis | [cis_nominal/tensorqtl_cis_nominal.wdl](cis_nominal/tensorqtl_cis_nominal.wdl) | [Inputs and runtime](cis_nominal/README.md) |
| Cis eQTL permutations | [cis_permutations/workflows/tensorqtl_cis_permutations.wdl](cis_permutations/workflows/tensorqtl_cis_permutations.wdl) | [Inputs and runtime](cis_permutations/README.md) |
| Cis sQTL permutations | [cis_permutations/workflows/tensorqtl_cis_permutations_sQTL.wdl](cis_permutations/workflows/tensorqtl_cis_permutations_sQTL.wdl) | [Inputs and runtime](cis_permutations/README.md) |
| Trans | [trans/tensorQTL_trans.wdl](trans/tensorQTL_trans.wdl) | [Inputs and runtime](trans/README.md) |

All four entry points use WDL 1.0. The root `.dockstore.yml` defines their descriptor paths and names.
Register this repository in Dockstore to import these workflows into Terra.
Keep the sQTL descriptor and its imported eQTL descriptor in the same directory.

The sQTL entry point uses the existing cis permutation task. It requires a groups file at
`tensorqtl_cis_permutations_workflow.phenotype_groups`.
Other task inputs keep the prefix `tensorqtl_cis_permutations_workflow.tensorqtl_cis_permutations.`.
The sQTL task now uses the eQTL task's L4 runtime defaults, input checks, safe quoting, and command logs.
Its optional machine, GPU, and zone inputs can override those defaults.

## Checks

The root GitHub Actions workflows run each pipeline's checks from its own directory.
They check WDL syntax, workflow-scope file writes, and simulated cloud-to-local file handling.
The cis permutation and trans checks also use the Cromwell WDL parser.
The nominal cis and trans checks run the existing tensorQTL CPU smoke tests in GitHub Actions.
No local Docker build is needed. The pipelines use their existing container references.
The source Docker build files are empty placeholders; their nested Actions files are inactive.

For local command tests, use Python 3.11 with `miniwdl==1.14.2`, `numpy==1.26.4`,
`pandas==2.2.3`, and `pyarrow==19.0.1`. Run each suite in a separate process:

```bash
(cd cis_nominal && python -m unittest discover -s tests)
(cd cis_permutations && python -m unittest discover -s tests)
(cd trans && python -m unittest discover -s tests)
```

The complete workflows have not been tested on Terra during this migration.
The command tests and CPU smoke tests do not confirm GPU execution or Terra localization.

## Source repositories

The import commits preserve the default-branch Git history from these repositories:

- [tensor_qtl_cis](https://github.com/AoU-Multiomics-Analysis/tensor_qtl_cis) → `cis_nominal/`
- [tensorQTL_cis_permutations](https://github.com/AoU-Multiomics-Analysis/tensorQTL_cis_permutations) → `cis_permutations/`
- [tensorQTL_trans](https://github.com/AoU-Multiomics-Analysis/tensorQTL_trans) → `trans/`

The source repositories remain available. Source tags, other branches, issues, releases,
and repository settings are not part of this code import.
Each folder contains its source documentation and configuration files.
GitHub Actions and Dockstore use the configuration files at the root of this repository.
