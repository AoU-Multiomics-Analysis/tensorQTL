# Phenotype-count sharding implementation plan

The user requested removal of the variant-count limit. The selected design
uses `max_phenotypes_per_shard` as the only shard-size target. Keep complete
groups together. Preserve chromosome ordering, complete cis windows, variant
counts in audit output, and the default phenotype target of 250. A target of
5000 is an input choice, not a new default. This narrow change needs no new
analysis method or cloud run.

- [x] Write failing tests for count-based packing, a 5000-phenotype target
  with more than 100000 cis variants, final partial shards, and intact groups.
  Check that WDL and CLI interfaces no longer expose the variant limit.
- [x] Remove the variant limit from the planner, CLI, WDL, and manifest's
  oversized condition. Count variants only for extraction validation and audit.
- [x] Update real PLINK2/tensorQTL smoke tests, input examples, and user docs.
  Document removal of the old input and the need for the updated image.
- [x] Run local tests, WDL validation, and an independent code review.

Delivery: push a review branch and monitor the existing GitHub Actions
checks. Record remote results in the pull request. Build Docker only in
Actions. Do not merge or submit Terra jobs without authorization.

## Design choices

Fixed row chunks would split groups. The selected approach keeps group-aware
packing and removes the variant-count condition. The final shard can be
smaller than the target. A group larger than the target is an oversized
shard. The earlier preparation memory issue is outside this change; it
still reads the full phenotype matrix before creating shards.

## Interfaces

`plan_shards(rows, index, window, max_phenotypes)` returns coordinate-sorted
phenotype shards. `index` and `window` remain necessary for group validation.
The WDL input `max_phenotypes_per_shard` and CLI argument
`--max-phenotypes-per-shard` retain their names and defaults. Remove
`max_variants_per_shard` and `--max-variants-per-shard` from active interfaces.

## Verification before push

- New CLI, WDL, and planner checks failed before implementation.
- All 82 local tests passed: 43 permutation, 11 nominal cis, and 28 trans.
- miniwdl accepts the sharded WDL, with the existing Terra runtime and
  optional-default warnings. actionlint and `git diff --check` pass.
- The example Terra JSON passes WDL input type checks with a target of 5000.
- Independent code review found no actionable defects.
- This change does not establish a successful complete Terra run.
