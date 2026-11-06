# Shard output collection fix

The log confirms that preparation wrote all 163 shards. Cromwell then
failed to download a shard from GCS. Its Batch backend builds an output
transfer list before execution. At that stage, `read_json("shards.json")`
exposes only the JSON file. The File fields inside it do not cause uploads.

Use five literal `glob` outputs with type `Array[File]`. Scatter over their
indices and pass the File values directly. Fixed-width shard names give
the arrays the same Bash sort order. Derive the prefix from the PGEN
basename. Keep the Python script and published image unchanged. This fix
is separate from the phenotype-count change in PR #8.

- [x] Add failing tests for static upload discovery and cloud path changes.
  Include three shards, shard 158, and empty groups files.
- [x] Replace the output and scatter declarations. Update type checks and
  documentation. Preserve WDL 1.0 and task-only file writes.
- [x] Run tests and WDL checks. Obtain independent review.

Delivery: push a fix PR and monitor GitHub Actions. Do not build Docker
locally or submit a Terra job.

Review focus: no File paths from dynamic JSON may enter the scatter; glob
outputs stay typed File; all five families use the same shard index; empty
groups files remain outputs; the existing image remains compatible.

Local evidence: both new tests failed on the old output declaration and
passed after the fix. All 80 tests passed: 41 permutation, 11 nominal cis,
and 28 trans. miniwdl accepts the WDL with the existing runtime and optional
default warnings. Actual Terra bucket contents could not be checked because
the available credentials require reauthentication. The uploaded preparation
log and Cromwell 92 source establish the file discovery failure mechanism.
Independent review found no code defect. It repeated all 41 permutation
tests, miniwdl validation, and the diff whitespace check successfully.
