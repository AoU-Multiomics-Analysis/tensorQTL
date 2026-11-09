# Merge file list localization fix

The Terra merge command tries to open a generated `write_lines` file at a
GCS URI. `sub(write_lines(results), ...)` converts the File to String before
Cromwell applies its command File mapper. The input result files can be
localized correctly while the generated list path remains a cloud URI.

Create the newline list inside the task command with a quoted heredoc.
Interpolate `results` directly as an Array[File] with a newline separator.
This keeps File types until command rendering and preserves shell special
characters as literal text. Reject an empty array before computation.
Keep the merge CLI, image, and calculation unchanged.

- [x] Reproduce the reported failure with a test standard library that
  assigns cloud paths to generated command files. Test one and several
  results, path order, shell characters, and an empty array.
- [x] Replace the generated list expression. Keep input readability and
  unresolved-URI checks. Update the documentation.
- [x] Run the test suites, WDL checks, and an independent review.

Delivery: push a fix PR and monitor the existing GitHub CPU smoke test.

Review focus: direct File-array interpolation, safe heredoc quoting,
empty input rejection, no workflow-scope writes, and existing-image
compatibility. Do not build Docker locally or submit Terra jobs. Local
and GitHub checks do not establish a complete successful Terra run.

Evidence before delivery: the cloud-name regression failed against the old
WDL with the reported missing local file error. It passes after the fix for
one and three result files. Empty results fail before the script starts.
All 82 local tests pass: 43 permutation, 11 nominal cis, and 28 trans. The
trans suite emits the existing sandbox CPU-query warnings. miniwdl passes
with three existing runtime/optional-input warnings. actionlint and the
diff whitespace check pass. Preparation, mapping, and image scripts are
unchanged, so their successful calls remain eligible for cache reuse.
Independent review found no defects. It repeated the five WDL command
tests and checked Cromwell 92's indentation, newline separator, and
per-element File mapping before interpolation.
