# Variant QC before cis sharding

The user requested MAC and variant missingness filters before sharding.
The selected defaults are `min_mac=10` and `max_missing_fraction=0.05`.

Use the existing CPU preparation task. Read the phenotype sample subset
from PGEN in bounded batches, with the same hard-call or dosage mode used
by mapping. Calculate MAC from observed diploid-coded values, before mean
imputation. Retain equality at both cutoffs. Count all analysis samples,
including nonfounders. These counts do not use the imputed output MAC.

Write a retained-variant list and a filtered position index. Plan shards
from that index. Intersect each complete cis window with the retained list
during PLINK extraction. Do not create another full genotype dataset.
Report original, retained, low-MAC, high-missingness, and overlapping failure
counts in a task File output. Record phenotypes whose windows become empty.
Fail before scatter if no valid phenotypes remain.

## Approaches

- Direct PGEN QC matches mapping mode and avoids another genotype copy.
  This is the selected implementation approach.
- A separate filtered PGEN is simpler but requires another large genotype
  copy on the preparation disk.
- Plain PLINK MAC filtering can count dosages while mapping uses hard calls.
  Do not introduce this mismatch.

## Constraints and validation

Keep WDL 1.0 and Terra File localization. Use named CLI arguments and stage
logs. Keep existing mapping-time MAF filtering. Preserve the script-only
image rebuild trigger. Add boundary, sample-subset, dosage, empty-window,
and cloud-localization tests. Extend the existing GitHub Actions smoke test.
Do not build a local Docker image or submit a Terra run.
