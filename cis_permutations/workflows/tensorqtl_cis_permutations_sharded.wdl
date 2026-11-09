version 1.0

workflow tensorqtl_cis_permutations_sharded {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        File? phenotype_groups
        String prefix
        String docker_image
        Int cis_window = 1000000
        Int max_phenotypes_per_shard = 250
        Int max_variants_per_shard = 100000
        Int min_mac = 10
        Float max_missing_fraction = 0.05
        Int permutations = 10000
        Float maf_threshold = 0.0
        Int seed = 1
        Boolean dosages = false
        Float fdr = 0.05
        Float? qvalue_lambda
        Int prepare_cpu = 4
        Int prepare_memory_gb = 32
        Int prepare_disk_gb = 500
        Int map_cpu = 16
        Int map_memory_gb = 64
        Int map_disk_gb = 100
        String? machine_type = "g2-standard-16"
        String gpu_type = "nvidia-l4"
        Int gpu_count = 1
        Array[String] zones = ["us-central1-c"]
        Int merge_memory_gb = 8
        Int merge_disk_gb = 20
        Int preemptible = 2
    }

    call prepare_cis {
        input:
            plink_pgen = plink_pgen, plink_pvar = plink_pvar, plink_psam = plink_psam,
            phenotype_bed = phenotype_bed, covariates = covariates,
            phenotype_groups = phenotype_groups, docker_image = docker_image,
            cis_window = cis_window, max_phenotypes_per_shard = max_phenotypes_per_shard,
            max_variants_per_shard = max_variants_per_shard,
            min_mac = min_mac, max_missing_fraction = max_missing_fraction, dosages = dosages,
            cpu = prepare_cpu, memory_gb = prepare_memory_gb, disk_gb = prepare_disk_gb,
            zones = zones, preemptible = preemptible
    }
    # The five glob outputs use the same fixed-width shard names and Bash order.
    # Consume the returned Files: Cromwell can put each glob in a separate directory.
    scatter (shard_index in range(length(prepare_cis.plink_pgens))) {
        call map_cis {
            input:
                plink_pgen = prepare_cis.plink_pgens[shard_index],
                plink_pvar = prepare_cis.plink_pvars[shard_index],
                plink_psam = prepare_cis.plink_psams[shard_index],
                phenotype_bed = prepare_cis.phenotype_beds[shard_index],
                phenotype_groups = prepare_cis.phenotype_groups_files[shard_index],
                covariates = covariates,
                prefix = basename(prepare_cis.plink_pgens[shard_index], ".pgen"),
                docker_image = docker_image,
                cis_window = cis_window, permutations = permutations,
                maf_threshold = maf_threshold, seed = seed, dosages = dosages,
                cpu = map_cpu, memory_gb = map_memory_gb, disk_gb = map_disk_gb,
                machine_type = machine_type, gpu_type = gpu_type, gpu_count = gpu_count,
                zones = zones, preemptible = preemptible
        }
    }
    call merge_cis {
        input:
            results = map_cis.cis_qtl, prefix = prefix, docker_image = docker_image,
            fdr = fdr, qvalue_lambda = qvalue_lambda,
            memory_gb = merge_memory_gb, disk_gb = merge_disk_gb,
            zones = zones, preemptible = preemptible
    }

    output {
        File cis_qtl = merge_cis.cis_qtl
        File manifest = prepare_cis.manifest
        File variant_qc = prepare_cis.variant_qc
        File preparation_skipped = prepare_cis.skipped
        File preparation_log = prepare_cis.log
        File merge_log = merge_cis.log
        Array[File] shard_results = map_cis.cis_qtl
        Array[File] shard_logs = map_cis.log
        Array[File] mapping_skipped = map_cis.skipped
    }
    parameter_meta {
        docker_image: "Published sharded-cis image. Use an immutable image digest or commit tag."
        cis_window: "Distance in bp on each side of a phenotype's one-based start and end."
        max_variants_per_shard: "Target variant count in the union of complete cis windows. A single phenotype or group can exceed it."
        max_phenotypes_per_shard: "Target phenotype count. Groups are kept complete."
        min_mac: "Minimum observed diploid-coded minor allele count before sharding and imputation, across all phenotype samples."
        max_missing_fraction: "Maximum fraction of phenotype samples with missing genotypes before sharding and imputation."
        prepare_disk_gb: "Space for the source PLINK triplet, all extracted shards, and temporary files."
        dosages: "Use PGEN dosages for preparation QC and mapping instead of hard calls. PLINK extraction retains dosage data."
    }
    meta {
        description: "Prepare complete cis windows, scatter phenotype shards, then calculate global Storey q-values."
    }
}

task prepare_cis {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        File? phenotype_groups
        String docker_image
        Int cis_window = 1000000
        Int max_phenotypes_per_shard = 250
        Int max_variants_per_shard = 100000
        Int min_mac = 10
        Float max_missing_fraction = 0.05
        Boolean dosages = false
        Int cpu = 4
        Int memory_gb = 32
        Int disk_gb = 500
        Array[String] zones = ["us-central1-c"]
        Int preemptible = 2
    }
    command <<<
        set -euo pipefail
        {
        log() { echo "[$(date -u +%FT%TZ)] stage=prepare $*" >&2; }
        check_file() {
            case "$1" in *://*) log "Input localization error: unresolved URI $1"; exit 1;; esac
            [[ -r "$1" ]] || { log "Input localization error: unreadable file $1"; exit 1; }
        }
        pgen='~{sub(plink_pgen, "'", "'\"'\"'")}'
        pvar='~{sub(plink_pvar, "'", "'\"'\"'")}'
        psam='~{sub(plink_psam, "'", "'\"'\"'")}'
        bed='~{sub(phenotype_bed, "'", "'\"'\"'")}'
        cov='~{sub(covariates, "'", "'\"'\"'")}'
        groups='~{sub(select_first([phenotype_groups, ""]), "'", "'\"'\"'")}'
        for path in "$pgen" "$pvar" "$psam" "$bed" "$cov"; do check_file "$path"; done
        args=(prepare --plink-pgen "$pgen" --plink-pvar "$pvar" --plink-psam "$psam"
              --phenotype-bed "$bed" --covariates "$cov" --cis-window ~{cis_window}
              --max-phenotypes-per-shard ~{max_phenotypes_per_shard}
              --max-variants-per-shard ~{max_variants_per_shard}
              --min-mac ~{min_mac} --max-missing-fraction ~{max_missing_fraction}
              --threads ~{cpu} --plink-memory-mb ~{memory_gb * 512})
        if ~{dosages}; then args+=(--dosages); fi
        if [[ -n "$groups" ]]; then check_file "$groups"; args+=(--phenotype-groups "$groups"); fi
        log "Start phenotype and genotype preparation."
        python3 -u /opt/tensorqtl/sharded_cis.py "${args[@]}"
        log "Preparation completed."
        } 2>&1 | tee preparation.log
    >>>
    output {
        # Google Batch must discover the upload patterns before the command runs.
        # Files named inside task-created JSON are not discoverable at that point.
        Array[File] plink_pgens = glob("shards/*.pgen")
        Array[File] plink_pvars = glob("shards/*.pvar")
        Array[File] plink_psams = glob("shards/*.psam")
        Array[File] phenotype_beds = glob("shards/*.phenotypes.bed.gz")
        Array[File] phenotype_groups_files = glob("shards/*.groups.tsv")
        File manifest = "shards.manifest.tsv"
        File variant_qc = "variant_qc.tsv"
        File skipped = "preparation.skipped.tsv"
        File log = "preparation.log"
    }
    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{disk_gb} SSD"
        bootDiskSizeGb: 30
        zones: zones
        preemptible: preemptible
    }
}

task map_cis {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        File phenotype_groups
        String prefix
        String docker_image
        Int cis_window = 1000000
        Int permutations = 10000
        Float maf_threshold = 0.0
        Int seed = 1
        Boolean dosages = false
        Int cpu = 16
        Int memory_gb = 64
        Int disk_gb = 100
        String? machine_type = "g2-standard-16"
        String gpu_type = "nvidia-l4"
        Int gpu_count = 1
        Array[String] zones = ["us-central1-c"]
        Int preemptible = 2
    }
    command <<<
        set -euo pipefail
        {
        log() { echo "[$(date -u +%FT%TZ)] stage=map $*" >&2; }
        check_file() {
            case "$1" in *://*) log "Input localization error: unresolved URI $1"; exit 1;; esac
            [[ -r "$1" ]] || { log "Input localization error: unreadable file $1"; exit 1; }
        }
        pgen='~{sub(plink_pgen, "'", "'\"'\"'")}'
        pvar='~{sub(plink_pvar, "'", "'\"'\"'")}'
        psam='~{sub(plink_psam, "'", "'\"'\"'")}'
        bed='~{sub(phenotype_bed, "'", "'\"'\"'")}'
        cov='~{sub(covariates, "'", "'\"'\"'")}'
        groups='~{sub(phenotype_groups, "'", "'\"'\"'")}'
        prefix='~{sub(prefix, "'", "'\"'\"'")}'
        for path in "$pgen" "$pvar" "$psam" "$bed" "$cov" "$groups"; do check_file "$path"; done
        export OMP_NUM_THREADS=~{cpu} MKL_NUM_THREADS=~{cpu} OPENBLAS_NUM_THREADS=~{cpu}
        args=(map --plink-pgen "$pgen" --plink-pvar "$pvar" --plink-psam "$psam"
              --phenotype-bed "$bed" --covariates "$cov" --phenotype-groups "$groups"
              --prefix "$prefix" --cis-window ~{cis_window} --permutations ~{permutations}
              --maf-threshold ~{maf_threshold} --seed ~{seed} --threads ~{cpu})
        if ~{dosages}; then args+=(--dosages); fi
        if [[ ~{gpu_count} -gt 0 ]]; then args+=(--require-gpu); fi
        log "Start cis permutations; window=~{cis_window}; permutations=~{permutations}."
        python3 -u /opt/tensorqtl/sharded_cis.py "${args[@]}"
        log "Mapping completed."
        } 2>&1 | tee mapping.log
    >>>
    output {
        File cis_qtl = prefix + ".cis_qtl.txt.gz"
        File skipped = prefix + ".skipped.tsv"
        File log = "mapping.log"
    }
    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{disk_gb} SSD"
        bootDiskSizeGb: 30
        predefinedMachineType: select_first([machine_type, ""])
        gpuType: gpu_type
        gpuCount: gpu_count
        zones: zones
        preemptible: preemptible
    }
}

task merge_cis {
    input {
        Array[File] results
        String prefix
        String docker_image
        Float fdr = 0.05
        Float? qvalue_lambda
        Int memory_gb = 8
        Int disk_gb = 20
        Array[String] zones = ["us-central1-c"]
        Int preemptible = 2
    }
    command <<<
        set -euo pipefail
        {
        log() { echo "[$(date -u +%FT%TZ)] stage=merge $*" >&2; }
        if [[ ~{length(results)} -eq 0 ]]; then
            log "Input error: no shard result files were supplied."
            exit 1
        fi
        # Keep the File array typed until interpolation. A quoted heredoc
        # writes the localized paths without expanding shell characters.
        log "Create local list for ~{length(results)} shard result files."
        results_list='shard_results.list'
        cat > "$results_list" <<'TENSORQTL_LOCALIZED_RESULTS'
        ~{sep="\n" results}
        TENSORQTL_LOCALIZED_RESULTS
        while IFS= read -r path; do
            case "$path" in *://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done < "$results_list"
        prefix='~{sub(prefix, "'", "'\"'\"'")}'
        args=(merge --results-list "$results_list" --prefix "$prefix" --fdr ~{fdr})
        lambda='~{default="" qvalue_lambda}'
        if [[ -n "$lambda" ]]; then args+=(--qvalue-lambda "$lambda"); fi
        log "Calculate global Storey q-values from all permutation summaries."
        python3 -u /opt/tensorqtl/sharded_cis.py "${args[@]}"
        log "Merge completed."
        } 2>&1 | tee merge.log
    >>>
    output {
        File cis_qtl = prefix + ".cis_qtl.txt.gz"
        File log = "merge.log"
    }
    runtime {
        docker: docker_image
        cpu: 2
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{disk_gb} SSD"
        bootDiskSizeGb: 30
        zones: zones
        preemptible: preemptible
    }
}
