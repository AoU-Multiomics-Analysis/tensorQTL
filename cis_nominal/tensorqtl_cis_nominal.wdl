version 1.0

task tensorqtl_cis_nominal {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        String prefix
        File? interaction
        File? phenotype_groups
        Int memory
        Int disk_space
        Int num_threads
        Int num_gpus = 1
        Int num_preempt
        String? machine_type = "g2-standard-16"
        String gpu_type = "nvidia-l4"
        Array[String] zones = ["us-central1-c"]
    }

    command <<<
        set -euo pipefail
        log() { echo "[$(date -u +%FT%TZ)] stage=tensorqtl_cis_nominal $*" >&2; }
        log "Validate localized inputs."
        plink_pgen='~{sub(plink_pgen, "'", "'\"'\"'")}'
        plink_pvar='~{sub(plink_pvar, "'", "'\"'\"'")}'
        plink_psam='~{sub(plink_psam, "'", "'\"'\"'")}'
        phenotype_bed='~{sub(phenotype_bed, "'", "'\"'\"'")}'
        covariates='~{sub(covariates, "'", "'\"'\"'")}'
        interaction='~{sub(select_first([interaction, ""]), "'", "'\"'\"'")}'
        phenotype_groups='~{sub(select_first([phenotype_groups, ""]), "'", "'\"'\"'")}'
        prefix='~{sub(prefix, "'", "'\"'\"'")}'
        for path in "$plink_pgen" "$plink_pvar" "$plink_psam" "$phenotype_bed" "$covariates"; do
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done
        for path in "$interaction" "$phenotype_groups"; do
            [[ -n "$path" ]] || continue
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done
        [[ "$prefix" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { log "Input error: prefix must be a filename prefix."; exit 1; }
        # Cromwell can localize the three PLINK files into different directories.
        mkdir -p localized_genotypes
        ln -s "$plink_pgen" localized_genotypes/input.pgen
        ln -s "$plink_pvar" localized_genotypes/input.pvar
        ln -s "$plink_psam" localized_genotypes/input.psam
        args=(localized_genotypes/input "$phenotype_bed" "$prefix" --mode cis_nominal
              --covariates "$covariates")
        if [[ -n "$interaction" ]]; then args+=(--interaction "$interaction"); fi
        if [[ -n "$phenotype_groups" ]]; then args+=(--phenotype_groups "$phenotype_groups"); fi
        log "Start nominal cis mapping; interaction enabled: ~{defined(interaction)}."
        python3 -m tensorqtl "${args[@]}"
        log "Mapping completed."
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        bootDiskSizeGb: 25
        cpu: num_threads
        preemptible: num_preempt
        # Batch requires String. An empty string lets Batch select the machine.
        predefinedMachineType: select_first([machine_type, ""])
        gpuType: gpu_type
        gpuCount: num_gpus
        zones: zones
    }


    output {
        Array[File] chr_parquet = glob(prefix + "*.parquet")
        File log = glob(prefix + "*.log")[0]
    }
    parameter_meta {
        machine_type: "GCP machine type. Default: g2-standard-16 with L4. Set null or an empty string to let Batch select the machine from cpu, memory, gpu_type, and num_gpus."
        gpu_type: "GCP GPU type. Default: nvidia-l4. The GPU type and count must be compatible with the selected machine."
        zones: "GCP zones that support the selected machine and GPU. Default: us-central1-c."
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow tensorqtl_cis_nominal_workflow {
    call tensorqtl_cis_nominal
    output {
        Array[File] chr_parquet = tensorqtl_cis_nominal.chr_parquet
        File log = tensorqtl_cis_nominal.log
    }
}
