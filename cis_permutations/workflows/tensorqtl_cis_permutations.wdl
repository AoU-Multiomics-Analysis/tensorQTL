version 1.0

task tensorqtl_cis_permutations {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        String prefix
        File? phenotype_groups
        Float? fdr
        Float? qvalue_lambda
        Float? pval_threshold
        Int? seed
        String? flags
        Int cis_window = 1000000
        Int memory
        Int disk_space
        Int num_threads
        Int num_gpus
        Int num_preempt
        String? machine_type = "g2-standard-16"
        String gpu_type = "nvidia-l4"
        Array[String] zones = ["us-central1-c"]
    }

    command <<<
        set -euo pipefail
        log() { echo "[$(date -u +%FT%TZ)] stage=tensorqtl_cis_permutations $*" >&2; }
        log "Validate localized inputs."
        plink_pgen='~{sub(plink_pgen, "'", "'\"'\"'")}'
        plink_pvar='~{sub(plink_pvar, "'", "'\"'\"'")}'
        plink_psam='~{sub(plink_psam, "'", "'\"'\"'")}'
        phenotype_bed='~{sub(phenotype_bed, "'", "'\"'\"'")}'
        covariates='~{sub(covariates, "'", "'\"'\"'")}'
        phenotype_groups='~{sub(select_first([phenotype_groups, ""]), "'", "'\"'\"'")}'
        prefix='~{sub(prefix, "'", "'\"'\"'")}'
        for path in "$plink_pgen" "$plink_pvar" "$plink_psam" "$phenotype_bed" "$covariates"; do
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done
        if [[ -n "$phenotype_groups" ]]; then
            path="$phenotype_groups"
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        fi
        [[ "$prefix" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { log "Input error: prefix must be a filename prefix."; exit 1; }
        # Cromwell can localize the three PLINK files into different directories.
        mkdir -p localized_genotypes
        ln -s "$plink_pgen" localized_genotypes/input.pgen
        ln -s "$plink_pvar" localized_genotypes/input.pvar
        ln -s "$plink_psam" localized_genotypes/input.psam
        args=(localized_genotypes/input "$phenotype_bed" "$prefix" --mode cis
              --covariates "$covariates" --window ~{cis_window})
        if [[ -n "$phenotype_groups" ]]; then args+=(--phenotype_groups "$phenotype_groups"); fi
        fdr='~{default="" fdr}'
        qvalue_lambda='~{default="" qvalue_lambda}'
        pval_threshold='~{default="" pval_threshold}'
        seed='~{default="" seed}'
        if [[ -n "$fdr" ]]; then args+=(--fdr "$fdr"); fi
        if [[ -n "$qvalue_lambda" ]]; then args+=(--qvalue_lambda "$qvalue_lambda"); fi
        if [[ -n "$pval_threshold" ]]; then args+=(--pval_threshold "$pval_threshold"); fi
        if [[ -n "$seed" ]]; then args+=(--seed "$seed"); fi
        flags='~{sub(select_first([flags, ""]), "'", "'\"'\"'")}'
        if [[ -n "$flags" ]]; then
            # Split shell-style arguments without executing shell expressions.
            python3 -c 'import shlex, sys; sys.stdout.buffer.write(b"".join(arg.encode() + b"\0" for arg in shlex.split(sys.argv[1])))' "$flags" > extra_flags.list
            while IFS= read -r -d '' flag; do args+=("$flag"); done < extra_flags.list
        fi
        log "Start cis permutation mapping."
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
        File cis_qtl = prefix + ".cis_qtl.txt.gz"
        File log = prefix + ".tensorQTL.cis.log"
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

workflow tensorqtl_cis_permutations_workflow {
    call tensorqtl_cis_permutations
}
