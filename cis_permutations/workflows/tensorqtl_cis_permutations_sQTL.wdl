version 1.0

import "tensorqtl_cis_permutations.wdl" as cis

workflow tensorqtl_cis_permutations_workflow {
    input {
        File phenotype_groups
    }

    call cis.tensorqtl_cis_permutations {
        input:
            phenotype_groups = phenotype_groups
    }

    output {
        File cis_qtl = tensorqtl_cis_permutations.cis_qtl
        File log = tensorqtl_cis_permutations.log
    }
}
