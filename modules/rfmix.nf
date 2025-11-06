process run_rfmix {
    tag "RFMix Ancestry"

    input:
        tuple val(chromosome), file(phased_vcf), file(reference_vcf), file(reference_vcf_index), file(genetic_map), file(sample_map)

    output:
        path "rfmix_${chromosome}_results/*"

    script:
        """
        mkdir -p rfmix_${chromosome}_results
        rfmix \
            -f ${phased_vcf} \
            -r ${reference_vcf} \
            -m ${sample_map} \
            -g ${genetic_map} \
            -o ${chromosome} \
            --chromosome=${chromosome}
        """
}
