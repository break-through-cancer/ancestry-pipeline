process run_rfmix {
    tag "RFMix Ancestry"

    input:
        path phased_vcf
        path ref_vcf
        path sample_map
        path genetic_map
        val chromosome

    output:
        path "rfmix_${chromosome}_results/*"

    script:
        """
        mkdir -p rfmix_${chromosome}_results
        rfmix \
            -f ${phased_vcf} \
            -r ${ref_vcf} \
            -m ${sample_map} \
            -g ${genetic_map} \
            --chromosome=${chromosome}
        """
}
