process run_rfmix {
    tag "RFMix Ancestry"

    container 'your-dockerhub-username/eagle-rfmix:latest'

    input:
    path phased_vcf
    path ref_vcf
    path sample_map
    path genetic_map
    val chromosome
    val output_prefix

    output:
    path "${output_prefix}_rfmix/"

    script:
    """
    mkdir -p ${output_prefix}_rfmix
    rfmix \
        -f ${phased_vcf} \
        -r ${ref_vcf} \
        -m ${sample_map} \
        -g ${genetic_map} \
        --chromosome=${chromosome} \
        -o ${output_prefix}_rfmix/results
    """
}
