process phase_with_eagle {
    tag "Eagle Phasing"

    container 'your-dockerhub-username/eagle-rfmix:latest'

    input:
    path input_vcf
    path ref_vcf
    path genetic_map
    val chromosome
    val output_prefix

    output:
    path "${output_prefix}_phased.vcf.gz"

    script:
    """
    eagle \
        --vcf=${input_vcf} \
        --ref=${ref_vcf} \
        --geneticMapFile=${genetic_map} \
        --chrom=${chromosome} \
        --outPrefix=${output_prefix}_phased
    """
}
