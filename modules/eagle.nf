process phase_with_eagle {
    tag "Eagle Phasing"

    input:
        path input_vcf
        path ref_vcf
        path genetic_map
        val chromosome

    output:
        path "${output_prefix}_phased.vcf.gz"

    script:
        """
        eagle \
            --vcf=${input_vcf} \
            --ref=${ref_vcf} \
            --geneticMapFile=${genetic_map} \
            --chrom=${chromosome}
        """
}
