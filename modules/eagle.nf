process phase_with_eagle {
    tag "Eagle Phasing"

    input:
        path input_genotype
        path reference_vcf
        path genetic_map
        val chromosome

    output:
        path "*.vcf.gz", emit: phased_vcf

    script:
        """
        eagle \
            --vcf=${input_genotype} \
            --ref=${reference_vcf} \
            --geneticMapFile=${genetic_map} \
            --chrom=${chromosome}
        """
}
z