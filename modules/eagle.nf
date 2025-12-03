// process phase_with_eagle {
//     tag "Eagle Phasing"

//     input:
//         tuple path(input_genotype), path(input_genotype_index), path(genetic_map), val(chromosome)

//     output:
//         path "*.vcf.gz", emit: phased_vcf

//     script:
//         """
//         eagle \
//             --vcfTarget=${input_genotype} \
//             --geneticMapFile=${genetic_map} \
//             --chrom=${chromosome} \
//             --outPrefix=test \
//             2>&1 | tee eagle.log
//         """
// }
process phase_with_eagle {
    tag "Eagle Phasing"

    input:
        tuple path(input_genotype), path(input_genotype_index), path(reference_vcf), path(ref_vcf_index), path(genetic_map), val(chromosome)

    output:
        path "*.vcf.gz", emit: phased_vcf

    script:
        """
        eagle \
            --vcfTarget=${input_genotype} \
            --vcfRef=${reference_vcf} \
            --geneticMapFile=${genetic_map} \
            --chrom=${chromosome} \
            --outPrefix=test \
            2>&1 | tee eagle.log
        """
}