process phase_with_eagle {
    tag "Eagle Phasing ${chromosome}"

    input:
        tuple path(input_genotype),
              path(input_genotype_index),
              path(reference_vcf),
              path(ref_vcf_index),
              path(genetic_map),
              val(chromosome)

    output:
        tuple val(chromosome),
              path("phased_${chromosome}.vcf.gz"),
              path("phased_${chromosome}.vcf.gz.tbi"),
              emit: phased_vcf

    script:
        """
        set -euo pipefail

        echo "Running Eagle for ${chromosome}"

        eagle \
            --vcfTarget=${input_genotype} \
            --vcfRef=${reference_vcf} \
            --geneticMapFile=${genetic_map} \
            --chrom=${chromosome} \
            --outPrefix=phased_${chromosome} \
            > eagle_${chromosome}.log 2>&1

        # Hard failure if Eagle didn't produce output
        if [ ! -f phased_${chromosome}.vcf.gz ]; then
            echo "ERROR: Eagle did not produce phased_${chromosome}.vcf.gz"
            cat eagle_${chromosome}.log
            exit 1
        fi

        # Index output
        tabix -f -p vcf phased_${chromosome}.vcf.gz

        echo "Eagle completed successfully for ${chromosome}"
        """
}