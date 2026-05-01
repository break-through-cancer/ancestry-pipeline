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

        echo "=== Eagle input check ==="
        echo "chromosome=${chromosome}"
        echo "Target VCF: ${input_genotype}"
        echo "Reference VCF: ${reference_vcf}"

        echo "Target chromosomes:"
        bcftools query -f '%CHROM\\n' ${input_genotype} | sort -u | head -30 || true

        echo "Reference chromosomes:"
        bcftools query -f '%CHROM\\n' ${reference_vcf} | sort -u | head -30 || true

        eagle \\
            --vcfTarget=${input_genotype} \\
            --vcfRef=${reference_vcf} \\
            --geneticMapFile=${genetic_map} \\
            --chrom=${chromosome} \\
            --outPrefix=phased_${chromosome} \\
            2>&1 | tee eagle_${chromosome}.log

        tabix -f -p vcf phased_${chromosome}.vcf.gz

        echo "=== Eagle output check ==="
        ls -lh phased_${chromosome}.vcf.gz*
        bcftools query -f '%CHROM\\n' phased_${chromosome}.vcf.gz | head
        """
}