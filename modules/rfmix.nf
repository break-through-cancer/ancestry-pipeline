process run_rfmix {
    tag "RFMix ${chromosome}"

    input:
        tuple val(chromosome),
              path(phased_vcf),
              path(phased_vcf_index),
              path(reference_vcf),
              path(reference_vcf_index),
              path(genetic_map),
              path(sample_map)

    output:
        path "${chromosome}.*"

    script:
        """
        set -euo pipefail

        echo "=== DEBUG INFO for chromosome ${chromosome} ==="
        echo "Phased VCF: ${phased_vcf}"
        echo "Reference VCF: ${reference_vcf}"
        echo "Genetic map: ${genetic_map}"
        echo "Sample map: ${sample_map}"

        echo "--- Chromosomes in phased VCF ---"
        bcftools query -f '%CHROM\\n' ${phased_vcf} | sort -u | head

        echo "--- Chromosomes in reference VCF ---"
        bcftools query -f '%CHROM\\n' ${reference_vcf} | sort -u | head

        echo "--- First lines of genetic map ---"
        head ${genetic_map}

        echo "--- Running RFMix ---"

        rfmix \\
            -f ${phased_vcf} \\
            -r ${reference_vcf} \\
            -m ${sample_map} \\
            -g ${genetic_map} \\
            -o ${chromosome} \\
            --chromosome=${chromosome}
        """
}