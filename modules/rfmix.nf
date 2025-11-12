process run_rfmix {
    tag "RFMix Ancestry"

    input:
        tuple val(chromosome), file(phased_vcf), file(reference_vcf), file(reference_vcf_index), file(genetic_map), file(sample_map)

    output:
        path "rfmix_${chromosome}_results/*"

    script:
        """
        #!/bin/bash

        echo "=== DEBUG INFO for chromosome ${chromosome} ==="
        echo "Phased VCF: ${phased_vcf}"
        echo "Reference VCF: ${reference_vcf}"
        echo "Genetic map: ${genetic_map}"
        echo "Sample map: ${sample_map}"

        echo "--- Head of genetic map ---"
        if [[ ${genetic_map} == *.gz ]]; then
            zcat ${genetic_map} | head -n 5
        else
            head -n 5 ${genetic_map}
        fi

        echo "--- Head of sample map ---"
        head -n 5 ${sample_map}

        echo "--- Running RFMix ---"
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
