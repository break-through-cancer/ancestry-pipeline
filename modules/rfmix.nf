process run_rfmix {
    tag "RFMix Ancestry"

    input:
        tuple val(chromosome), file(phased_vcf), file(reference_vcf), file(reference_vcf_index), file(genetic_map), file(sample_map)

    output:
        path "${chromosome}.*"

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

        echo "--- Check all chromosomes in genetic map ---"
        if [[ ${genetic_map} == *.gz ]]; then
            zcat ${genetic_map} | awk '{print \$1}' | sort | uniq
        else
            awk '{print \$1}' ${genetic_map} | sort | uniq
        fi

        echo "--- Genetic map file size ---"
        ls -lh ${genetic_map}

        echo "--- Check for carriage returns in genetic map ---"
        zcat ${genetic_map} | od -c | grep '\r' || echo "No carriage returns found"

        echo "--- Check delimiter / special characters in first 10 lines ---"
        zcat ${genetic_map} | head -n 10 | sed -n l


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
