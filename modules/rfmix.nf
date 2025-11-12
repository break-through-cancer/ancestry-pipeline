process run_rfmix {
    tag "RFMix Ancestry"

    input:
        tuple val(chromosome), file(phased_vcf), file(reference_vcf), file(reference_vcf_index), file(genetic_map), file(sample_map)

    output:
        path "rfmix_${chromosome}_results/*"

    script:
        """
        echo "=== DEBUG INFO for chromosome ${chromosome} ==="
        echo "Phased VCF: ${phased_vcf}"
        echo "Reference VCF: ${reference_vcf}"
        echo "Genetic map: ${genetic_map}"
        echo "Sample map: ${sample_map}"

        print_head() {
            local file="$1"
            echo "--- Showing head of: ${file} ---"
            if [[ ${file} == *.gz ]]; then
                echo "(gzipped file — printing first few lines as text)"
                zcat "${file}" 2>/dev/null | head || echo "⚠️ could not read gzipped file"
            else
                file -b "${file}" | grep -q "text" && head "${file}" || echo "(binary file — skipped)"
            fi
            echo
        }

        print_head "${genetic_map}"
        print_head "${sample_map}"

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
