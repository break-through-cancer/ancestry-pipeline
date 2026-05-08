#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { phase_with_eagle } from './modules/eagle'
include { run_rfmix } from './modules/rfmix'

if (!params.input_genotype) {
    exit 1, 'Please, provide an input genotype data !'
}

if (!params.input_genotype_index) {
    exit 1, 'Please, provide an input genotype index !'
}

process download_genetic_map {

    output:
        path "genetic_map_all_chr.txt", emit: map

    script:
    """
    set -euo pipefail

    echo "Downloading PLINK GRCh38 maps..."
    curl -fL --retry 5 --retry-delay 10 \\
        -o plink.GRCh38.map.zip \\
        https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip

    unzip -o plink.GRCh38.map.zip

    for i in {1..22}; do
        awk -v chr="chr\${i}" '{
            print chr, \$4, \$3
        }' OFS=' ' no_chr_in_chrom_field/plink.chr\${i}.GRCh38.map
    done > genetic_map_all_chr.txt

    echo "RFMix genetic map ready:"
    head genetic_map_all_chr.txt
    """
}

process download_genetic_map_eagle {

    output:
        path "genetic_map_hg38_withX.txt.gz", emit: map

    script:
    """
    set -euo pipefail

    curl -fL --retry 5 --retry-delay 10 \\
        -o genetic_map_hg38_withX.txt.gz \\
        https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz

    ls -lh genetic_map_hg38_withX.txt.gz
    """
}

process download_sample_map {

    output:
        path "rfmix_sample_map.txt", emit: sample_map

    script:
    """
    set -euo pipefail

    curl -fL --retry 5 --retry-delay 10 \\
        -o full_sample_map.txt \\
        https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

    awk 'NR > 1 {print \$2, \$6}' OFS='\\t' full_sample_map.txt > rfmix_sample_map.txt

    echo "RFMix sample map:"
    head rfmix_sample_map.txt
    echo "Sample map line count:"
    wc -l rfmix_sample_map.txt
    """
}

process normalize_chrom_names {
    tag "$vcf_file"

    input:
        tuple path(vcf_file), path(vcf_index)

    output:
        tuple path("normalized_${vcf_file.simpleName}.vcf.gz"),
              path("normalized_${vcf_file.simpleName}.vcf.gz.tbi"),
              emit: normalized

    script:
    """
    set -euo pipefail

    OUT_VCF="normalized_${vcf_file.simpleName}.vcf.gz"

    echo "Normalizing chromosome names for ${vcf_file}..."
    echo "Input index: ${vcf_index}"

    FIRST_CHROM=\$(bcftools view -H ${vcf_file} | head -1 | cut -f1 || true)
    echo "First chromosome found: \${FIRST_CHROM}"

    if [[ "\${FIRST_CHROM}" == chr* ]]; then
        echo "Chromosomes already chr-prefixed."
        bcftools view -Oz -o \${OUT_VCF} ${vcf_file}
    else
        echo "Converting chromosomes to chr-prefixed format..."

        bcftools query -f '%CHROM\\n' ${vcf_file} | sort -u | while read chrom; do
            if [[ "\${chrom}" == "MT" ]]; then
                echo -e "MT\\tchrM"
            elif [[ "\${chrom}" == "M" ]]; then
                echo -e "M\\tchrM"
            elif [[ "\${chrom}" == chr* ]]; then
                echo -e "\${chrom}\\t\${chrom}"
            else
                echo -e "\${chrom}\\tchr\${chrom}"
            fi
        done > reheader.txt

        echo "Chromosome rename map:"
        cat reheader.txt

        bcftools annotate \\
            --rename-chrs reheader.txt \\
            ${vcf_file} \\
            -Oz -o \${OUT_VCF}
    fi

    tabix -f -p vcf \${OUT_VCF}

    echo "Normalized output:"
    ls -lh \${OUT_VCF}*
    echo "Chromosomes in normalized output:"
    bcftools query -f '%CHROM\\n' \${OUT_VCF} | sort -u | head -30
    """
}

process download_1000g_phased_reference {
    tag "$chr_name"

    input:
        val chr_name

    output:
        tuple val(chr_name),
              path("ref_${chr_name}.vcf.gz"),
              path("ref_${chr_name}.vcf.gz.tbi"),
              emit: phased_ref

    script:
    def chr_num = chr_name.toString().replaceFirst(/^chr/, "")
    def base_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
    def remote_prefix = "CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr_num}.filtered.shapeit2-duohmm-phased.vcf.gz"

    """
    set -euo pipefail

    echo "=== Downloading already-phased 1000G hg38 reference for ${chr_name} ==="

    REF_URL="${base_url}/${remote_prefix}"
    TBI_URL="${base_url}/${remote_prefix}.tbi"

    echo "VCF URL: \${REF_URL}"
    echo "TBI URL: \${TBI_URL}"

    wget -c \
        --tries=20 \
        --waitretry=30 \
        --read-timeout=60 \
        --timeout=60 \
        -O ref_${chr_name}.vcf.gz \
        "${REF_URL}"

    wget -c \
        --tries=20 \
        --waitretry=30 \
        --read-timeout=60 \
        --timeout=60 \
        -O ref_${chr_name}.vcf.gz.tbi \
        "${TBI_URL}"

    echo "Downloaded reference:"
    ls -lh ref_${chr_name}.vcf.gz ref_${chr_name}.vcf.gz.tbi

    echo "Reference chromosomes:"
    bcftools query -f '%CHROM\\n' ref_${chr_name}.vcf.gz | sort -u | head -20

    echo "Reference sample count:"
    bcftools query -l ref_${chr_name}.vcf.gz | wc -l
    """
}

process convert_rfmix_inputs_to_bcf {
    tag "$chr_name"

    input:
        tuple val(chr_name),
              path(target_vcf),
              path(target_tbi),
              path(ref_vcf),
              path(ref_tbi)

    output:
        tuple val(chr_name),
              path("target_${chr_name}.bcf"),
              path("target_${chr_name}.bcf.csi"),
              path("ref_${chr_name}.bcf"),
              path("ref_${chr_name}.bcf.csi"),
              emit: bcf_inputs

    script:
    """
    set -euo pipefail

    echo "=== Converting RFMix inputs to BCF for ${chr_name} ==="

    echo "Target VCF: ${target_vcf}"
    echo "Reference VCF: ${ref_vcf}"

    echo "--- Target VCF chromosomes ---"
    bcftools query -f '%CHROM\\n' ${target_vcf} | sort -u | head

    echo "--- Reference VCF chromosomes ---"
    bcftools query -f '%CHROM\\n' ${ref_vcf} | sort -u | head

    echo "--- Converting target phased VCF to BCF ---"
    bcftools view -Ob -o target_${chr_name}.bcf ${target_vcf}
    bcftools index -f target_${chr_name}.bcf

    echo "--- Converting reference phased VCF to BCF ---"
    bcftools view -Ob -o ref_${chr_name}.bcf ${ref_vcf}
    bcftools index -f ref_${chr_name}.bcf

    echo "BCF outputs:"
    ls -lh target_${chr_name}.bcf target_${chr_name}.bcf.csi ref_${chr_name}.bcf ref_${chr_name}.bcf.csi

    echo "Target BCF chromosomes:"
    bcftools query -f '%CHROM\\n' target_${chr_name}.bcf | sort -u | head

    echo "Reference BCF chromosomes:"
    bcftools query -f '%CHROM\\n' ref_${chr_name}.bcf | sort -u | head

    echo "Target BCF sample count:"
    bcftools query -l target_${chr_name}.bcf | wc -l

    echo "Reference BCF sample count:"
    bcftools query -l ref_${chr_name}.bcf | wc -l
    """
}

workflow ancestry_pipeline {

    chr_ch = Channel.from(1..22).map { chr_num -> "chr${chr_num}" }

    download_genetic_map()
    download_genetic_map_eagle()
    download_sample_map()

    map_file_ch = download_genetic_map.out.map.flatten()
    map_file_eagle_ch = download_genetic_map_eagle.out.map.flatten()
    sample_file_ch = download_sample_map.out.sample_map.flatten()

    input_vcf_ch = Channel.of(
        tuple(
            file(params.input_genotype),
            file(params.input_genotype_index)
        )
    )

    normalized_vcf_ch = input_vcf_ch | normalize_chrom_names

    phased_ref_ch = download_1000g_phased_reference(chr_ch).phased_ref

    eagle_inputs_ch = phased_ref_ch
        .combine(normalized_vcf_ch.normalized)
        .combine(map_file_eagle_ch)
        .map { chr_name, ref_vcf, ref_tbi, target_vcf, target_tbi, eagle_map ->

            println "Preparing Eagle target phasing inputs for chromosome ${chr_name}"
            println "Target VCF: ${target_vcf}"
            println "Phased 1000G reference VCF: ${ref_vcf}"
            println "Eagle map: ${eagle_map}"

            tuple(
                file(target_vcf),
                file(target_tbi),
                file(ref_vcf),
                file(ref_tbi),
                file(eagle_map),
                chr_name
            )
        }

    phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

    rfmix_vcf_inputs_ch = phased_vcf_ch.phased_vcf
        .join(phased_ref_ch)
        .map { chr_name, target_phased_vcf, target_phased_tbi, ref_phased_vcf, ref_phased_tbi ->

            println "Preparing BCF conversion inputs for chromosome ${chr_name}"
            println "Target phased VCF: ${target_phased_vcf}"
            println "Reference phased VCF: ${ref_phased_vcf}"

            tuple(
                chr_name,
                file(target_phased_vcf),
                file(target_phased_tbi),
                file(ref_phased_vcf),
                file(ref_phased_tbi)
            )
        }

    rfmix_bcf_ch = convert_rfmix_inputs_to_bcf(rfmix_vcf_inputs_ch).bcf_inputs

    rfmix_inputs_ch = rfmix_bcf_ch
        .combine(map_file_ch)
        .combine(sample_file_ch)
        .map { chr_name, target_bcf, target_csi, ref_bcf, ref_csi, rfmix_map, sample_map ->

            println "Preparing RFMix BCF inputs for chromosome ${chr_name}"
            println "Target phased BCF: ${target_bcf}"
            println "Reference phased BCF: ${ref_bcf}"
            println "RFMix map: ${rfmix_map}"
            println "Sample map: ${sample_map}"

            tuple(
                chr_name,
                file(target_bcf),
                file(target_csi),
                file(ref_bcf),
                file(ref_csi),
                file(rfmix_map),
                file(sample_map)
            )
        }

    rfmix_results = run_rfmix(rfmix_inputs_ch)

    emit:
        rfmix_results
}

workflow {
    ancestry_pipeline()
}