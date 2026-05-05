#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { phase_with_eagle } from './modules/eagle'
include { run_rfmix } from './modules/rfmix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input_genotype) {
    input_genotype = params.input_genotype
} else {
    exit 1, 'Please, provide an input genotype data !'
}

if (params.input_genotype_index) {
    input_genotype_input = params.input_genotype_index
} else {
    exit 1, 'Please, provide an input genotype index !'
}


process download_genetic_map {

    output:
        path "genetic_map_all_chr.txt", emit: map

    script:
    """
    echo "Downloading PLINK GRCh38 maps..."
    if [ ! -f plink.GRCh38.map.zip ]; then
        wget -q https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
    fi

    echo "Unzipping..."
    unzip -o plink.GRCh38.map.zip

    echo "Building combined genetic map with chr-prefixed chromosomes..."

    # PLINK map columns: chr snp cM pos
    # RFMix wants: chr pos cM
    # Force chr-prefixed chromosome names here.
    for i in {1..22}; do
        awk -v chr="chr\${i}" '{
            print chr, \$4, \$3
        }' OFS=' ' no_chr_in_chrom_field/plink.chr\${i}.GRCh38.map
    done > genetic_map_all_chr.txt

    echo "Genetic map ready with chr-prefixed chromosomes:"
    head genetic_map_all_chr.txt

    # DO NOT compress - RFMix v2.03 has issues reading compressed genetic maps
    echo "Genetic map ready (uncompressed for RFMix compatibility)"
    """
}


process download_genetic_map_eagle {

    output:
        path "genetic_map_hg38_withX.txt.gz", emit: map

    script:
    """
    echo "=== Starting download_genetic_map_eagle process ==="

    if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
        echo "Downloading Eagle genetic map with curl..."
        curl -s -L -o genetic_map_hg38_withX.txt.gz \\
            https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
    else
        echo "Eagle genetic map already exists, skipping download."
    fi
    """
}

process download_sample_map {

    output:
        path "rfmix_sample_map.txt", emit: sample_map

    script:
    """
    set -euo pipefail

    curl -s -L -o full_sample_map.txt \\
        https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/20130606_g1k_3202_samples_ped_population.txt

    awk 'NR > 1 {print \$2, \$6}' OFS='\\t' full_sample_map.txt > rfmix_sample_map.txt

    echo "RFMix sample map:"
    head rfmix_sample_map.txt
    wc -l rfmix_sample_map.txt
    """
}

process normalize_chrom_names {
    tag "$vcf_file"

    input:
        path vcf_file

    output:
        tuple path("normalized_${vcf_file.simpleName}.vcf.gz"),
              path("normalized_${vcf_file.simpleName}.vcf.gz.tbi"),
              emit: normalized

    script:
    """
    set -euo pipefail

    OUT_VCF="normalized_${vcf_file.simpleName}.vcf.gz"

    echo "Normalizing chromosome names for ${vcf_file}..."
    echo "Goal: all output chromosomes should use chr-prefix format, e.g. chr1, chr2, chr21."

    FIRST_CHROM=\$(zcat ${vcf_file} 2>/dev/null | grep -v '^#' | head -1 | cut -f1 || true)
    echo "First chromosome found: \${FIRST_CHROM}"

    if [[ "\${FIRST_CHROM}" == chr* ]]; then
        echo "Chromosomes already start with chr. Copying VCF without renaming..."
        bcftools view -Oz -o \${OUT_VCF} ${vcf_file}
    else
        echo "Chromosomes are numeric/non-chr format. Converting to chr-prefixed format..."

        # Build chromosome renaming map:
        # 1  -> chr1
        # 2  -> chr2
        # 21 -> chr21
        # MT -> chrM
        (zcat ${vcf_file} | grep -v '^#' | cut -f1 | sort -u || true) | while read chrom; do
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

        echo "Chromosome mapping:"
        cat reheader.txt

        bcftools annotate \\
            --rename-chrs reheader.txt \\
            ${vcf_file} \\
            -Oz -o \${OUT_VCF}
    fi

    tabix -f -p vcf \${OUT_VCF}

    echo "Output files created:"
    ls -lh \${OUT_VCF}* || true

    echo "Chromosomes in normalized output:"
    bcftools query -f '%CHROM\\n' \${OUT_VCF} | sort -u | head -30 || true

    echo "Contig lines in header:"
    bcftools view -h \${OUT_VCF} | grep '^##contig' | head -20 || true
    """
}


workflow ancestry_pipeline {

    chr_ch = Channel.from(1..22)

    download_genetic_map()
    download_genetic_map_eagle()
    download_sample_map()

    map_file_ch = download_genetic_map.out.map.flatten()
    map_file_eagle_ch = download_genetic_map_eagle.out.map.flatten()
    sample_file_ch = download_sample_map.out.sample_map.flatten()

    map_file_ch.view { println "Map for RFMix: $it" }

    input_vcf_ch = Channel.fromPath(params.input_genotype)

    normalized_vcf_ch = input_vcf_ch | normalize_chrom_names

    eagle_inputs_ch = chr_ch
        .combine(normalized_vcf_ch.normalized)
        .combine(map_file_eagle_ch)
        .map { chr_num, vcf_file, vcf_index, map_file ->

            def chr_name = "chr${chr_num}"

            println "Preparing Eagle inputs for chromosome ${chr_name} with genetic map ${map_file}"

            def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr_num}.recalibrated_variants.vcf.gz"
            def ref_vcf_index = "${ref_vcf}.tbi"
            tuple(
                file(vcf_file),
                file(vcf_index),
                file(ref_vcf),
                file(ref_vcf_index),
                file(map_file),
                chr_name
            )
        }

    phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

    rfmix_inputs_ch = phased_vcf_ch.phased_vcf
        .combine(map_file_ch)
        .combine(sample_file_ch)
        .map { chr_name, phased_vcf, phased_tbi, map_file, sample_map ->

            def chr_num = chr_name.toString().replaceFirst(/^chr/, "")

            def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr_num}.filtered.shapeit2-duohmm-phased.vcf.gz"
            def ref_vcf_index = "${ref_vcf}.tbi"

            tuple(
                chr_name,
                phased_vcf,
                phased_tbi,
                file(ref_vcf),
                file(ref_vcf_index),
                file(map_file),
                file(sample_map)
            )
        }

    rfmix_results = run_rfmix(rfmix_inputs_ch)

    emit:
        rfmix_results
}


// Default workflow for Cirro to run
workflow {
    ancestry_pipeline()
}