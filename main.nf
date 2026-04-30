#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { phase_with_eagle } from './modules/eagle'
include { run_rfmix } from './modules/rfmix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input_genotype) { input_genotype = params.input_genotype } else { exit 1, 'Please, provide an input genotype data !' }
if (params.input_genotype_index) { input_genotype_input = params.input_genotype_index } else { exit 1, 'Please, provide an input genotype index !' }
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

    echo "Building combined genetic map..."

    # Combine all chromosomes into one file
    for i in {1..22}; do
        # PLINK map columns: chr snp cM pos
        # RFMix wants: chr pos cM
        awk '{
            print \$1, \$4, \$3
        }' OFS=' ' no_chr_in_chrom_field/plink.chr\${i}.GRCh38.map
    done > genetic_map_all_chr.txt

    # DO NOT compress - RFMix v2.03 has issues reading compressed genetic maps
    echo "Genetic map ready (uncompressed for RFMix compatibility)"
    """
}



process download_genetic_map_eagle {
    
    output:
    path "genetic_map_hg38_withX.txt.gz", emit: map
    script:
    """
    echo "=== Starting download_genetic_map process ==="
    if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
        echo "Downloading genetic map with curl..."
        curl -s -L -o genetic_map_hg38_withX.txt.gz \
        https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
    else
        echo "Genetic map already exists, skipping download."
    fi
    """
}

process download_sample_map {

    output:
        path "rfmix_sample_map.txt"

    script:
    """
    echo "=== Starting download_sample_map process ==="
    if [ ! -f integrated_call_samples_v3.20130502.ALL.panel ]; then
        echo "Downloading 1000 Genomes sample metadata..."
        curl -s -L -o integrated_call_samples_v3.20130502.ALL.panel \\
        https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
    else
        echo "Panel file already exists, skipping download."
    fi

    python3 - <<'EOF'
    import pandas as pd

    panel_file = "integrated_call_samples_v3.20130502.ALL.panel"
    df = pd.read_csv(panel_file, sep='\\t')
    sample_map = df[['sample', 'pop']]
    sample_map.to_csv("rfmix_sample_map.txt", sep='\\t', index=False, header=False)
    print("✅ Sample map saved to rfmix_sample_map.txt")
    EOF
        """
    }


process normalize_chrom_names {
    tag "$vcf_file"
    
    input:
    path vcf_file
    
    output:
    tuple path("normalized_${vcf_file.simpleName}.vcf.gz"), path("normalized_${vcf_file.simpleName}.vcf.gz.tbi"), emit: normalized
    
    script:
    """
    set -e
    
    OUT_VCF="normalized_${vcf_file.simpleName}.vcf.gz"
    
    echo "Normalizing chromosome names for ${vcf_file}..."
    
    # Check first data line for chr prefix
    FIRST_CHROM=\$(zcat ${vcf_file} 2>/dev/null | grep -v '^#' | head -1 | cut -f1 || true)
    echo "First chromosome found: \${FIRST_CHROM}"
    
    if [[ "\${FIRST_CHROM}" == chr* ]]; then
        echo "Chromosomes start with 'chr', converting to numbers..."
        
        # Create chromosome renaming map
        (zcat ${vcf_file} | grep -v '^#' | cut -f1 | sort -u | grep -v '^chr\$' || true) | while read chrom; do
            new_chrom=\$(echo "\${chrom}" | sed 's/^chr//' | sed 's/^M\$/MT/')
            if [[ -n "\${new_chrom}" ]]; then
                echo -e "\${chrom}\\t\${new_chrom}"
            fi
        done > reheader.txt
        
        echo "Chromosome mapping:"
        cat reheader.txt
        
        # Extract and fix header
        bcftools view -h ${vcf_file} | \\
        sed 's/^##contig=<ID=chr/##contig=<ID=/' | \\
        sed 's/^##contig=<ID=M,/##contig=<ID=MT,/' > temp_header.vcf
        
        # Extract data with renamed chromosomes
        bcftools annotate --rename-chrs reheader.txt ${vcf_file} 2>/dev/null | \\
        bcftools view -H >> temp_header.vcf
        
        # Compress the result
        bgzip -c temp_header.vcf > \${OUT_VCF}
        
        rm temp_header.vcf
    else
        echo "Chromosomes are already numeric (found: \${FIRST_CHROM}), just copying..."
        bcftools view -Oz -o \${OUT_VCF} ${vcf_file}
    fi
    
    # Index with tabix
    tabix -p vcf \${OUT_VCF}
    
    # Verify outputs
    echo "Output files created:"
    ls -lh \${OUT_VCF}* || true
    
    echo "Chromosomes in normalized output:"
    (zcat \${OUT_VCF} 2>/dev/null | grep -v '^#' | cut -f1 | sort -u | head -20) || true
    
    echo "Contig lines in header:"
    (zcat \${OUT_VCF} 2>/dev/null | grep '^##contig' | head -10) || true
    """
}
workflow ancestry_pipeline {

 
    chr_ch = Channel.from(1..22)
    download_genetic_map()
    download_genetic_map_eagle()
    download_sample_map()
    // ensure the process emits a usable file path
    map_file_ch = download_genetic_map.out.flatten()
    map_file_eagle_ch = download_genetic_map_eagle.out.flatten()
    sample_file_ch = download_sample_map.out.flatten()
    map_file_ch.view { println "Map for RFMix: $it" }

    // Preprocess input VCF once
    // Make input genotype a channel
    input_vcf_ch = Channel.fromPath(params.input_genotype)

    // Run the normalization process
    normalized_vcf_ch = input_vcf_ch | normalize_chrom_names
    eagle_inputs_ch = chr_ch
        .combine(normalized_vcf_ch.normalized)
        .combine(map_file_eagle_ch)
        .map { chr, vcf_file, vcf_index, map_file ->
            println "Preparing Eagle inputs for chromosome ${chr} with genetic map ${map_file}"
            def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
            def ref_vcf_index = "${ref_vcf}.tbi"
            tuple(
                file(vcf_file),
                file(vcf_index),
                file(ref_vcf),
                file(ref_vcf_index),
                file(map_file), 
                chr
            )
        }

    phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

    // phased_vcf_with_chr_ch = phased_vcf_ch
    //     .combine(eagle_inputs_ch.map { it[5] })  // combine with the chromosomes
    //     .map { phased_vcf_file, chr ->
    //         println "Eagle finished for chromosome ${chr}: phased VCF = ${phased_vcf_file}"
    //         tuple(chr, phased_vcf_file)
    //     }
    rfmix_inputs_ch = phased_vcf_ch.phased_vcf
        .combine(map_file_ch)
        .combine(sample_file_ch)
        .map { chr, phased_vcf, phased_tbi, map_file, sample_map ->

            def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
            def ref_vcf_index = "${ref_vcf}.tbi"

            tuple(
                chr,
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
workflow { ancestry_pipeline() }

