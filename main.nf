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
        path "genetic_map_all_chr.txt"  // Changed from .txt.gz to .txt

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
            print "chr"\$1, \$4, \$3
        }' OFS=' ' no_chr_in_chrom_field/plink.chr\${i}.GRCh38.map
    done > genetic_map_all_chr.txt

    # DO NOT compress - RFMix v2.03 has issues reading compressed genetic maps
    echo "Genetic map ready (uncompressed for RFMix compatibility)"
    """
}

// workflow download_only {
//     download_genetic_map()
// }



// ========== below is the cirro code



process download_genetic_map_eagle {
    
    output:
    path "genetic_map_hg38_withX.txt.gz"

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



workflow ancestry_pipeline {

    chr_ch = Channel.from(21..21)
    download_genetic_map()
    download_genetic_map_eagle()
    download_sample_map()
    // ensure the process emits a usable file path
    map_file_ch = download_genetic_map.out.flatten()
    map_file_eagle_ch = download_genetic_map_eagle.out.flatten()
    sample_file_ch = download_sample_map.out.flatten()

    // now combine chromosomes with the actual file
    eagle_inputs_ch = chr_ch.combine(map_file_eagle_ch).map { chr, map_file ->
        println "Preparing Eagle inputs for chromosome ${chr} with genetic map ${map_file}"
        def vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
        def vcf_index = "${vcf}.tbi"
        [
            file(params.input_genotype),
            file(params.input_genotype_index),
            file(vcf),                    // reference VCF
            file(vcf_index),              // reference VCF index
            file(map_file), 
            chr
        ]
    }

    phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

    phased_vcf_with_chr_ch = phased_vcf_ch
        .combine(eagle_inputs_ch.map { it[5] })  // combine with the chromosomes
        .map { phased_vcf_file, chr ->
            println "Eagle finished for chromosome ${chr}: phased VCF = ${phased_vcf_file}"
            tuple(chr, phased_vcf_file)
        }
    rfmix_inputs_ch = phased_vcf_with_chr_ch
        .combine(map_file_ch)
        .combine(sample_file_ch)
        .map { it.flatten() }   // now each element is [chr, phased_vcf, map_file, sample_file]
        .map { combined ->
            def chr = combined[0]
            def phased_vcf = combined[1]
            def map_file = combined[2]
            def sample_map = combined[3]

            def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
            def ref_vcf_index = "${ref_vcf}.tbi"

            tuple(
                chr,
                file(phased_vcf),
                file(ref_vcf),
                file(ref_vcf_index),
                file(map_file, copy: true),
                file(sample_map, copy: true)
            )
        }
   
    rfmix_results = run_rfmix(rfmix_inputs_ch)

    emit:
        rfmix_results
    }



// Default workflow for Cirro to run
workflow { ancestry_pipeline() }

