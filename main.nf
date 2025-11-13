#!/usr/bin/env nextflow
// nextflow.enable.dsl = 2
include { phase_with_eagle } from './modules/eagle'
include { run_rfmix } from './modules/rfmix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input_genotype) { input_genotype = params.input_genotype } else { exit 1, 'Please, provide an input genotype data !' }
if (params.input_genotype_index) { input_genotype = params.input_genotype_index } else { exit 1, 'Please, provide an input genotype index !' }
// if (params.reference_vcf) { reference_vcf = params.reference_vcf } else { exit 1, 'Please, provide a reference vcf!' }
//if (params.genetic_map) { genetic_map = params.genetic_map } else { exit 1, 'Please, provide a genetic map !' }
if (params.sample_map) { sample_map = params.sample_map } else { exit 1, 'Please provide a sample map file' }
// if (params.chromosome) { chromosome = params.chromosome } else { exit 1, ' Please provide a chromosome to analyze via --chromosome <chr1|chr2|...>' }
//if (params.output_prefix) { output_prefix = params.output_prefix } else { output_prefix = "output" }

process download_genetic_map {
    
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
    echo "Adding chr prefix to genetic map..."
    zcat genetic_map_hg38_withX.txt.gz | sed -E '1!s/^([0-9XYM])/chr\1/' | gzip > genetic_map_chr.txt.gz
    echo "Genetic map prepared: genetic_map_chr.txt.gz"
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

    chr_ch = Channel.from(1..22)
    download_genetic_map()
    download_sample_map()
    // ensure the process emits a usable file path
    map_file_ch = download_genetic_map.out.flatten()
    sample_file_ch = download_sample_map.out.flatten()

    // now combine chromosomes with the actual file
    eagle_inputs_ch = chr_ch.combine(map_file_ch).map { chr, map_file ->
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

    rfmix_inputs_ch = phased_vcf_with_chr_ch.flatMap { chr, phased_vcf ->

        def tuples = []

        // collect the single genetic map file(s) into a list
        map_file_ch.toList().each { map_file ->

            // collect the single sample map file(s) into a list
            sample_file_ch.toList().each { sample_map ->

                def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
                def ref_vcf_index = "${ref_vcf}.tbi"

                tuples << tuple(
                    chr,
                    file(phased_vcf),
                    file(ref_vcf),
                    file(ref_vcf_index),
                    file(map_file),
                    file(sample_map)
                )
            }
        }
    }

    // return tuples
    // phased_vcf_with_chr_ch
    //     .combine(map_file_ch)
    //     .combine(sample_file_ch)
    //     .map { item, sample_map ->
    //         def (chr, phased_vcf, map_file) = item
    //         def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
    //         def ref_vcf_index = "${ref_vcf}.tbi"
    //         tuple(
    //             chr,
    //             file(phased_vcf),
    //             file(ref_vcf),
    //             file(ref_vcf_index),
    //             file(map_file),
    //             file(sample_map)
    //         )
    //     }


    rfmix_results = run_rfmix(rfmix_inputs_ch)

    emit:
        rfmix_results
    }



// Default workflow for Cirro to run
workflow { ancestry_pipeline() }


