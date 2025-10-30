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
// if (params.reference_vcf) { reference_vcf = params.reference_vcf } else { exit 1, 'Please, provide a reference vcf!' }
//if (params.genetic_map) { genetic_map = params.genetic_map } else { exit 1, 'Please, provide a genetic map !' }
// if (params.sample_map) { sample_map = params.sample_map } else { exit 1, 'Please provide a sample map file' }
// if (params.chromosome) { chromosome = params.chromosome } else { exit 1, ' Please provide a chromosome to analyze via --chromosome <chr1|chr2|...>' }
//if (params.output_prefix) { output_prefix = params.output_prefix } else { output_prefix = "output" }

process download_genetic_map {
    
    output:
    path "genetic_map_hg38_withX.txt.gz"

    script:
    """
    if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
        echo "Downloading genetic map with curl..."
        curl -s -L -o genetic_map_hg38_withX.txt.gz \
        https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
    else
        echo "Genetic map already exists, skipping download."
    fi
    """
}


workflow ancestry_pipeline {

    chr_ch = Channel.from(1..22)
    download_genetic_map()          
    // ensure the process emits a usable file path
    map_file_ch = download_genetic_map.out.flatten()

    // now combine chromosomes with the actual file
    eagle_inputs_ch = chr_ch.combine(map_file_ch).map { chr, map_file ->
        [
            file(params.input_genotype),
            file("s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"),
            file(map_file),  // ✅ now a resolved path, not a DataflowVariable
            chr
        ]
    }

    phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

    rfmix_results = run_rfmix(phased_vcf_ch, file(params.sample_map), map_file_ch)

    emit:
        rfmix_results
}



// Default workflow for Cirro to run
workflow { ancestry_pipeline() }


