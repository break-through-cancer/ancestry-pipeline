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
if (params.genetic_map) { genetic_map = params.genetic_map } else { exit 1, 'Please, provide a genetic map !' }
// if (params.sample_map) { sample_map = params.sample_map } else { exit 1, 'Please provide a sample map file' }
// if (params.chromosome) { chromosome = params.chromosome } else { exit 1, ' Please provide a chromosome to analyze via --chromosome <chr1|chr2|...>' }
//if (params.output_prefix) { output_prefix = params.output_prefix } else { output_prefix = "output" }


workflow ancestry_pipeline {

    chr_ch = Channel.from(1..22)

    eagle_inputs_ch = chr_ch.map { chr ->
        [
            file(params.input_genotype),
            file("s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"),
            file(params.genetic_map),
            chr
        ]
    }

    phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

    rfmix_results = run_rfmix(phased_vcf_ch, file(params.sample_map), file(params.genetic_map))

    emit:
        rfmix_results
}


// Default workflow for Cirro to run
workflow { ancestry_pipeline() }
