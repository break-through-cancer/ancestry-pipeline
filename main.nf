#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { phase_with_eagle } from '/modules/eagle'
include { run_rfmix } from '/modules/rfmix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input_genotype) { input_genotype = params.input_genotype } else { exit 1, 'Please, provide an input genotype data !' }
if (params.reference_vcf) { reference_vcf = params.reference_vcf } else { exit 1, 'Please, provide a reference vcf!' }
if (params.genetic_map) { genetic_map = params.genetic_map } else { exit 1, 'Please, provide a genetic map !' }
if (params.sample_map) { sample_map = params.sample_map } else { exit 1, 'Please provide a sample map file' }
if (params.chromosome) { chromosome = params.chromosome } else { exit 1, ' Please provide a chromosome to analyze via --chromosome <chr1|chr2|...>' }
if (params.output_prefix) { output_prefix = params.output_prefix } else { output_prefix = "output" }

workflow start{

    main:
        //ch_samplesheet = Channel.fromPath(samplesheet, checkIfExists: true)
        //phase w eagle
        phased_vcf = phase_with_eagle(input_genotype, reference_vcf, genetic_map, chromosome, output_prefix)
        //run rfmix
        rfmix_out = run_rfmix(phased_vcf, reference_vcf, sample_map, genetic_map, chromosome, output_prefix)

    emit:
    rfmix_out
}