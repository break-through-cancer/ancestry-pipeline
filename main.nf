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

    input_genotype_ch = Channel.fromPath(params.input_genotype)
    genetic_map_ch = Channel.fromPath(params.genetic_map)
    sample_map_ch = Channel.fromPath(params.sample_map)

    chr_ch = Channel.from(1..22)

    rfmix_results = chr_ch.flatMap { chr ->

        // VCF path as channel
        ref_vcf_path = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
        reference_vcf_ch = Channel.fromPath(ref_vcf_path)

        // Call Eagle process via channels
        phased_vcf_ch = phase_with_eagle(
            input_genotype_ch,
            reference_vcf_ch,
            genetic_map_ch,
            Channel.value(chr)      // wrap single value as a channel
        )

        // Call RFMix process via channels
        rfmix_out_ch = run_rfmix(
            phased_vcf_ch,
            reference_vcf_ch,
            sample_map_ch,
            genetic_map_ch,
            Channel.value(chr)
        )

        return rfmix_out_ch
    }

    emit:
        rfmix_results
}

// Default workflow for Cirro to run
workflow { ancestry_pipeline() }
