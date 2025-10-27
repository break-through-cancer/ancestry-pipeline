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
// if (params.reference_vcf) { reference_vcf = params.reference_vcf } else { exit 1, 'Please, provide a reference vcf!' }
if (params.genetic_map) { genetic_map = params.genetic_map } else { exit 1, 'Please, provide a genetic map !' }
// if (params.sample_map) { sample_map = params.sample_map } else { exit 1, 'Please provide a sample map file' }
// if (params.chromosome) { chromosome = params.chromosome } else { exit 1, ' Please provide a chromosome to analyze via --chromosome <chr1|chr2|...>' }
//if (params.output_prefix) { output_prefix = params.output_prefix } else { output_prefix = "output" }

workflow start {

    chr_ch = Channel.from(1..22)
    
    rfmix_results = chr_ch.flatMap { chr ->

        ref_vcf_path = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
        reference_vcf_for_eagle = Channel.fromPath(ref_vcf_path)
        reference_vcf_for_rfmix = Channel.fromPath(ref_vcf_path)

        phased_vcf = phase_with_eagle(input_genotype, reference_vcf_for_eagle, genetic_map, chr)
        rfmix_out = run_rfmix(phased_vcf, reference_vcf_for_rfmix, sample_map, genetic_map, chr)

        return rfmix_out
    }

    emit:
        rfmix_results
}

// workflow start {

//     main:
//         rfmix_results = Channel.empty() 

//         for chr in 1..22 {
//             echo "Processing chr${chr}..."
            
//             ref_vcf_path = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"

//             reference_vcf_for_eagle = Channel.fromPath(ref_vcf_path)
//             reference_vcf_for_rfmix = Channel.fromPath(ref_vcf_path)
//             //phase w eagle
//             phased_vcf = phase_with_eagle(input_genotype, reference_vcf_for_eagle, genetic_map, chr)
//             //run rfmix
//             rfmix_out = run_rfmix(phased_vcf, reference_vcf_for_rfmix, sample_map, genetic_map, chr)
//             rfmix_results = rfmix_results.mix(rfmix_out) // append each iteration
//         }

//     emit:
//         rfmix_results
// }