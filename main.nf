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
if (params.reference_vcf) { reference_vcf = params.reference_vcf } else { exit 1, 'Please, provide a reference vcf!' }
if (params.reference_vcf_index) { reference_vcf_index = params.reference_vcf_index } else { exit 1, 'Please, provide a reference vcf ref!' }
//if (params.genetic_map) { genetic_map = params.genetic_map } else { exit 1, 'Please, provide a genetic map !' }
//if (params.sample_map) { sample_map = params.sample_map } else { exit 1, 'Please provide a sample map file' }
// if (params.chromosome) { chromosome = params.chromosome } else { exit 1, ' Please provide a chromosome to analyze via --chromosome <chr1|chr2|...>' }
//if (params.output_prefix) { output_prefix = params.output_prefix } else { output_prefix = "output" }
process download_genetic_map {

    output:
        path "genetic_map_chr_cleaned.txt.gz"

    script:
    """
    echo "=== Starting download_genetic_map process ==="

    # -------------------------------------------------------
    # 1. Download original genetic map if missing
    # -------------------------------------------------------
    if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
        echo "Downloading hg38 genetic map..."
        curl -s -L -o genetic_map_hg38_withX.txt.gz \
            https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
    fi

    echo "Decompressing..."
    gunzip -c genetic_map_hg38_withX.txt.gz > genetic_map_raw.txt

    # -------------------------------------------------------
    # 2. Convert chromosome numbers → chr1, chr2, ..., chrX
    # -------------------------------------------------------
    echo "Running Python chromosome conversion..."
    python3 - << 'EOF'
import sys

INPUT  = "genetic_map_raw.txt"
OUTPUT = "genetic_map_chr.txt"

VALID_AUTOSOMES = {str(i) for i in range(1, 23)}
SPECIAL = {"23": "X", "24": "Y"}

def convert(line):
    parts = line.strip().split()
    if not parts:
        return line

    chrom = parts[0]

    # Header line?
    if chrom.lower() in ["chromosome", "chr", "position", "pos"]:
        return line

    # Chromosome mapping
    if chrom in SPECIAL:
        parts[0] = "chr" + SPECIAL[chrom]
    elif chrom in VALID_AUTOSOMES:
        parts[0] = "chr" + chrom

    return "\\t".join(parts)

with open(INPUT) as fin, open(OUTPUT, "w") as fout:
    for line in fin:
        fout.write(convert(line) + "\\n")
EOF


    # -------------------------------------------------------
    # 3. Enforce strictly increasing columns (RFMix requirement)
    # -------------------------------------------------------
    echo "Cleaning map to ensure strictly increasing positions..."

    python3 - << 'EOF'
import pandas as pd

df = pd.read_csv("genetic_map_chr.txt", sep="\\t", header=None)

# Expect structure: chr, pos, rate(cM/Mb), genetic_map(cM)
# Columns 2 and 3 (indexes 2,3) must be strictly increasing

epsilon = 1e-6
prev2, prev3 = df.iloc[0, 2], df.iloc[0, 3]

for i in range(1, len(df)):
    if df.iloc[i, 2] <= prev2:
        df.iloc[i, 2] = prev2 + epsilon
    if df.iloc[i, 3] <= prev3:
        df.iloc[i, 3] = prev3 + epsilon
    prev2, prev3 = df.iloc[i, 2], df.iloc[i, 3]

df.to_csv("genetic_map_chr_cleaned.txt", sep="\\t", index=False, header=False, float_format="%.12f")
print("Cleaning done!")
EOF


    # -------------------------------------------------------
    # 4. Compress final cleaned map
    # -------------------------------------------------------
    echo "Compressing final cleaned map..."
    gzip -c genetic_map_chr_cleaned.txt > genetic_map_chr_cleaned.txt.gz

    echo "=== Done! Preview first 10 lines ==="
    zcat genetic_map_chr_cleaned.txt.gz | head
    """
}

// process download_genetic_map {

//     output:
//         path "genetic_map_chr.txt.gz"

//     script:
//     """
//     echo "=== Starting download_genetic_map process ==="

//     # 1. Download original genetic map if missing
//     if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
//         echo "Downloading genetic map..."
//         curl -s -L -o genetic_map_hg38_withX.txt.gz \
//             https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
//     fi

//     echo "Decompressing..."
//     gunzip -c genetic_map_hg38_withX.txt.gz > genetic_map_raw.txt

//     echo "Running Python chromosome conversion..."
//     python3 - << 'EOF'
// import sys

// INPUT = "genetic_map_raw.txt"
// OUTPUT = "genetic_map_chr.txt"

// VALID_AUTOSOMES = {str(i) for i in range(1, 23)}
// SPECIAL = {"23": "X", "24": "Y"}

// def convert(line):
//     parts = line.strip().split()
//     if not parts:
//         return line

//     chrom = parts[0]

//     # header?
//     if chrom.lower() in ["chromosome", "chr", "position", "pos"]:
//         return line

//     if chrom in SPECIAL:
//         parts[0] = "chr" + SPECIAL[chrom]
//         return "\t".join(parts)

//     if chrom in VALID_AUTOSOMES:
//         parts[0] = "chr" + chrom
//         return "\t".join(parts)

//     return line

// with open(INPUT) as fin, open(OUTPUT, "w") as fout:
//     for line in fin:
//         fout.write(convert(line) + "\\n")
// EOF

//     echo "Compressing final map..."
//     gzip -c genetic_map_chr.txt > genetic_map_chr.txt.gz

//     echo "=== Preview first 10 lines ==="
//     zcat genetic_map_chr.txt.gz
//     """
// }


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
    map_file_ch.view { println "Map for RFMix: $it" }

    // now combine chromosomes with the actual file
    eagle_inputs_ch = chr_ch.combine(map_file_eagle_ch).map { chr, map_file ->
        println "Preparing Eagle inputs for chromosome ${chr} with genetic map ${map_file}"
        // def vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
        // def vcf_index = "${vcf}.tbi"
        [
            file(params.input_genotype),
            file(params.input_genotype_index),
            file(params.reference_vcf),                    // reference VCF
            file(params.reference_vcf_index),              // reference VCF index
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

            // def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
            // def ref_vcf_index = "${ref_vcf}.tbi"

            tuple(
                chr,
                file(phased_vcf),
                file(params.reference_vcf),
                file(params.reference_vcf_index),
                file(map_file, copy: true),
                file(sample_map, copy: true)
            )
        }
    // phased_vcf_with_chr_ch.view { println "phased_vcf_with_chr_ch = $it" }
    // map_file_ch.view { println "map_file_ch = $it" }
    // sample_file_ch.view { println "sample_file_ch = $it" }

    // // rfmix_inputs_ch = phased_vcf_with_chr_ch
    // //     .combine(map_file_ch)
    // //     .combine(sample_file_ch)

    // rfmix_inputs_ch.view { println "rfmix_inputs_ch = $it" }



    // rfmix_inputs_ch = phased_vcf_with_chr_ch.flatMap { chr, phased_vcf ->

    //     def tuples = []

    //     // collect the single genetic map file(s) into a list
    //     map_file_ch.toList().each { map_file ->

    //         // collect the single sample map file(s) into a list
    //         sample_file_ch.toList().each { sample_map ->

    //             def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
    //             def ref_vcf_index = "${ref_vcf}.tbi"

    //             tuples << tuple(
    //                 chr,
    //                 file(phased_vcf),
    //                 file(ref_vcf),
    //                 file(ref_vcf_index),
    //                 file(map_file),
    //                 file(sample_map)
    //             )
    //         }
    //     }
    //     return tuples
    // }

    // 
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


// workflow run_download_genetic_map {
//     download_genetic_map()
// }

// ========== below is the cirro code


// #!/usr/bin/env nextflow
// // nextflow.enable.dsl = 2
// include { phase_with_eagle } from './modules/eagle'
// include { run_rfmix } from './modules/rfmix'

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Check mandatory parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

// if (params.input_genotype) { input_genotype = params.input_genotype } else { exit 1, 'Please, provide an input genotype data !' }
// if (params.input_genotype_index) { input_genotype_input = params.input_genotype_index } else { exit 1, 'Please, provide an input genotype index !' }
// // if (params.reference_vcf) { reference_vcf = params.reference_vcf } else { exit 1, 'Please, provide a reference vcf!' }
// //if (params.genetic_map) { genetic_map = params.genetic_map } else { exit 1, 'Please, provide a genetic map !' }
// if (params.sample_map) { sample_map = params.sample_map } else { exit 1, 'Please provide a sample map file' }
// // if (params.chromosome) { chromosome = params.chromosome } else { exit 1, ' Please provide a chromosome to analyze via --chromosome <chr1|chr2|...>' }
// //if (params.output_prefix) { output_prefix = params.output_prefix } else { output_prefix = "output" }

// process download_genetic_map {

//     output:
//         path "genetic_map_chr.txt.gz"

//     script:
//     """
//     echo "=== Starting download_genetic_map process ==="

//     # 1. Download original genetic map if missing
//     if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
//         echo "Downloading genetic map..."
//         curl -s -L -o genetic_map_hg38_withX.txt.gz \
//             https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
//     fi

//     echo "Decompressing..."
//     gunzip -c genetic_map_hg38_withX.txt.gz > genetic_map_raw.txt

//     echo "Running Python chromosome conversion..."
//     python3 - << 'EOF'
// import sys

// INPUT = "genetic_map_raw.txt"
// OUTPUT = "genetic_map_chr.txt"

// VALID_AUTOSOMES = {str(i) for i in range(1, 23)}
// SPECIAL = {"23": "X", "24": "Y"}

// def convert(line):
//     parts = line.strip().split()
//     if not parts:
//         return line

//     chrom = parts[0]

//     # header?
//     if chrom.lower() in ["chromosome", "chr", "position", "pos"]:
//         return line

//     if chrom in SPECIAL:
//         parts[0] = "chr" + SPECIAL[chrom]
//         return "\t".join(parts)

//     if chrom in VALID_AUTOSOMES:
//         parts[0] = "chr" + chrom
//         return "\t".join(parts)

//     return line

// with open(INPUT) as fin, open(OUTPUT, "w") as fout:
//     for line in fin:
//         fout.write(convert(line) + "\\n")
// EOF

//     echo "Compressing final map..."
//     gzip -c genetic_map_chr.txt > genetic_map_chr.txt.gz

//     echo "=== Preview first 10 lines ==="
//     zcat genetic_map_chr.txt.gz
//     """
// }


// process download_genetic_map_eagle {
    
//     output:
//     path "genetic_map_hg38_withX.txt.gz"

//     script:
//     """
//     echo "=== Starting download_genetic_map process ==="
//     if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
//         echo "Downloading genetic map with curl..."
//         curl -s -L -o genetic_map_hg38_withX.txt.gz \
//         https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
//     else
//         echo "Genetic map already exists, skipping download."
//     fi
//     """
// }

// process download_sample_map {

//     output:
//         path "rfmix_sample_map.txt"

//     script:
//     """
//     echo "=== Starting download_sample_map process ==="
//     if [ ! -f integrated_call_samples_v3.20130502.ALL.panel ]; then
//         echo "Downloading 1000 Genomes sample metadata..."
//         curl -s -L -o integrated_call_samples_v3.20130502.ALL.panel \\
//         https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
//     else
//         echo "Panel file already exists, skipping download."
//     fi

//     python3 - <<'EOF'
//     import pandas as pd

//     panel_file = "integrated_call_samples_v3.20130502.ALL.panel"
//     df = pd.read_csv(panel_file, sep='\\t')
//     sample_map = df[['sample', 'pop']]
//     sample_map.to_csv("rfmix_sample_map.txt", sep='\\t', index=False, header=False)
//     print("✅ Sample map saved to rfmix_sample_map.txt")
//     EOF
//         """
//     }



// workflow ancestry_pipeline {

//     chr_ch = Channel.from(21..21)
//     download_genetic_map()
//     download_genetic_map_eagle()
//     download_sample_map()
//     // ensure the process emits a usable file path
//     map_file_ch = download_genetic_map.out.flatten()
//     map_file_eagle_ch = download_genetic_map_eagle.out.flatten()
//     sample_file_ch = download_sample_map.out.flatten()

//     // now combine chromosomes with the actual file
//     eagle_inputs_ch = chr_ch.combine(map_file_eagle_ch).map { chr, map_file ->
//         println "Preparing Eagle inputs for chromosome ${chr} with genetic map ${map_file}"
//         def vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
//         def vcf_index = "${vcf}.tbi"
//         [
//             file(params.input_genotype),
//             file(params.input_genotype_index),
//             file(vcf),                    // reference VCF
//             file(vcf_index),              // reference VCF index
//             file(map_file), 
//             chr
//         ]
//     }

//     phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

//     phased_vcf_with_chr_ch = phased_vcf_ch
//         .combine(eagle_inputs_ch.map { it[5] })  // combine with the chromosomes
//         .map { phased_vcf_file, chr ->
//             println "Eagle finished for chromosome ${chr}: phased VCF = ${phased_vcf_file}"
//             tuple(chr, phased_vcf_file)
//         }
//     rfmix_inputs_ch = phased_vcf_with_chr_ch
//         .combine(map_file_ch)
//         .combine(sample_file_ch)
//         .map { it.flatten() }   // now each element is [chr, phased_vcf, map_file, sample_file]
//         .map { combined ->
//             def chr = combined[0]
//             def phased_vcf = combined[1]
//             def map_file = combined[2]
//             def sample_map = combined[3]

//             def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
//             def ref_vcf_index = "${ref_vcf}.tbi"

//             tuple(
//                 chr,
//                 file(phased_vcf),
//                 file(ref_vcf),
//                 file(ref_vcf_index),
//                 file(map_file),
//                 file(sample_map)
//             )
//         }
//     // phased_vcf_with_chr_ch.view { println "phased_vcf_with_chr_ch = $it" }
//     // map_file_ch.view { println "map_file_ch = $it" }
//     // sample_file_ch.view { println "sample_file_ch = $it" }

//     // // rfmix_inputs_ch = phased_vcf_with_chr_ch
//     // //     .combine(map_file_ch)
//     // //     .combine(sample_file_ch)

//     // rfmix_inputs_ch.view { println "rfmix_inputs_ch = $it" }



//     // rfmix_inputs_ch = phased_vcf_with_chr_ch.flatMap { chr, phased_vcf ->

//     //     def tuples = []

//     //     // collect the single genetic map file(s) into a list
//     //     map_file_ch.toList().each { map_file ->

//     //         // collect the single sample map file(s) into a list
//     //         sample_file_ch.toList().each { sample_map ->

//     //             def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
//     //             def ref_vcf_index = "${ref_vcf}.tbi"

//     //             tuples << tuple(
//     //                 chr,
//     //                 file(phased_vcf),
//     //                 file(ref_vcf),
//     //                 file(ref_vcf_index),
//     //                 file(map_file),
//     //                 file(sample_map)
//     //             )
//     //         }
//     //     }
//     //     return tuples
//     // }

//     // 
//     // phased_vcf_with_chr_ch
//     //     .combine(map_file_ch)
//     //     .combine(sample_file_ch)
//     //     .map { item, sample_map ->
//     //         def (chr, phased_vcf, map_file) = item
//     //         def ref_vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
//     //         def ref_vcf_index = "${ref_vcf}.tbi"
//     //         tuple(
//     //             chr,
//     //             file(phased_vcf),
//     //             file(ref_vcf),
//     //             file(ref_vcf_index),
//     //             file(map_file),
//     //             file(sample_map)
//     //         )
//     //     }


//     rfmix_results = run_rfmix(rfmix_inputs_ch)

//     emit:
//         rfmix_results
//     }



// // Default workflow for Cirro to run
// workflow { ancestry_pipeline() }


// // workflow run_download_genetic_map {
// //     download_genetic_map()
// // }
