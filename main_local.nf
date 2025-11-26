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
        curl -s -L -o genetic_map_hg38_withX.txt.gz \\
            https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
    fi

    echo "Decompressing..."
    gunzip -c genetic_map_hg38_withX.txt.gz > genetic_map_raw.txt

    # -------------------------------------------------------
    # 2. Keep only chr, pos, mapCentiMorgan (column 4)
    # -------------------------------------------------------
    echo "Extracting 3 required columns..."
    python3 - << 'EOF'
INPUT  = "genetic_map_raw.txt"
OUTPUT = "genetic_map_chr.txt"

with open(INPUT) as fin, open(OUTPUT, "w") as fout:
    for line in fin:
        parts = line.strip().split()
        if not parts:
            continue

        # Skip Eagle header line
        if parts[0].lower() in ["chr", "chromosome"] or parts[1].lower() in ["pos", "position"]:
            continue

        if re.fullmatch(r"[0-9]+", chrom):
            num = int(chrom)
            if 1 <= num <= 22:
                parts[0] = f"chr{num}"
            elif num == 23:
                parts[0] = "chrX"
            elif num == 24:
                parts[0] = "chrY"
        chr_, pos, rate, mapcM = parts[:4]
        fout.write(f"{chr_} {pos} {mapcM}\\n")
EOF

    # -------------------------------------------------------
    # 3. Clean + enforce strict monotonic genetic map per row
    # -------------------------------------------------------
    echo "Cleaning map (strictly increasing cM)..."
    python3 - << 'EOF'
import pandas as pd
import numpy as np

INPUT = "genetic_map_chr.txt"
OUTPUT = "genetic_map_chr_cleaned.txt"

epsilon = 1e-6

print(f"Reading {INPUT} ...")

# No header in file → header=None, define names manually
df = pd.read_csv(INPUT, sep="\\s+", header=None, names=["chr", "pos", "map"])

# Enforce monotonic map (required for RFMix)
df["map"] = np.maximum.accumulate(df["map"].values + np.arange(len(df)) * epsilon)

# Save whitespace-separated, no header
df.to_csv(OUTPUT, sep=" ", index=False, header=False, float_format="%.12f")

print(f"Saved cleaned map to {OUTPUT}")
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

// workflow download_only {
//     download_genetic_map()
// }

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

    rfmix_results = run_rfmix(rfmix_inputs_ch)

    emit:
        rfmix_results
    }



// Default workflow for Cirro to run
workflow { ancestry_pipeline() }




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

