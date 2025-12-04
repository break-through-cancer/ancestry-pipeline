#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { phase_with_eagle } from './modules/eagle'
include { run_rfmix } from './modules/rfmix'
include { phase_with_eagle as phase_with_eagle_ref } from './modules/eagle'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input_genotype) { input_genotype = params.input_genotype } else { exit 1, 'Please, provide an input genotype data !' }
if (params.input_genotype_index) { input_genotype_input = params.input_genotype_index } else { exit 1, 'Please, provide an input genotype index !' }
if (params.reference_vcf) { reference_vcf = params.reference_vcf } else { exit 1, 'Please, provide a reference vcf!' }
if (params.reference_vcf_index) { reference_vcf_index = params.reference_vcf_index } else { exit 1, 'Please, provide a reference vcf ref!' }

// process download_genetic_map {

//     output:
//         path "genetic_map_chr_cleaned.txt.gz"

//     script:
//     """
//     echo "=== Starting download_genetic_map process ==="

//     # -------------------------------------------------------
//     # 1. Download original genetic map if missing
//     # -------------------------------------------------------
//     if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
//         echo "Downloading hg38 genetic map..."
//         curl -s -L -o genetic_map_hg38_withX.txt.gz \\
//             https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
//     fi

//     echo "Decompressing..."
//     gunzip -c genetic_map_hg38_withX.txt.gz > genetic_map_raw.txt

//     # -------------------------------------------------------
//     # 2. Keep only chr, pos, mapCentiMorgan (column 4)
//     # -------------------------------------------------------
//     echo "Extracting 3 required columns..."
//     python3 - << 'EOF'
// INPUT  = "genetic_map_raw.txt"
// OUTPUT = "genetic_map_chr.txt"

// with open(INPUT) as fin, open(OUTPUT, "w") as fout:
//     for line in fin:
//         parts = line.strip().split()
//         if not parts:
//             continue

//         # Skip Eagle header line
//         if parts[0].lower() in ["chr", "chromosome"] or parts[1].lower() in ["pos", "position"]:
//             continue

//         if re.fullmatch(r"[0-9]+", chrom):
//             num = int(chrom)
//             if 1 <= num <= 22:
//                 parts[0] = f"chr{num}"
//             elif num == 23:
//                 parts[0] = "chrX"
//             elif num == 24:
//                 parts[0] = "chrY"
//         chr_, pos, rate, mapcM = parts[:4]
//         fout.write(f"{chr_} {pos} {mapcM}\\n")
// EOF

//     # -------------------------------------------------------
//     # 3. Clean + enforce strict monotonic genetic map per row
//     # -------------------------------------------------------
//     echo "Cleaning map (strictly increasing cM)..."
//     python3 - << 'EOF'
// import pandas as pd
// import numpy as np

// INPUT = "genetic_map_chr.txt"
// OUTPUT = "genetic_map_chr_cleaned.txt"

// epsilon = 1e-6

// print(f"Reading {INPUT} ...")

// # No header in file → header=None, define names manually
// df = pd.read_csv(INPUT, sep="\\s+", header=None, names=["chr", "pos", "map"])

// # Enforce monotonic map (required for RFMix)
// df["map"] = np.maximum.accumulate(df["map"].values + np.arange(len(df)) * epsilon)

// # Save whitespace-separated, no header
// df.to_csv(OUTPUT, sep=" ", index=False, header=False, float_format="%.12f")

// print(f"Saved cleaned map to {OUTPUT}")
// EOF

//     # -------------------------------------------------------
//     # 4. Compress final cleaned map
//     # -------------------------------------------------------
//     echo "Compressing final cleaned map..."
//     gzip -c genetic_map_chr_cleaned.txt > genetic_map_chr_cleaned.txt.gz

//     echo "=== Done! Preview first 10 lines ==="
//     zcat genetic_map_chr_cleaned.txt.gz | head
//     """
// }
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

process download_genetic_map_eagle_plain {

    output:
        path "genetic_map_hg38_withX.txt",emit: map
    script:
    """
    echo "Downloading Eagle GRCh38 genetic map..."
    if [ ! -f genetic_map_hg38_withX.txt.gz ]; then
        curl -s -L -o genetic_map_hg38_withX.txt.gz \
            https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
    else
        echo "Genetic map already exists, skipping download."
    fi

    echo "Decompressing to plain text..."
    gunzip -c genetic_map_hg38_withX.txt.gz > genetic_map_hg38_withX.txt
    """
}

// process normalize_chrom_names {

//     tag "$vcf_file"

//     input:
//     path vcf_file

//     output:
//         tuple path("normalized_${vcf_file.simpleName}.vcf.gz"), path("normalized_${vcf_file.simpleName}.vcf.gz.tbi"), emit: normalized

//     script:
//     """
//     echo "Normalizing chromosome names for ${vcf_file}..."

//     # Check if CHROM column contains 'chr'
//     if zgrep -v '^##' ${vcf_file} | cut -f1 | grep -q '^chr'; then
//         echo "Chromosomes start with 'chr', converting to numbers..."

//         # Create temporary reheader file
//         zgrep -v '^##' ${vcf_file} | cut -f1 | sort -u | \\
//         awk '{ gsub("chr",""); if(\$1=="M") print \$1 "\\tMT"; else print \$1 "\\t" \$1 }' > reheader.txt

//         # Apply reheader and compress
//         bcftools annotate --rename-chrs reheader.txt -Oz -o normalized_${vcf_file.baseName}.vcf.gz ${vcf_file}

//     else
//         echo "Chromosomes are already numeric, just compressing..."
//         bgzip -c ${vcf_file} > normalized_${vcf_file.baseName}.vcf.gz
//     fi

//     # Index the output VCF
//     bcftools index normalized_${vcf_file.baseName}.vcf.gz
//     """
// }
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
            tuple(
                file(vcf_file),
                file(vcf_index),
                file(params.reference_vcf),
                file(params.reference_vcf_index),
                file(map_file), 
                chr
            )
        }

    phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

    phased_vcf_with_chr_ch = phased_vcf_ch
        .combine(eagle_inputs_ch.map { it[5] })  // combine with the chromosomes
        .map { phased_vcf_file, chr ->
            println "Eagle finished for chromosome ${chr}: phased VCF = ${phased_vcf_file}"
            tuple(chr, phased_vcf_file)
        }
    rfmix_inputs_ch = phased_vcf_with_chr_ch
        .combine(map_file_ch)        // RFMix map
        .combine(sample_file_ch)     // Sample map
        .combine(normalized_vcf_ch.normalized)  // Normalized reference VCF (tuple: [vcf.gz, .tbi])
        .map { items ->
            // Items is a list containing all combined elements
            // The last element is the tuple [vcf.gz, .tbi] from normalized_vcf_ch
            def chr = items[0]
            def phased_vcf = items[1]
            def map_file = items[2]
            def sample_map = items[3]
           
            println "Preparing RFMix for chr ${chr}"
            
            tuple(
                chr,
                file(phased_vcf),
                file(params.reference_vcf),
                file(params.reference_vcf_index),
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
        path "rfmix_sample_map.txt", emit: map
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


    // // now combine chromosomes with the actual file
    // eagle_inputs_ch = chr_ch.combine(map_file_eagle_ch).map { chr, map_file ->
    //     println "Preparing Eagle inputs for chromosome ${chr} with genetic map ${map_file}"
    //     // def vcf = "s3://1000genomes/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz"
    //     // def vcf_index = "${vcf}.tbi"
    //     [
    //         file(params.input_genotype),
    //         file(params.input_genotype_index),
    //         file(map_file), 
    //         chr
    //     ]
    // }

    // // phase the reference files also
    // reference_inputs_eagle_phasing = chr_ch.combine(map_file_eagle_ch).map { chr, map_file ->
    //     println "Preparing reference inputs for chromosome ${chr} through phasing"
    //     def ref_vcf = "s3://1000genomes/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    //     def ref_vcf_index = "${ref_vcf}.tbi"
    //     [
    //         file(params.reference_vcf),                    // reference VCF
    //         file(params.reference_vcf_index),              // reference VCF index
            
    //         file(map_file), 
    //         chr
    //     ]
    // }

 

    // phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)
    // phased_ref_ch = phase_with_eagle_ref(reference_inputs_eagle_phasing)

    // phased_vcf_with_chr_ch = phased_vcf_ch
    //     .combine(eagle_inputs_ch.map { it[3] })  // combine with the chromosomes
    //     .map { phased_vcf_file, chr ->
    //         println "Eagle finished for chromosome ${chr}: phased VCF = ${phased_vcf_file}"
    //         tuple(chr, phased_vcf_file)
    //     }

    // phased_ref_with_chr_ch = phased_ref_ch
    //     .combine(reference_inputs_eagle_phasing.map { it[3] })  // chromosome and the index
    //     .map { phased_ref_file, chr ->
    //         tuple(chr, phased_ref_file)
    //     }

    // rfmix_inputs_ch = phased_vcf_with_chr_ch
    //     .join(phased_ref_with_chr_ch)
    //     .combine(map_file_ch)
    //     .combine(sample_file_ch)
    //     .map { it.flatten() }   // now each element is [chr, phased_vcf, map_file, sample_file]
    //     .map { combined ->
    //         def chr = combined[0][0]
    //         def phased_input_vcf = combined[0][1]
    //         def phased_ref_vcf = combined[1][1]
    //         def map_file = combined[2]
    //         def sample_map = combined[3]
    //         def ref_vcf = "s3://1000genomes/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    //         def ref_vcf_index = "${ref_vcf}.tbi"


    //         tuple(
    //             chr,
    //             file(phased_input_vcf),
    //             file(phased_ref_vcf),
    //             file(ref_vcf_index),
    //             file(map_file, copy: true),
    //             file(sample_map, copy: true)
    //         )
    //     }
    // now combine chromosomes with the actual file