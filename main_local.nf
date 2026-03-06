#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { phase_with_eagle } from './modules/eagle'
include { run_rfmix } from './modules/rfmix'

if( !params.input_vcf ) {
  log.error "Usage: nextflow run main_local.nf --input_vcf /path/to/input.vcf.gz"
  System.exit(1)
}

def INPUT_VCF = params.input_vcf

def inputFileObj = new File(INPUT_VCF)
if( !inputFileObj.exists() ) {
  log.error "Input VCF not found: ${INPUT_VCF}"
  System.exit(1)
}
/*
 * ----------------------------------------------------------------------------
 * Processes (your existing ones below are unchanged)
 * ----------------------------------------------------------------------------
 */

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

  for i in {1..22}; do
      awk '{
          print \$1, \$4, \$3
      }' OFS=' ' no_chr_in_chrom_field/plink.chr\${i}.GRCh38.map
  done > genetic_map_all_chr.txt

  echo "Genetic map ready (uncompressed for RFMix compatibility)"
  """
}

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

process normalize_chrom_names {
  tag "$vcf_file"

  input:
    path vcf_file

  output:
    tuple path("normalized_${vcf_file.simpleName}.vcf.gz"), path("normalized_${vcf_file.simpleName}.vcf.gz.tbi"), emit: normalized

  script:
  """
  set -euo pipefail

  OUT_VCF="normalized_${vcf_file.simpleName}.vcf.gz"

  echo "Normalizing chromosome names for ${vcf_file}..."

  FIRST_CHROM=\$(zcat ${vcf_file} 2>/dev/null | grep -v '^#' | head -1 | cut -f1 || true)
  echo "First chromosome found: \${FIRST_CHROM}"

  if [[ "\${FIRST_CHROM}" == chr* ]]; then
      echo "Chromosomes start with 'chr', converting to numbers..."

      (zcat ${vcf_file} | grep -v '^#' | cut -f1 | sort -u | grep -v '^chr\$' || true) | while read chrom; do
          new_chrom=\$(echo "\${chrom}" | sed 's/^chr//' | sed 's/^M\$/MT/')
          if [[ -n "\${new_chrom}" ]]; then
              echo -e "\${chrom}\\t\${new_chrom}"
          fi
      done > reheader.txt

      echo "Chromosome mapping:"
      cat reheader.txt

      bcftools view -h ${vcf_file} | \\
      sed 's/^##contig=<ID=chr/##contig=<ID=/' | \\
      sed 's/^##contig=<ID=M,/##contig=<ID=MT,/' > temp_header.vcf

      bcftools annotate --rename-chrs reheader.txt ${vcf_file} 2>/dev/null | \\
      bcftools view -H >> temp_header.vcf

      bgzip -c temp_header.vcf > \${OUT_VCF}
      rm temp_header.vcf
  else
      echo "Chromosomes are already numeric (found: \${FIRST_CHROM}), just copying..."
      bcftools view -Oz -o \${OUT_VCF} ${vcf_file}
  fi

  tabix -p vcf \${OUT_VCF}
  """
}

/*
 * ----------------------------------------------------------------------------
 * Workflow
 * ----------------------------------------------------------------------------
 */
workflow ancestry_pipeline {

  def chr_ch = Channel.from(1..22)

  download_genetic_map()
  download_genetic_map_eagle()
  download_sample_map()

  def map_file_ch        = download_genetic_map.out.map
  def map_file_eagle_ch  = download_genetic_map_eagle.out.map
  def sample_file_ch     = download_sample_map.out

  // Local input VCF
  def input_vcf_ch = Channel.fromPath(INPUT_VCF, checkIfExists: true)

  // Normalize once
  def normalized_vcf_ch = input_vcf_ch | normalize_chrom_names

  // Build Eagle inputs per chromosome
  def eagle_inputs_ch = chr_ch
    .combine(normalized_vcf_ch.normalized)
    .combine(map_file_eagle_ch)
    .map { chr, vcf_file, vcf_index, map_file ->
      def ref_vcf = "s3://1000genomes/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      def ref_vcf_index = "${ref_vcf}.tbi"
      tuple(
        file(vcf_file),
        file(vcf_index),
        file(ref_vcf),
        file(ref_vcf_index),
        file(map_file),
        chr
      )
    }

  /*
   * IMPORTANT: I’m assuming phase_with_eagle emits tuples like:
   *   tuple(chr, phased_vcf)
   * If your module emits just "phased_vcf" without chr, tell me and I’ll adjust.
   */
  def phased_vcf_ch = phase_with_eagle(eagle_inputs_ch)

  // Build RFMix inputs (broadcast map + sample map)
  def rfmix_inputs_ch = phased_vcf_ch
    .combine(map_file_ch)
    .combine(sample_file_ch)
    .map { phased_tuple, map_file, sample_map ->
      def chr        = phased_tuple[0]
      def phased_vcf = phased_tuple[1]

      def ref_vcf = "s3://1000genomes/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      def ref_vcf_index = "${ref_vcf}.tbi"

      tuple(
        chr,
        file(phased_vcf),
        file(ref_vcf),
        file(ref_vcf_index),
        file(map_file),
        file(sample_map)
      )
    }

    run_rfmix(rfmix_inputs_ch)

    emit:
        rfmix_results = run_rfmix.out
}

workflow { ancestry_pipeline() }