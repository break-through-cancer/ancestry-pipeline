#!/usr/bin/env python3

INPUT = "genetic_map_hg38_withX.txt"     # plain text input
OUTPUT = "genetic_map_chr_test.txt"      # plain text output

# Valid chromosome names for autosomes
VALID_AUTOSOMES = {str(i) for i in range(1, 23)}
SPECIAL_MAP = {
    "23": "X",
    "24": "Y"
}

def fix_line(line):
    """
    Adds 'chr' prefix to the chromosome column.
    Converts '23' → 'chrX' because RFMix expects X.
    Converts '24' → 'chrY' if present.
    """
    parts = line.strip().split()
    if not parts:
        return line

    chrom = parts[0]

    # If header, keep unchanged
    if chrom.lower() in ["chromosome", "chr", "pos", "position"]:
        return line

    # Convert special chromosomes 23→X, 24→Y
    if chrom in SPECIAL_MAP:
        chrom_new = SPECIAL_MAP[chrom]
        parts[0] = "chr" + chrom_new
        return "\t".join(parts) + "\n"

    # Convert autosomes 1..22
    if chrom in VALID_AUTOSOMES:
        parts[0] = "chr" + chrom
        return "\t".join(parts) + "\n"

    # Otherwise leave unchanged
    return line


def main():
    print(f"Reading: {INPUT}")
    print(f"Writing cleaned file to: {OUTPUT}")

    with open(INPUT, "r") as fin, open(OUTPUT, "w") as fout:
        for line in fin:
            fout.write(fix_line(line))

    print("Done! Inspect the output with:")
    print(f"  head {OUTPUT}")


if __name__ == "__main__":
    main()
