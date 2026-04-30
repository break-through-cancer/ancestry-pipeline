import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re

files = sorted(
    glob.glob("cirro_p2_run/run/*.msp.tsv"),
    key=lambda x: int(os.path.basename(x).split(".")[0])
)

cols = ["chrom", "start_bp", "end_bp", "start_cm", "end_cm", "num_snps", "hap1", "hap2"]

# =========================
# Extract RFMix ancestry mapping
# =========================
code_maps = {}

for f in files:
    with open(f) as fh:
        for line in fh:
            if line.startswith("#Subpopulation order/codes"):
                code_maps[f] = line.strip()
                break

unique_maps = set(code_maps.values())

if len(unique_maps) != 1:
    print("Different ancestry mappings found!")
    for f, m in code_maps.items():
        print(f, "→", m)
    raise ValueError("Inconsistent ancestry mappings across files")

mapping_line = list(unique_maps)[0]
pairs = re.findall(r"(\w+)=(\d+)", mapping_line)
code_to_pop = {int(idx): pop for pop, idx in pairs}

print("Population code mapping:")
print(code_to_pop)

# =========================
# Collapse populations into groups
# =========================
pop_to_group = {
    "ACB": "African", "ASW": "African", "ESN": "African", "GWD": "African",
    "LWK": "African", "MSL": "African", "YRI": "African",

    "CEU": "European", "FIN": "European", "GBR": "European",
    "IBS": "European", "TSI": "European",

    "CDX": "East Asian", "CHB": "East Asian", "CHS": "East Asian",
    "JPT": "East Asian", "KHV": "East Asian",

    "BEB": "South Asian", "GIH": "South Asian", "ITU": "South Asian",
    "PJL": "South Asian", "STU": "South Asian",

    "CLM": "Admixed American", "MXL": "Admixed American",
    "PEL": "Admixed American", "PUR": "Admixed American",
}

code_to_group = {
    code: pop_to_group.get(pop, "Other")
    for code, pop in code_to_pop.items()
}

group_colors = {
    "African": "tab:blue",
    "European": "tab:orange",
    "East Asian": "tab:green",
    "South Asian": "tab:red",
    "Admixed American": "tab:purple",
    "Other": "tab:gray",
}

# =========================
# Genome-wide local ancestry plot
# =========================
fig, ax = plt.subplots(figsize=(16, 10))

y_offset = 0
yticks = []
yticklabels = []
groups_used = set()
max_end_bp = 0

for f in files:
    chrom = int(os.path.basename(f).split(".")[0])
    df = pd.read_csv(f, sep="\t", comment="#", header=None, names=cols)

    max_end_bp = max(max_end_bp, df["end_bp"].max())
    first_start = df["start_bp"].min()

    for _, row in df.iterrows():
        for hap_i, hap in enumerate(["hap1", "hap2"]):
            anc_code = int(row[hap])
            group = code_to_group.get(anc_code, "Other")
            groups_used.add(group)

            xmin = 0 if row["start_bp"] == first_start else row["start_bp"]

            ax.hlines(
                y=y_offset + hap_i,
                xmin=xmin,
                xmax=row["end_bp"],
                color=group_colors[group],
                linewidth=6,
            )

    yticks.append(y_offset + 0.5)
    yticklabels.append(str(chrom))
    y_offset += 3

ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels, fontsize=8)
ax.set_ylabel("Chromosome")
ax.set_xlabel("Genomic position (bp)")
ax.set_title("Genome-wide local ancestry painting", pad=12)
ax.set_xlim(0, max_end_bp)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

handles = [
    plt.Line2D([0], [0], color=group_colors[g], lw=6, label=g)
    for g in ["African", "European", "East Asian", "South Asian", "Admixed American", "Other"]
    if g in groups_used
]

# ax.legend(
#     handles=handles,
#     title="Ancestry group",
#     loc="upper right",
#     fontsize=8,
#     title_fontsize=9,
#     frameon=True,
# )

plt.tight_layout()
plt.savefig("genome_local_ancestry_grouped.png", dpi=300, bbox_inches="tight")
plt.show()

# =========================
# Pie chart: genome-wide ancestry proportions
# =========================
group_lengths = {}

for f in files:
    df = pd.read_csv(f, sep="\t", comment="#", header=None, names=cols)

    for _, row in df.iterrows():
        segment_len = row["end_bp"] - row["start_bp"]

        for hap in ["hap1", "hap2"]:
            anc_code = int(row[hap])
            group = code_to_group.get(anc_code, "Other")
            group_lengths[group] = group_lengths.get(group, 0) + segment_len

total_len = sum(group_lengths.values())

group_percents = {
    group: 100 * length / total_len
    for group, length in group_lengths.items()
}

print("\nGenome-wide ancestry proportions:")
for group, pct in sorted(group_percents.items(), key=lambda x: x[1], reverse=True):
    print(f"{group}: {pct:.2f}%")

sorted_groups = sorted(group_percents.items(), key=lambda x: x[1], reverse=True)

pie_labels = [f"{group}\n{pct:.1f}%" for group, pct in sorted_groups]
pie_sizes = [pct for group, pct in sorted_groups]
pie_colors = [group_colors[group] for group, pct in sorted_groups]

plt.figure(figsize=(7, 7))
plt.pie(
    pie_sizes,
    labels=pie_labels,
    colors=pie_colors,
    startangle=90,
    counterclock=False,
)
plt.title("Genome-wide ancestry proportions by ancestry group")
plt.tight_layout()
plt.savefig("genome_ancestry_continent_pie.png", dpi=300, bbox_inches="tight")
plt.show()