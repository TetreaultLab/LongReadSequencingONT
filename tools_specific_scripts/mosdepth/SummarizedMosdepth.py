#!/usr/bin/env python3
"""
SummaryTable.py (region-aware, robust) + integrated coverage plot

Usage:
  python3 SummaryTable.py --prefix test --thresholds 1,5,10,20,30
  Optional plotting flags:
    --bins 500            (bin size used when mosdepth created regions; only used for plot label)
    --plot-out <path>     (output PNG for coverage plot; default: <prefix>_coverage_plot.png)

Produces:
  - <prefix>_coverage_summary.csv
  - <prefix>_coverage_summary.png  (table image)
  - <prefix>_coverage_plot.png     (coverage plot)  <-- only if <prefix>.regions.bed.gz exists
"""
import argparse
import pandas as pd
import numpy as np
import sys
import re
import matplotlib.pyplot as plt

## Define functions
# Account for different reference format
def is_canonical(chrom):
    ch = chrom
    if isinstance(ch, str) and ch.startswith("chr"):
        ch = ch[3:]
    if ch == "MT":
        ch = "M"
    return ch in set([str(i) for i in range(1,23)] + ["X","Y","M"])

# Define which chromosomes are the main interest
def normalize_chrom_list(use_chr_prefix):
    base = [str(i) for i in range(1,23)] + ["X","Y","M"]
    return (["chr" + c for c in base] if use_chr_prefix else base)

# Correctly parse the thresholds (helper for sorting)
def col_sort_key(colname):
    m = re.search(r'(\d+)', colname)
    return int(m.group(1)) if m else 9999

## Command line definition
parser = argparse.ArgumentParser()
parser.add_argument("--prefix", "-p", required=True, help="mosdepth prefix (e.g. 'Sample1')")
parser.add_argument("--thresholds", "-T", default="1,5,10,20,30",
                    help="comma-separated thresholds in same order it was given to mosdepth (default 1,5,10,20,30)")
parser.add_argument("--out", "-o", default=None, help="output CSV (default: <prefix>_coverage_summary.csv)")
parser.add_argument("--patch-prefixes", default="KI,GL",
                    help="comma-separated contig prefixes to merge as 'patches' (default: KI,GL)")
# Plotting options 
parser.add_argument("--bins", type=int, default=50000, help="bin size used when producing regions (for plot title)")
parser.add_argument("--plot-out", default=None, help="path for coverage plot PNG (default: <prefix>_coverage_plot.png)")
args = parser.parse_args()

## Variable definition (expected from mosdepth)
prefix = args.prefix
thr_list_user = [int(x) for x in args.thresholds.split(",")] if args.thresholds else []
summary_file = f"{prefix}.mosdepth.summary.txt"
thr_file = f"{prefix}.thresholds.bed.gz"
patch_prefixes = tuple(p.strip() for p in args.patch_prefixes.split(","))

# --- Read summary ---
# Fail if the summary is inexistant or has the wrong name
try:
    summary_raw = pd.read_csv(summary_file, sep="\t", header=0, dtype={"chrom":str}, low_memory=False)
except FileNotFoundError:
    print(f"ERROR: cannot find {summary_file}", file=sys.stderr)
    sys.exit(1)

# Drop duplicated "_region" rows mosdepth emits (in order to not count twice)
summary_raw = summary_raw[~summary_raw['chrom'].str.endswith("_region")].copy()
summary = summary_raw.set_index("chrom").copy()
summary['length'] = pd.to_numeric(summary['length'], errors='coerce').fillna(0).astype(np.int64)
summary['bases']  = pd.to_numeric(summary['bases'],  errors='coerce').fillna(0).astype(np.int64)
summary['mean']   = pd.to_numeric(summary['mean'],   errors='coerce').fillna(0.0).astype(float)

# --- Read thresholds file robustly (skip header lines that start with '#') ---
# Fail if the file is inexistant or has the wrong name
try:
    thr = pd.read_csv(thr_file, sep="\t", header=None, comment="#",
                      compression='infer', low_memory=False)
except FileNotFoundError:
    print(f"ERROR: cannot find {thr_file}", file=sys.stderr)
    sys.exit(1)

# Validate the file is in the expected format (usually depend on number of thresholds)
ncols = thr.shape[1]
if ncols < 4:
    raise SystemExit(f"Unexpected thresholds file format (need >=4 cols): {thr_file}")

# Detect whether 4th column is a non-numeric 'region' column (otherwise it shifts everything)
has_region = False
if ncols >= 4:
    # attempt to coerce the 4th column to numeric; if many values fail, treat as region
    col4_coerced = pd.to_numeric(thr.iloc[:, 3], errors='coerce')
    # if more than 50% are non-numeric -> treat as region column
    if col4_coerced.isna().sum() > (len(thr) * 0.5):
        has_region = True

threshold_count = ncols - (4 if has_region else 3)

# Build column names in output
if len(thr_list_user) == threshold_count:
    # If the user inputted thresholds match the threshold files
    if has_region:
        thr_names = ["chrom","start","end","region"] + [f"count_{v}X" for v in thr_list_user]
    else:
        thr_names = ["chrom","start","end"] + [f"count_{v}X" for v in thr_list_user]
else:
    # Otherwise just use values from file but warn user it doesn't match
    if has_region:
        thr_names = ["chrom","start","end","region"] + [f"count_col{i+1}" for i in range(threshold_count)]
    else:
        thr_names = ["chrom","start","end"] + [f"count_col{i+1}" for i in range(threshold_count)]
    print(f"WARNING: provided --thresholds length ({len(thr_list_user)}) does not match "
          f"thresholds file columns ({threshold_count}). Using file-order names: {thr_names[3:]}", file=sys.stderr)

thr.columns = thr_names

# This just takes care of the numeric columns, keep 'region' as string if present
for c in thr.columns:
    if c == "chrom":
        thr[c] = thr[c].astype(str)
    elif c in ("start","end"):
        thr[c] = pd.to_numeric(thr[c], errors='coerce').fillna(0).astype(np.int64)
    elif c == "region":
        thr[c] = thr[c].astype(str)
    else:
        thr[c] = pd.to_numeric(thr[c], errors='coerce').fillna(0).astype(np.int64)

# Sum thresholds by chrom (numeric_only to ignore 'region')
drop_cols = [c for c in ("start","end","region") if c in thr.columns]
thr_agg = thr.groupby("chrom").sum(numeric_only=True).drop(columns=drop_cols, errors='ignore')

# --- Join summary and thresholds ---
# From left to right addition
merged = summary.join(thr_agg, how='left').fillna(0)

# Identify threshold cols present (from previously defined formula)
thr_cols = [c for c in merged.columns if isinstance(c, str) and (c.startswith("count_"))]
thr_cols = sorted(thr_cols, key=col_sort_key)  # sort by number inside column name if present

# Compute proportions (bases >= X) / contig_length
length = merged['length'].replace({0: np.nan})
for col in thr_cols:
    x = re.search(r'(\d+)X', col).group(1)  # extract "1", "5", etc.
    merged[f"Proportion {x}X"] = (merged[col] / length).fillna(0.0)

# Determine chr prefix usage
sample_chroms = list(merged.index)
use_chr_prefix = any(isinstance(x, str) and x.startswith("chr") for x in sample_chroms if isinstance(x, str))
canonical_order = normalize_chrom_list(use_chr_prefix)

# Build main rows (canonical chromosomes)
# Note: is_canonical handles both 'chr1' and '1' forms
main_present = [c for c in canonical_order if c in merged.index]
# If some canonical names are present but not in canonical_order (unexpected), add them
extra_main = [c for c in merged.index if is_canonical(c) and c not in main_present]
main_present.extend(sorted(extra_main))

# Build patch rows: ONLY those starting with the specified patch prefixes:
# Explicitly exclude any 'total' or '*_region' entries (part of normal output)
patch_rows = [c for c in merged.index
              if isinstance(c, str)
              and c not in ("total",)
              and not c.endswith("_region")
              and (not is_canonical(c))
              and c.startswith(patch_prefixes)]

# --- Assemble output rows in the desired order ---
# One row per canonical & patch + total
out_rows = []
for c in main_present:
    row = merged.loc[c]
    out_rows.append((c, row))

# Aggregate patches if any
if len(patch_rows) > 0:
    p_len = int(merged.loc[patch_rows, "length"].sum())
    p_bases = int(merged.loc[patch_rows, "bases"].sum()) if "bases" in merged.columns else 0
    p_mean = float(p_bases / p_len) if p_len > 0 else 0.0
    patch_dict = {"length": p_len, "bases": p_bases, "mean": p_mean}
    for col in thr_cols:  # sum raw counts
        col_sum = int(merged.loc[patch_rows, col].sum())
        patch_dict[col] = col_sum
        x = re.search(r'(\d+)X', col).group(1)
        patch_dict[f"Proportion {x}X"] = (col_sum / p_len) if p_len > 0 else 0.0
    out_rows.append(("Patches (Combined)", pd.Series(patch_dict)))

# Ensure we're not double-counting 'total' or 'total_region'
clean_rows = merged.loc[~merged.index.str.lower().str.contains("total")].copy()

# First: compute per-contig proportions (if not already present)
for col in thr_cols:
    x = re.search(r'(\d+)X', col).group(1)
    prop_col = f"Proportion {x}X"
    clean_rows[prop_col] = clean_rows[col] / clean_rows['length']

# Now compute total stats
total_len = int(clean_rows['length'].sum())
total_bases = int(clean_rows['bases'].sum()) if "bases" in clean_rows.columns else 0
total_mean = float(total_bases / total_len) if total_len > 0 else 0.0

total_dict = {
    "length": total_len,
    "bases": total_bases,
    "mean": total_mean,
}

# Include raw threshold sums (not really needed for weighted prop, but OK to keep)
for col in thr_cols:
    col_sum = int(clean_rows[col].sum())
    total_dict[col] = col_sum
    x = re.search(r'(\d+)X', col).group(1)
    prop_col = f"Proportion {x}X"
    weighted_prop_sum = (clean_rows[prop_col] * clean_rows['length']).sum()
    total_dict[prop_col] = weighted_prop_sum / total_len if total_len > 0 else 0.0

# Add total row
out_rows.append(("Total", pd.Series(total_dict)))

# Total row
#total_len = int(merged['length'].sum())
#total_bases = int(merged['bases'].sum()) if "bases" in merged.columns else 0
#total_mean = float(total_bases / total_len) if total_len > 0 else 0.0
#
#total_dict = {"length": total_len, "bases": total_bases, "mean": total_mean}
#
## Weighted average for proportions
#for col in thr_cols:
#    x = re.search(r'(\d+)X', col).group(1)
#    prop_col = f"Proportion {x}X"
#    # Weighted by length
#    weighted_num = (merged[prop_col] * merged['length']).sum()
#    total_dict[prop_col] = (weighted_num / total_len) if total_len > 0 else 0.0
#
#out_rows.append(("Total", pd.Series(total_dict)))

# --- Build DataFrame, round values to 3 decimals ---
# At more than 3 decimals it's messy, but also the average_coverage stops at 2 decimals
rows = []
for name, s in out_rows:
    d = {
        "Chromomosome": name,
        "Length (bp)": int(s.get("length", 0)),
        "Avg Coverage (X)": round(float(s.get("mean", 0.0)), 3)
    }
    for col in thr_cols:
        x = re.search(r'(\d+)X', col).group(1)
        prop_col = f"Proportion {x}X"
        d[prop_col] = round(float(s.get(prop_col, 0.0)), 3)
    rows.append(d)
  

# Output file format
df_out = pd.DataFrame(rows)
out_file = args.out if args.out else f"{prefix}_coverage_summary.csv"
df_out.to_csv(out_file, index=False)
print(f"Wrote summary to {out_file}")

# ---- Make visual output for the user (table image) ----
def save_table_as_png(df, output_file):
    # Set up the plot and axes
    fig, ax = plt.subplots(figsize=(12, 6))  # Adjust size as needed
    ax.axis('tight')
    ax.axis('off')

    # Create the table
    table = ax.table(cellText=df.values,
                    colLabels=df.columns,
                    cellLoc='center',
                    loc='center')

    table.auto_set_font_size(False)
    table.set_fontsize(12)

    for (i, j), cell in table.get_celld().items():
        if i == 0 or j == 0:  # Header row or header column
            cell.set_fontsize(14)
            cell.set_text_props(color='white', weight='bold')
            cell.set_facecolor('#285780')
        else:
            cell.set_fontsize(10)
            cell.set_text_props(color='black')
            cell.set_facecolor('#f9f9f9')
            if i % 2 == 0:
                cell.set_facecolor('#e9e9e9')
        cell.set_edgecolor('black')

    table.scale(1.2, 1.5)
    table.auto_set_column_width(col=list(range(len(df.columns))))

    png_file = output_file.replace('.csv', '.png').replace('.tsv', '.png')
    plt.savefig(png_file, format='png', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Saved table image to {png_file}")

# Call the function to save the PNG of the table
save_table_as_png(df_out, out_file)

# ----- Read regions file and then plot: scale based on Total average and highlight outliers -----
regions_file = f"{prefix}.regions.bed.gz"
try:
    regs = pd.read_csv(regions_file, sep="\t", header=None, compression='infer',
                       names=["chrom","start","end","mean_cov"], low_memory=False)
except FileNotFoundError:
    print(f"Regions file '{regions_file}' not found — skipping coverage plot.", file=sys.stderr)
    # exit 0 so main CSV+table PNG are considered successful
    sys.exit(0)

# Drop any _region duplicate rows if present
regs = regs[~regs['chrom'].str.endswith("_region")].copy()

# Ensure numeric types
regs['start'] = pd.to_numeric(regs['start'], errors='coerce').fillna(0).astype(np.int64)
regs['end'] = pd.to_numeric(regs['end'], errors='coerce').fillna(0).astype(np.int64)
regs['mean_cov'] = pd.to_numeric(regs['mean_cov'], errors='coerce').fillna(0.0).astype(float)

# Detect chr prefix usage in regions file
sample_chroms = set(regs['chrom'].unique())
use_chr_prefix_regs = any(str(x).startswith("chr") for x in sample_chroms)

# Build canonical order for those present (uses your normalize_chrom_list)
canon = normalize_chrom_list(use_chr_prefix_regs)
present_canon = [c for c in canon if c in sample_chroms]

if not present_canon:
    print("No canonical chromosomes found in regions file. Skipping coverage plot.", file=sys.stderr)
    sys.exit(0)

# Compute chromosome sizes from regions (max end)
chrom_sizes = regs.groupby('chrom')['end'].max().to_dict()

# Build offsets cumulatively (linear coordinate)
offset = {}
cum = 0
for c in present_canon:
    offset[c] = cum
    cum += int(chrom_sizes.get(c, 0))

# Filter to canonical present and compute linear coordinate (fast vectorized mapping)
regs = regs[regs['chrom'].isin(present_canon)].copy()
regs['coord'] = regs['start'].astype(int) + regs['chrom'].map(offset).astype(int)

# Sort by linear coordinate to avoid diagonal connectors between chromosomes
regs.sort_values('coord', inplace=True)

# ----- Determine y-limits from "Total" average in the summary table -----
if "Chromomosome" in df_out.columns and "Total" in df_out["Chromomosome"].values:
    avg_total = float(df_out.loc[df_out["Chromomosome"] == "Total", "Avg Coverage (X)"].iloc[0])
else:
    # fallback: use merged row 'total' or genome mean
    avg_total = float(merged.loc["total", "mean"]) if "total" in merged.index else float(merged["mean"].mean())

# Y limits rule: ymin = max(0, avg_total - 30), ymax = avg_total + 30
ymin = max(0.0, avg_total - 30.0)
ymax = avg_total + 30.0

# Create clipped display values so axis is not dominated by spikes.
display_vals = regs['mean_cov'].clip(lower=ymin, upper=ymax)

# Color points by their true value relative to the y-limits (green = above, red = below, blue = within)
colors = np.where(regs['mean_cov'] > ymax, 'green',
         np.where(regs['mean_cov'] < ymin, 'red', 'tab:blue'))

# Plot: clipped line + colored points for outliers
fig, ax = plt.subplots(figsize=(14,4))
ax.plot(regs['coord'], display_vals, lw=0.4, color='tab:blue', zorder=1)
ax.scatter(regs['coord'], display_vals, c=colors, s=2, linewidths=0, zorder=2)

# Horizontal reference lines: average (dashed), ymin/ymax (dotted)
ax.axhline(avg_total, color='gray', linestyle='--', linewidth=1, label=f'avg={avg_total:.2f}')
ax.axhline(ymin, color='red', linestyle=':', linewidth=0.8)
ax.axhline(ymax, color='green', linestyle=':', linewidth=0.8)

# Set Y limits explicitly
ax.set_ylim(ymin, ymax)

# Labels + title
ax.set_xlabel("Genome (concatenated chromosomes)")
ax.set_ylabel(f"Mean depth per {args.bins}bp bin")
ax.set_title(f"{prefix} coverage ({args.bins}bp bins)")

# Chromosome tick positions at midpoint (same logic as before)
ticks = []
labels = []
for c in present_canon:
    size = chrom_sizes.get(c, 0)
    if size > 0:
        ticks.append(offset[c] + size/2)
        lab = c[3:] if use_chr_prefix_regs else c
        if lab == "MT":
            lab = "M"
        labels.append(lab)
ax.set_xticks(ticks)
ax.set_xticklabels(labels, rotation=90, fontsize=8)

# Legend
ax.legend(loc='upper right', fontsize=8)

plt.tight_layout()
outpath = args.plot_out if args.plot_out else f"{prefix}_coverage_plot.png"
plt.savefig(outpath, dpi=150)
plt.close(fig)
print(f"Wrote coverage plot to {outpath}")
