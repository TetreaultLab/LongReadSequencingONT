#!/usr/bin/env python3
"""
mosdepth_report.py

Create a standalone self-contained HTML report that embeds PNGs as base64.
File naming convention expected (by default):
    <sample>_<slot>.png
e.g. FR-0003_coverage_plot.png  FR-0003_coverage_summary.png

Usage:
    python3 mosdepth_report.py --indir ./qc --out report.html
Options:
    --indir DIR        directory containing PNGs (default: current dir)
    --out PATH         output html path (default ./report.html)
    --delimiter STR    delimiter between sample and slot in filename (default "_")
    --ext EXT          image extension to include (default "png")
"""
import os
import argparse
import glob 
import base64
import re
import html

# Script parameters
parser = argparse.ArgumentParser()
parser.add_argument("--indir", "-i", default=".", help="directory with PNG files (default: current dir)")
parser.add_argument("--out", "-o", default="report.html", help="output html file")
parser.add_argument("--delimiter", "-d", default="_", help="delimiter between sample and slot (default '_')")
parser.add_argument("--ext", default="png", help="image extension to look for (default png)")
args = parser.parse_args()

indir = os.path.abspath(args.indir)
delim = args.delimiter
ext = args.ext.lstrip(".")

# If user doesn't give a report name with --out (default is "report.html") but gives standard input directory:
# The parent folder should be a PoolName (i.e. ./qc -> ../../Pool_FXN3-8-20_MT02/qc).
# Write the report into the same name so the report is self-contained.
if args.out == "report.html" or not args.out:
    pool_name = os.path.basename(os.path.dirname(indir)) or "report"
    outpath = os.path.join(indir, f"{pool_name}_mosdepth_report.html")
else:
    outpath = os.path.abspath(args.out)

# Find files - Assumes naming from previous scripts
pattern = os.path.join(indir, f"*.{ext}")
files = sorted(glob.glob(pattern))
if not files:
    print(f"No .{ext} files found in {indir}. Nothing to do.")
    raise SystemExit(1)

# Parse files into structure: slots -> sample -> filepath
slots = {}   # { slot_name: { sample: filepath } }
samples_set = set()

for fp in files:
    base = os.path.basename(fp)
    name, _ = os.path.splitext(base)
    if delim in name:
        sample, slot = name.split(delim, 1)
    else:
        # Fallback: whole name is sample, slot "image"
        sample, slot = name, "image"
    samples_set.add(sample)
    # Normalize slot label: keep raw slot but also a nice display name
    slot_key = slot  # Used as internal id
    slots.setdefault(slot_key, {})[sample] = fp

samples = sorted(samples_set)

# Helper functions
def safe_id(text):
    # Make a string safe for HTML id/class (letters, digits, underscore, dash)
    return re.sub(r'[^0-9A-Za-z_-]', '_', text)

def nice_label(slot):
    # Turn "coverage_plot" -> "Coverage Plot"
    s = slot.replace('_', ' ').replace('-', ' ')
    return s.title()

def embed_b64(path):
    with open(path, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode("ascii")

# Build HTML
html_parts = []
html_parts.append("<!doctype html>")
html_parts.append("<html><head><meta charset='utf-8'>")
html_parts.append("<title>Sequencing QC Report</title>")
html_parts.append("<style>")
html_parts.append("""
body{font-family: Arial, Helvetica, sans-serif; margin:20px;}
h1,h2{margin:8px 0;}
.tabs {margin:6px 0 12px 0;}
.tab {display:inline-block; padding:6px 10px; margin:2px; border:1px solid #bbb; cursor:pointer; border-radius:6px; background:#f7f7f7;}
.tab.active {background:#dfe9f5; border-color:#6f9ed8;}
.panel {display:none; text-align:center; margin-bottom:20px;}
.panel img {max-width:95%; height:auto; border:1px solid #ddd; box-shadow: 2px 2px 6px rgba(0,0,0,0.08);}
.slot-section {margin-bottom: 34px; padding-bottom:8px; border-bottom:1px solid #eee;}
.meta {font-size:0.9em;color:#666;margin-top:6px;}
.note {font-size:0.9em;color:#444;margin-bottom:8px;}
""")
html_parts.append("</style>")

# JS to show/hide panels within a slot
html_parts.append("<script>")
html_parts.append("""
function show(slot, sample) {
  // hide all panels for this slot
  var panels = document.getElementsByClassName('panel_slot_' + slot);
  for (var i=0; i<panels.length; i++) panels[i].style.display = 'none';
  // show selected
  var id = 'panel_' + slot + '_' + sample;
  var el = document.getElementById(id);
  if (el) el.style.display = 'block';
  // update tab active states
  var tabs = document.getElementsByClassName('tab_slot_' + slot);
  for (var i=0; i<tabs.length; i++) tabs[i].classList.remove('active');
  var tid = 'tab_' + slot + '_' + sample;
  var t = document.getElementById(tid);
  if (t) t.classList.add('active');
}
</script>
""")
html_parts.append("</head><body>")
html_parts.append("<h1>Sequencing QC Report</h1>")
html_parts.append(f"<div class='note'>Generated from images in <strong>{html.escape(indir)}</strong>.<br>"
                  "This HTML was generated only for two plots but more could be added.</div>")

# For each slot, emit tabs and panels
first_slot = True
for slot_key in sorted(slots.keys()):
    slot_safe = safe_id(slot_key)
    slot_label = nice_label(slot_key)
    html_parts.append(f"<div class='slot-section'>")
    html_parts.append(f"<h2>{html.escape(slot_label)}</h2>")

    # Tabs (only for samples that have this slot)
    html_parts.append(f"<div class='tabs' id='tabs_{slot_safe}'>")
    slot_samples = sorted(slots[slot_key].keys())
    if not slot_samples:
        html_parts.append("<em>No images found for this slot.</em>")
    else:
        for i, s in enumerate(slot_samples):
            s_safe = safe_id(s)
            active = " active" if i == 0 else ""
            html_parts.append(f"<span id='tab_{slot_safe}_{s_safe}' class='tab tab_slot_{slot_safe}{active}' "
                              f" onclick=\"show('{slot_safe}','{s_safe}')\">{html.escape(s)}</span>")
    html_parts.append("</div>")  # end tabs

    # Panels with embedded images
    for i, s in enumerate(slot_samples):
        s_safe = safe_id(s)
        panel_id = f"panel_{slot_safe}_{s_safe}"
        style = "display:block" if i == 0 else "display:none"
        html_parts.append(f"<div id='{panel_id}' class='panel panel_slot_{slot_safe}' style='{style}'>")
        img_path = slots[slot_key].get(s)
        if img_path and os.path.exists(img_path):
            b64 = embed_b64(img_path)
            html_parts.append(f"<div class='meta'>File: {html.escape(os.path.basename(img_path))} &nbsp; | &nbsp; Size: {os.path.getsize(img_path):,} bytes</div>")
            html_parts.append(f"<img src='data:image/{ext};base64,{b64}' alt='{html.escape(s)} - {html.escape(slot_label)}'/>")
        else:
            html_parts.append("<div><em>Image file missing</em></div>")
        html_parts.append("</div>")  # end panel

    html_parts.append("</div>")  # end slot-section

html_parts.append("<div style='font-size:0.9em;color:#666;margin-top:18px'>Report generated with mosdepth_report.py</div>")
html_parts.append("</body></html>")

with open(outpath, "w") as fh:
    fh.write("\n".join(html_parts))

print(f"Wrote standalone report to: {outpath}  (contains {len(files)} images embedded)")
