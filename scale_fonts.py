#!/usr/bin/env python3
"""
Scale up all font sizes in analysis scripts for improved readability.
- Hard-coded fontsize=N → fontsize=int(N * SCALE)
- rcParams font.size → increased
- rcParams axes.titlesize → increased
- DPI → 300 for high-resolution output
"""

import re
import os

SCALE = 1.5  # 50% larger text
DPI_NEW = 300

SCRIPTS = [
    'analysis.py',
    'pathway_enrichment.py',
    'seasonal_maturity_progression.py',
    'biomarker_discovery.py',
    'gwas_qtl_crossref.py',
    'stratified_vs_pooled_comparison.py',
    'environment_titer_integration.py',
    'heat_stress_signatures.py',
    'analysis_schematic.py',
    'grbv_biology_diagram.py',
]

WORK_DIR = '/Users/prempratapsingh/Data_Analysis'


def scale_fontsize(match):
    """Scale a fontsize=N parameter."""
    prefix = match.group(1)  # 'fontsize=' or 'font_size='
    value = match.group(2)
    try:
        num = float(value)
        scaled = max(8, int(round(num * SCALE)))  # minimum 8pt
        return f'{prefix}{scaled}'
    except ValueError:
        return match.group(0)


def scale_pt_size(match):
    """Scale Pt(N) in presentation/pdf scripts."""
    value = match.group(1)
    try:
        num = float(value)
        scaled = max(8, int(round(num * SCALE)))
        return f'Pt({scaled})'
    except ValueError:
        return match.group(0)


def update_rcparams(content):
    """Update rcParams font sizes and DPI."""
    # font.size
    content = re.sub(
        r"'font\.size'\s*:\s*(\d+\.?\d*)",
        lambda m: f"'font.size': {max(12, int(round(float(m.group(1)) * SCALE)))}",
        content
    )
    # axes.titlesize
    content = re.sub(
        r"'axes\.titlesize'\s*:\s*(\d+\.?\d*)",
        lambda m: f"'axes.titlesize': {max(14, int(round(float(m.group(1)) * SCALE)))}",
        content
    )
    # axes.labelsize
    content = re.sub(
        r"'axes\.labelsize'\s*:\s*(\d+\.?\d*)",
        lambda m: f"'axes.labelsize': {max(12, int(round(float(m.group(1)) * SCALE)))}",
        content
    )
    # xtick.labelsize
    content = re.sub(
        r"'xtick\.labelsize'\s*:\s*(\d+\.?\d*)",
        lambda m: f"'xtick.labelsize': {max(10, int(round(float(m.group(1)) * SCALE)))}",
        content
    )
    # ytick.labelsize
    content = re.sub(
        r"'ytick\.labelsize'\s*:\s*(\d+\.?\d*)",
        lambda m: f"'ytick.labelsize': {max(10, int(round(float(m.group(1)) * SCALE)))}",
        content
    )
    # legend.fontsize
    content = re.sub(
        r"'legend\.fontsize'\s*:\s*(\d+\.?\d*)",
        lambda m: f"'legend.fontsize': {max(10, int(round(float(m.group(1)) * SCALE)))}",
        content
    )
    # figure.dpi
    content = re.sub(
        r"'figure\.dpi'\s*:\s*\d+",
        f"'figure.dpi': {DPI_NEW}",
        content
    )
    return content


def update_savefig_dpi(content):
    """Update dpi= in savefig calls, and add dpi where missing."""
    # Update existing dpi=N
    content = re.sub(
        r'(savefig\([^)]*?)dpi=\d+',
        f'\\1dpi={DPI_NEW}',
        content
    )
    # Add dpi to savefig calls that don't have it
    content = re.sub(
        r"(plt\.savefig\([^)]*)(bbox_inches='tight'\s*\))",
        f"\\1\\2".replace("bbox_inches='tight')", f"dpi={DPI_NEW}, bbox_inches='tight')"),
        content
    )
    return content


def process_script(filepath):
    """Process a single script file."""
    with open(filepath, 'r') as f:
        content = f.read()

    original = content

    # 1. Scale fontsize=N and font_size=N (but not inside rcParams which we handle separately)
    content = re.sub(r'(fontsize=)(\d+\.?\d*)', scale_fontsize, content)
    content = re.sub(r'(font_size=)(\d+\.?\d*)', scale_fontsize, content)

    # 2. Update rcParams
    content = update_rcparams(content)

    # 3. Update savefig DPI
    # First handle existing dpi= params
    content = re.sub(r'(savefig\([^)]*?)dpi=\d+', f'\\1dpi={DPI_NEW}', content)
    # Add dpi to savefig calls without it
    content = re.sub(
        r"(plt\.savefig\([^)]*?)(bbox_inches='tight'\))",
        lambda m: m.group(1) + f"dpi={DPI_NEW}, bbox_inches='tight')" if 'dpi=' not in m.group(1) else m.group(0),
        content
    )

    if content != original:
        with open(filepath, 'w') as f:
            f.write(content)
        return True
    return False


# Process all scripts
print("Scaling fonts across all analysis scripts...")
print(f"  Scale factor: {SCALE}x")
print(f"  Target DPI: {DPI_NEW}")
print()

for script in SCRIPTS:
    path = os.path.join(WORK_DIR, script)
    if os.path.exists(path):
        changed = process_script(path)
        status = "UPDATED" if changed else "no changes"
        print(f"  {script}: {status}")
    else:
        print(f"  {script}: NOT FOUND")

print("\nDone! All scripts updated. Run each script to regenerate figures.")
