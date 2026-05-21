import os
import toml
from pathlib import Path

cwd = os.getcwd()

toml_config = toml.load(f"{cwd}/scripts/config_final.toml")

samples = toml_config["general"]["samples"]
analyses = toml_config["general"]["analysis"]
flowcells = toml_config["general"]["fc_dir_names"]

output = toml_config["general"]["project_path"]
name = output.rstrip("/").split("/")[-2].split("_", 1)[1]

# Check alignments
## In flowcells. Folder "alignments" should not exist.
for f in flowcells:
    path = f"{cwd}/{f}/alignments"
    if not os.path.isdir(path):
        print(f"Directory '{f}/alignments' was previously removed to free space.")
    elif not os.listdir(path):
        print(f"Directory '{f}/alignments' exists but is empty.")
    elif path.is_dir() and any(path.iterdir()):
        print(
            f"WARNING! Directory '{f}/alignments' exists and is NOT empty! Please remove directory."
        )

## In alignments. Check bam for all sample
for s in samples:
    align_dir = Path(f"{cwd}/alignments")
    bam_file = align_dir / f"{s}_sorted.bam"
    bai_file = align_dir / f"{s}_sorted.bam.bai"

    if bam_file.is_file() and bai_file.is_file():
        print(f"{s} alignments files found!")
    else:
        print(f"WARNING! No alignemnt for {s}!")

# Check QC
mosdepth = Path(f"{cwd}/qc/{name}.html")
if mosdepth.is_file():
    print(f"Mosdepth summary found : 'qc/{name}.html'")

# Check results
for s in samples:
    print(f"Checking results for {s}")
    all_exist = True

    base_dir = Path(f"{cwd}/results")

    ## General EPI2ME
    if (
        "SNP" in toml_config["general"]["analysis"]
        or "SV" in toml_config["general"]["analysis"]
        or "CNV" in toml_config["general"]["analysis"]
        or "phasing" in toml_config["general"]["analysis"]
        or "repeats" in toml_config["general"]["analysis"]
        or "methylation" in toml_config["general"]["analysis"]
    ):
        epi2me_files = [
            f"{s}.flagstat.tsv",
            f"{s}.mosdepth.global.dist.txt",
            f"{s}.mosdepth.summary.txt",
            f"{s}.readstats.tsv.gz",
            f"{s}.regions.bed.gz",
            f"{s}.stats.json",
            f"{s}.thresholds.bed.gz",
            f"{s}.wf-human-alignment-report.html",
        ]

        for file_name in epi2me_files:
            matching_files = list(base_dir.glob(f"*/{s}/{file_name}"))
            if matching_files:
                continue
            else:
                print(f"\tWARNING! Missing: {matching_files}")
                all_exist = False

    ## SNP
    if "SNP" in toml_config["general"]["analysis"]:
        snp_files = [
            f"{s}.wf-human-snp-report.html",
            f"{s}.wf_snp_clinvar.vcf.gz",
            f"{s}.wf_snp_clinvar.vcf.gz.tbi",
            f"{s}.wf_snp.vcf.gz",
            f"{s}.wf_snp.vcf.gz.tbi",
            f"{s}_deepvariant.vcf",
            f"{s}_deepvariant.visual_report.html",
        ]

        for file_name in snp_files:
            matching_files = list(base_dir.glob(f"*/{s}/{file_name}"))
            if matching_files:
                continue
            else:
                print(f"\tWARNING! Missing: {matching_files}")
                all_exist = False

    ## SV
    if "SV" in toml_config["general"]["analysis"]:
        sv_files = [
            f"{s}.wf-human-sv-report.html",
            f"{s}.wf_sv.snf",
            f"{s}.wf_sv.vcf.gz",
            f"{s}.wf_sv.vcf.gz.tbi",
            f"{s}_cutesv.vcf",
        ]

        for file_name in sv_files:
            matching_files = list(base_dir.glob(f"*/{s}/{file_name}"))
            if matching_files:
                continue
            else:
                print(f"\tWARNING! Missing: {matching_files}")
                all_exist = False

    ## CNV
    if "CNV" in toml_config["general"]["analysis"]:
        cnv_files = [
            f"{s}.wf_cnv.vcf.gz",
            f"{s}.wf_cnv.vcf.gz.tbi",
            f"{s}.wf-human-cnv-report.html",
        ]

        for file_name in cnv_files:
            matching_files = list(base_dir.glob(f"*/{s}/{file_name}"))
            if matching_files:
                continue
            else:
                print(f"\tWARNING! Missing: {matching_files}")
                all_exist = False

    ## repeats
    if "repeats" in toml_config["general"]["analysis"]:
        repeats_files = [
            f"{s}.wf-human-str-report.html",
            f"{s}.wf_str.straglr.tsv",
            f"{s}.wf_str.vcf.gz",
            f"{s}.wf_str.vcf.gz.tbi",
            f"{s}.sorted.vcf.gz",
            f"{s}.sorted.vcf.gz.csi",
            f"{s}sample.spanning.sorted.bam",
            f"{s}.spanning.sorted.bam.bai",
            f"{s}_last.json",
            f"{s}_last.vcf",
        ]

        for file_name in repeats_files:
            matching_files = list(base_dir.glob(f"*/{s}/{file_name}"))
            if matching_files:
                continue
            else:
                print(f"\tWARNING! Missing: {matching_files}")
                all_exist = False

    ## methylation
    if "methylation" in toml_config["general"]["analysis"]:
        methylation_files = [
            f"{s}.wf_mods.5mC.bw",
            f"{s}.wf_mods.bedmethyl.gz",
        ]

        for file_name in methylation_files:
            matching_files = list(base_dir.glob(f"*/{s}/{file_name}"))
            if matching_files:
                continue
            else:
                print(f"\tWARNING! Missing: {matching_files}")
                all_exist = False

    ## phasing
    if "phasing" in toml_config["general"]["analysis"]:
        phasing_files = [
            f"{s}.haplotagged.cram",
            f"{s}.haplotagged.cram.crai",
            f"{s}_haplotype.txt",
            f"{s}_haplotype.phased.vcf.gz",
            f"{s}_haplotype.phased.vcf.gz.tbi",
        ]
        if "SNP" in toml_config["general"]["analysis"]:
            phasing_files.append(
                [
                    f"{s}.wf_snp.haploblocks.gtf",
                ]
            )
        if "methylation" in toml_config["general"]["analysis"]:
            phasing_files.append(
                [
                    f"{s}.wf_mods.1-5mC.bw",
                    f"{s}.wf_mods.1.bedmethyl.gz",
                    f"{s}.wf_mods.2-5mC.bw",
                    f"{s}.wf_mods.2.bedmethyl.gz",
                    "dmrs_table.bed",
                    "dmrs_table_annotated.bed",
                    "annotation_summary.log",
                    "dmr_status.log",
                    "dmr_summary_report.html",
                    "dmr_summary_stats.tsv",
                ]
            )

        for file_name in phasing_files:
            matching_files = list(base_dir.glob(f"*/{s}/{file_name}"))
            if matching_files:
                continue
            else:
                print(f"\tWARNING! Missing: {matching_files}")
                all_exist = False

    ## splicing
    if "splicing" in toml_config["general"]["analysis"]:
        splicing_files = [
            f"{s}_junctions.bed",
            f"{s}_intron_calls.tsv",
            f"{s}.isoform.read.map.txt",
            f"{s}.isoforms.bed",
            f"{s}.isoforms.fa",
            f"{s}.isoforms.gtf",
            f"{s}.counts.tsv",
            f"{s}.isoform.counts.txt",
        ]

        for file_name in splicing_files:
            matching_files = list(base_dir.glob(f"*/{s}/{file_name}"))
            if matching_files:
                continue
            else:
                print(f"\tWARNING! Missing: {matching_files}")
                all_exist = False

    if all_exist:
        print("\tAll results files are present!")
