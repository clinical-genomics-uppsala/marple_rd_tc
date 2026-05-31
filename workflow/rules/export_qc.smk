__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2023, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule export_qc_bedtools_intersect:
    input:
        left="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
        coverage_csi="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz.csi",
        right=config["reference"]["exon_bed"],
    output:
        results=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.txt"),
    params:
        extra=config.get("export_qc_bedtools_intersect", {}).get("extra", ""),
    log:
        "qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.log",
    benchmark:
        repeat(
            "qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.benchmark.tsv",
            config.get("export_qc_bedtools_intersect", {}).get("benchmark_repeats", 1),
        )
    threads: rule_resource("export_qc_bedtools_intersect", "threads")
    resources:
        mem_mb=rule_resource("export_qc_bedtools_intersect", "mem_mb"),
        mem_per_cpu=rule_resource("export_qc_bedtools_intersect", "mem_per_cpu"),
        partition=rule_resource("export_qc_bedtools_intersect", "partition"),
        threads=rule_resource("export_qc_bedtools_intersect", "threads"),
        time=rule_resource("export_qc_bedtools_intersect", "time"),
    container:
        rule_container("export_qc_bedtools_intersect")
    message:
        "{rule}: export low cov regions from {input.left} based on {input.right}"
    wrapper:
        "v1.32.0/bio/bedtools/intersect"


rule export_qc_bedtools_intersect_pgrs:
    input:
        left="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
        coverage_csi="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz.csi",
        right=config["reference"]["pgrs_bed"],
    output:
        results=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.pgrs_cov.txt"),
    params:
        extra=config.get("export_qc_bedtools_intersect", {}).get("extra", ""),
    log:
        "qc/mosdepth_bed/{sample}_{type}.mosdepth.pgrs_cov.log",
    benchmark:
        repeat(
            "qc/mosdepth_bed/{sample}_{type}.mosdepth.pgrs_cov.benchmark.tsv",
            config.get("export_qc_bedtools_intersect", {}).get("benchmark_repeats", 1),
        )
    threads: rule_resource("export_qc_bedtools_intersect", "threads")
    resources:
        mem_mb=rule_resource("export_qc_bedtools_intersect", "mem_mb"),
        mem_per_cpu=rule_resource("export_qc_bedtools_intersect", "mem_per_cpu"),
        partition=rule_resource("export_qc_bedtools_intersect", "partition"),
        threads=rule_resource("export_qc_bedtools_intersect", "threads"),
        time=rule_resource("export_qc_bedtools_intersect", "time"),
    container:
        rule_container("export_qc_bedtools_intersect")
    message:
        "{rule}: export low cov regions from {input.left} based on {input.right}"
    wrapper:
        "v1.32.0/bio/bedtools/intersect"


rule export_qc_xlsx_report:
    input:
        mosdepth_summary="qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt",
        mosdepth_thresholds="qc/mosdepth_bed/{sample}_{type}.thresholds.bed.gz",
        mosdepth_regions="qc/mosdepth_bed/{sample}_{type}.regions.bed.gz",
        mosdepth_perbase="qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.exon_bed.txt",
        picard_dup="qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt",
        pgrs_coverage="qc/mosdepth_bed/{sample}_{type}.mosdepth.pgrs_cov.txt",
        design_bed=config["reference"]["design_bed"],
        pgrs_bed=config["reference"]["pgrs_bed"],
        wanted_transcripts=config["export_qc_xlsx_report"]["wanted_transcripts"],
    output:
        results=temp("qc/xlsx_report/{sample}_{type}.xlsx"),
    params:
        coverage_thresholds=config["mosdepth_bed"]["thresholds"],
        sequenceid=config["sequenceid"],
    log:
        "qc/xlsx_report/{sample}_{type}.xlsx.log",
    benchmark:
        repeat(
            "qc/xlsx_report/{sample}_{type}.xlsx.benchmark.tsv",
            config.get("export_qc_xlsx_report", {}).get("benchmark_repeats", 1),
        )
    threads: rule_resource("export_qc_xlsx_report", "threads")
    resources:
        mem_mb=rule_resource("export_qc_xlsx_report", "mem_mb"),
        mem_per_cpu=rule_resource("export_qc_xlsx_report", "mem_per_cpu"),
        partition=rule_resource("export_qc_xlsx_report", "partition"),
        threads=rule_resource("export_qc_xlsx_report", "threads"),
        time=rule_resource("export_qc_xlsx_report", "time"),
    container:
        rule_container("export_qc_xlsx_report")
    message:
        "{rule}: collecting qc values into {output}"
    # localrule: True
    script:
        "../scripts/export_qc_xlsx_report.py"
