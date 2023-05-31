__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2023, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule export_qc_bedtools_intersect:
    input:
        a="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
        coverage_csi="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz.csi",
        b=config["reference"]["design_bed"],
    output:
        results=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.design_bed.txt"),
    params:
        extra=config.get("export_qc_bedtools_intersect", {}).get("extra", ""),
    log:
        "qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.design_bed.log",
    benchmark:
        repeat(
            "qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.design_bed.benchmark.tsv",
            config.get("export_qc_bedtools_intersect", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("export_qc_bedtools_intersect", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("export_qc_bedtools_intersect", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("export_qc_bedtools_intersect", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("export_qc_bedtools_intersect", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("export_qc_bedtools_intersect", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("export_qc_bedtools_intersect", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("export_qc_bedtools_intersect", {}).get("container", config["default_container"])
    message:
        "{rule}: export low cov regions from {input.a} based on {input.b}"
    wrapper:
        "v1.28.0/bio/bedtools/intersect"


rule export_qc_xlsx_report:
    input:
        mosdepth_summary="qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt",
        mosdepth_thresholds="qc/mosdepth_bed/{sample}_{type}.thresholds.bed.gz",
        mosdepth_regions="qc/mosdepth_bed/{sample}_{type}.regions.bed.gz",
        mosdepth_perbase="qc/mosdepth_bed/{sample}_{type}.mosdepth.per-base.design_bed.txt",
        picard_dup="qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt",
        design_bed=config["reference"]["design_bed"],
        pgrs_bed=config["reference"]["design_bed"],
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
    threads: config.get("export_qc_xlsx_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("export_qc_xlsx_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("export_qc_xlsx_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("export_qc_xlsx_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("export_qc_xlsx_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("export_qc_xlsx_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("export_qc_xlsx_report", {}).get("container", config["default_container"])
    message:
        "{rule}: collecting qc values into {output}"
    script:
        "../scripts/export_qc_xlsx_report.py"
