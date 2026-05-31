__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2023, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule sample_order_multiqc:
    output:
        replacement=temp("qc/multiqc/sample_replacement.tsv"),
        order=temp("qc/multiqc/sample_order.tsv"),
    params:
        filelist=[(u.sample, u.fastq1) for u in units[units.type == "N"].itertuples()],
    log:
        "qc/multiqc/sample_order.tsv.log",
    benchmark:
        repeat("qc/multiqc/sample_order.tsv.benchmark.tsv", config.get("sample_order_multiqc", {}).get("benchmark_repeats", 1))
    threads: rule_resource("sample_order_multiqc", "threads")
    resources:
        mem_mb=rule_resource("sample_order_multiqc", "mem_mb"),
        mem_per_cpu=rule_resource("sample_order_multiqc", "mem_per_cpu"),
        partition=rule_resource("sample_order_multiqc", "partition"),
        threads=rule_resource("sample_order_multiqc", "threads"),
        time=rule_resource("sample_order_multiqc", "time"),
    container:
        rule_container("sample_order_multiqc")
    message:
        "{rule}: Create a sample order tsv based on S_index in {params.filelist} for multiqc"
    script:
        "../scripts/sample_order_multiqc.py"
