__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule sample_prep:
    input:
        sample="config/samples.tsv",
    output:
        sex=temp("config/sample_sex.tsv"),
    log:
        temp("config/sample_prep.tsv.log"),
    threads: config.get("sample_prep", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sample_prep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sample_prep", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sample_prep", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sample_prep", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sample_prep", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sample_prep", {}).get("container", config["default_container"])
    message:
        "{rule}: Create tsv with sample and column with sex (all female) to run deepmosaic"
    shell:
        """(awk -F="\t" '{print $1, $2, "female"}' {input.samples} > {output.sex}) &> {log}"""
