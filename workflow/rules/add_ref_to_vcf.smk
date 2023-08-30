__author__ = "Jessika Nordin, Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule add_ref_to_vcf:
    input:
        vcf="snv_indels/haplotypecaller/{sample}_N.normalized.sorted.vep_annotated.vcf.gz",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/haplotypecaller/{sample}_N.normalized.sorted.vep_annotated.ref.vcf"),
    log:
        "snv_indels/haplotypecaller/{sample}_N.normalized.sorted.vep_annotated.ref.vcf.log",
    benchmark:
        repeat(
            "snv_indels/haplotypecaller/{sample}_N.normalized.sorted.vep_annotated.ref.vcf.benchmark.tsv",
            config.get("add_ref_to_vcf", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("add_ref_to_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("add_ref_to_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("add_ref_to_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("add_ref_to_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("add_ref_to_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("add_ref_to_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: Add reference to the header of the haplotypecaller vcf: {input.vcf}"
    script:
        "../scripts/add_ref_to_vcf.py"
