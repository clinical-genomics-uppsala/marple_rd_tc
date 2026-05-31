__author__ = "Jessika Nordin, Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule add_ref_to_vcf:
    input:
        vcf="snv_indels/deepvariant/{sample}_N.normalized.sorted.vep_annotated.vcf.gz",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/deepvariant/{sample}_N.normalized.sorted.vep_annotated.ref.vcf"),
    log:
        "snv_indels/deepvariant/{sample}_N.normalized.sorted.vep_annotated.ref.vcf.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_N.normalized.sorted.vep_annotated.ref.vcf.benchmark.tsv",
            config.get("add_ref_to_vcf", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=rule_resource("add_ref_to_vcf", "mem_mb"),
        mem_per_cpu=rule_resource("add_ref_to_vcf", "mem_per_cpu"),
        partition=rule_resource("add_ref_to_vcf", "partition"),
        threads=rule_resource("add_ref_to_vcf", "threads"),
        time=rule_resource("add_ref_to_vcf", "time"),
    container:
        rule_container("add_ref_to_vcf")
    message:
        "{rule}: Add reference to the header of the deepvariant vcf: {input.vcf}"
    script:
        "../scripts/add_ref_to_vcf.py"
