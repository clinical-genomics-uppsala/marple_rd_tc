__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule merge_fourgene:
    input:
        vcf=expand("Results/{sample}/{sample}.deepsomatic.4genes.vcf.gz", sample=samples),
        tbi=expand("Results/{sample}/{sample}.deepsomatic.4genes.vcf.gz.tbi", sample=samples),
    output:
        vcf="snv_indels/deepsomatic/deepsomatic.4genes.vcf.gz),
    log:
        "snv_indels/deepsomatic/4genes.vcf.log",
    benchmark:
        repeat(
            "snv_indels/deepsomatic/4genes.vcf.benchmark.tsv",
            config.get("merge", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merge all vcfs into one file"
    shell:
        "bcftools merge --merge both {input.vcf} -O z -o {output.vcf}"


rule merge_fourgene_filter:
    input:
        vcf=expand("Results/{sample}/{sample}.deepsomatic.4genes_filter.vcf.gz", sample=samples),
        tbi=expand("Results/{sample}/{sample}.deepsomatic.4genes_filter.vcf.gz.tbi", sample=samples),
    output:
        vcf="snv_indels/deepsomatic/deepsomatic.4genes_filter.vcf.gz),
    log:
        "snv_indels/deepsomatic/4genes.vcf.log",
    benchmark:
        repeat(
            "snv_indels/deepsomatic/4genes.vcf.benchmark.tsv",
            config.get("merge", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mergee", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merge all filtered vcfs into one file"
    shell:
        "bcftools merge --merge both {input.vcf} -O z -o {output.vcf}"

""" 
rule merge_mosaicforecast_phasing:
    input:
        txt=expand("Results/{sample}/{sample}.mosaicforecast.phasing", sample=samples),
    output:
        txt="Results/mosaicforecast.phasing),
    log:
        "snv_indels/mosaicforecast.phasing.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast.phasing.benchmark.tsv",
            config.get("merge", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merge all mosaicforecast phasing into one file"
    shell:
        "bcftools merge --merge both {input.txt} -O z -o {output.txt}"


rule merge_mosaicforecast_SNP:
    input:
        txt=expand("Results/{sample}/{sample}.mosaicforecast.SNP.predictions", sample=samples),
    output:
        txt="Results/mosaicforecast.SNP.predictions),
    log:
        "snv_indels/mosaicforecast.SNP.predictions.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast.SNP.predictions.benchmark.tsv",
            config.get("merge", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merge all mosaicforecast SNP.predictions into one file"
    shell:
        "bcftools merge --merge both {input.txt} -O z -o {output.txt}"


rule merge_mosaicforecast_DEL:
    input:
        txt=expand("Results/{sample}/{sample}.mosaicforecast.DEL.predictions", sample=samples),
    output:
        txt="Results/mosaicforecast.DEL.predictions),
    log:
        "snv_indels/mosaicforecast.DEL.predictions.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast.DEL.predictions.benchmark.tsv",
            config.get("merge", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merge all mosaicforecast DEL.predictions into one file"
    shell:
        "bcftools merge --merge both {input.txt} -O z -o {output.txt}"


rule merge_mosaicforecast_INS:
    input:
        txt=expand("Results/{sample}/{sample}.mosaicforecast.INS.predictions", sample=samples),
    output:
        txt="Results/mosaicforecast.INS.predictions),
    log:
        "snv_indels/mosaicforecast.INS.predictions.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast.INS.predictions.benchmark.tsv",
            config.get("merge", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merge all mosaicforecast INS.predictions into one file"
    shell:
        "bcftools merge --merge both {input.txt} -O z -o {output.txt}"


rule merge_deepmosaic:
    input:
        txt=expand("Results/{sample}/{sample}.deepmosaic.txt", sample=samples),
    output:
        txt="Results/deepmosaic.txt),
    log:
        "snv_indels/deepmosaic.log",
    benchmark:
        repeat(
            "snv_indels/deepmosaic.benchmark.tsv",
            config.get("merge", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("merge", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge", {}).get("container", config["default_container"])
    message:
        "{rule}: merge all mosaicforecast deepmosaic into one file"
    shell:
        "bcftools merge --merge both {input.txt} -O z -o {output.txt}"
 """
