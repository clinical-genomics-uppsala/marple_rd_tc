# :snake: Welcome to Marple - Inherited Cancer pipeline :woman_detective:

This pipeline is created to run on short-read Illumina data from a custom Twist Inherited Cancer panel, designed at [Clinical Genomics Uppsala](https://www.scilifelab.se/units/clinical-genomics-uppsala/#https://www.cgu.igp.uu.se).

This snakemake pipeline uses the module system from [Hydra Genetics](https://github.com/hydra-genetics/) to process `.fastq.gz` files. The pipeline produces a MultiQC `.html` report with QC-data, `.bam` alignment files, annotated `.vcf.gz` for SNVs and smaller indels, as well as `.txt` and `.aed` files for structural variants from Exomedepth. 

<br />
Marple :woman_detective: uses the following hydra genetics modules:

- [Alignment](https://github.com/hydra-genetics/alignment/tree/v0.4.0)
- [Annotation](https://github.com/hydra-genetics/annotation/tree/v0.3.0)
- [CNV](https://github.com/hydra-genetics/cnv_sv/tree/78f270c)
- [Prealignment](https://github.com/hydra-genetics/prealignment/tree/v1.0.0)
- [SNV indels](https://github.com/hydra-genetics/snv_indels/tree/v0.3.0)
- [QC](https://github.com/hydra-genetics/qc/tree/ca947b1)


## :judge: Rulegraph 
![dag plot](includes/rulegraph.svg){: style="height:100%;width:100%"}

---
# Hydra-genetics

We are an organization/community with the goal of making [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline development easier, faster, a bit more structured and of higher quality.

We do this by providing [snakemake modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules) that can be combined to create a complete analysis or included in already existing pipelines. All modules are subjected to extensive testing to make sure that new releases doesn't unexpectedly break existing pipeline or deviate from guidelines and best practices on how to write code.

There is also a small [tutorial](https://hydra-genetics.readthedocs.io/en/latest/simple_pipeline/) available to help you get started with Hydra-genetics.

# Snakemake
Marple and Hydra-genetics are snakemake bases pipeline/tools. The [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow management system is a tool to create reproducible and scalable data analyses. Workflows are described via a human readable, Python based language. They can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment. 

If Snakemake is new to you a good place to start is doing the [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) since this will help you setting Marple up.

