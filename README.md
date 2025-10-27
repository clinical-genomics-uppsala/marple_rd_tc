# :snake: clinical-genomics-uppsala/marple :female_detective:

#### Twist Cancer inherited hg38 hydra pipelines
---

<p align="center">
<a href="https://marple-rd-tc.readthedocs.io/en/latest/">https://marple-rd-tc.readthedocs.io/en/latest/</a>
</p>

This ReadMe is only a brief introduction, please refer to ReadTheDocs for the latest documentation. 

---
![Lint](https://github.com/clinical-genomics-uppsala/marple_rd_tc/actions/workflows/lint.yaml/badge.svg?branch=main)
![Snakefmt](https://github.com/clinical-genomics-uppsala/marple_rd_tc/actions/workflows/snakefmt.yaml/badge.svg?branch=main)
![snakemake dry run](https://github.com/clinical-genomics-uppsala/marple_rd_tc/actions/workflows/snakemake-dry-run.yaml/badge.svg?branch=main)
![integration test](https://github.com/clinical-genomics-uppsala/marple_rd_tc/actions/workflows/integration1.yaml/badge.svg?branch=main)

![pycodestyle](https://github.com/clinical-genomics-uppsala/marple_rd_tc/actions/workflows/pycodestyle.yaml/badge.svg?branch=main)
![pytest](https://github.com/clinical-genomics-uppsala/marple_rd_tc/actions/workflows/pytest.yaml/badge.svg?branch=main)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction
This pipeline is created to run on Illumina data from a custom Twist Inherited Cancer panel, designed at Clinical Genomics Uppsala.

## :judge: Rule Graph
![rule_graph](images/rulegraph.svg)

## :white_check_mark: Testing

The workflow repository contains a small test dataset (:exclamation: Todo: as of now dry-run only) `.tests/integration` which can be run like so:

```bash
$ cd .tests/integration
$ snakemake -n -s ../../workflow/Snakefile --configfiles ../../config/config.yaml config.yaml config_exomedepth.yaml --config sequenceid="990909_test" PATH_TO_REPO=/folder/containing/marple_rd_tc/
```
> **_NOTE:_**   If using the variable `PATH_TO_REPO` in the config-file this need to be defined in the commandline


## :rocket: [Usage](https://marple-rd-tc.readthedocs.io/en/latest/running/)

To run this pipeline `sample.tsv`, `units.tsv`, `resources.yaml`, and `config.yaml` files need to be available in the current directory (or otherwise specified in `config.yaml`). A config file specifying a path to an exomedepth reference file must also be given (`config_exomedepth_miseq.yaml` in the example below). You always need to specify the `config`-file and `sequenceid` variable in the command. To run the pipeline:

```bash
$ snakemake --profile snakemakeprofile --configfiles config.yaml config_exomedepth_miseq.yaml --config sequenceid="990909_test" -s /path/to/marple_rd_tc/workflow/Snakefile --config PATH_TO_REPO=/folder/containing/marple_rd_tc/
```

> **_NOTE:_**   If using the variable `PATH_TO_REPO` in the config this need to be defined in the commandline


## :books: [Output files](https://marple-rd-tc.readthedocs.io/en/latest/result_files/)

The following output files are located in `Results/`-folder:

| File | Format |Description |
|---|---|---|
|`multiqc_DNA.html` | html | Aggregated QC values for entire sequence run, open in browser |
|`{sample}/{sample}.xlsx`| xlsx | Excel file with QC stats (primarily coverage) for each sample|
|`{sample}/{sample}_N.cram"`| cram | Deduplicated alignment file |
|`{sample}/{sample}_N.cram.crai`| crai | Index for alignment file|
|`{sample}/{sample}.hard-filtered.vcf.gz`| vcf.gz | Compressed VCF-file decomposed, normalized and annotated with vep |
|`{sample}/{sample}.hard-filtered.vcf.gz.tbi`| tbi | Index for variant file |
|`{sample}/{sample}.genome.vcf.gz`| genome.vcf.gz | Compressed VCF-file for all positions in the design, not decomposed nor normalized |
|`{sample}/{sample}.genome.vcf.gz.tbi`| tbi | Index for genome VCF-file |
|`{sample}/{sample}_exomedepth_SV.txt`| txt | Nexus SV text file with structural variants from ExomeDepth |
|`{sample}/{sample}_exomedepth.aed`| aed | aed text file with structural variants from ExomeDepth |
|`{sample}/{sample}.cnv.vcf.gz`| vcf.gz | Compressed VCF-file with structural variants from ExomeDepth |
|`{sample}/{sample}.cnv.vcf.gz.tbi`| tbi | Index for variant file from ExomeDepth |
|`{sample}/mobile_elements/{sample}.ALU.vcf.gz`| vcf.gz | Compressed VCF-file with predicted ALU elements |
|`{sample}/mobile_elements/{sample}.LINE1.vcf.gz`| vcf.gz | Compressed VCF-file with predicted LINE1 elements |
|`{sample}/mobile_elements/{sample}.HERVK.vcf.gz`| vcf.gz | Compressed VCF-file with predicted HERVK elements |
|`{sample}/mobile_elements/{sample}.SVA.vcf.gz`| vcf.gz | Compressed VCF-file with predicted SVA elements |
|`{sample}/mosaic/{sample}.deepmosaic.txt`| tsv | Candidate variants and their predictions from DeepMosaic |
|`{sample}/mosaic/{sample}.deepsomatic.vcf.gz`| vcf.gz | Compressed VCF-file from DeepSomatic where PASS are possible mosaic variants|
|`{sample}/mosaic/{sample}.deepsomatic.vcf.gz.tbi`| vcf.gz | Index for genome VCF-file |
|`{sample}/mosaic/{sample}.mosaicforecast.phasing`| tsv | Candidate mosaic variants based on phasing from MosaicForecast |
|`{sample}/mosaic/{sample}.mosaicforecast.DEL.predictions`| tsv | Candidate deletion variants and their predictions from MosaicForecast |
|`{sample}/mosaic/{sample}.mosaicforecast.INS.predictions`| tsv | Candidate insertion variants and their predictions from MosaicForecast |
|`{sample}/mosaic/{sample}.mosaicforecast.SNP.predictions`| tsv | Candidate SNP variants and their predictions from MosaicForecast |
|`{sequenceid}_config.yaml`| yaml | yaml config-file with programversion and extra settings used |
|`{sequenceid}_config_exomedepth.yaml`| yaml | yaml config-file with which reference was used for ExomeDepth |
