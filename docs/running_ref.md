# Running the references pipeline
To run the reference pipeline you first have to [setup and run Marple](/running) until you have produced `.bam` files for all samples you want to use in your normal pool. 

### :rocket: Marple run command 
To generate `.bam` **and** `.bai`-files for all samples you need to run Marple until a rule that uses the `alignment/samtools_merge_bam/{sample}_{type}.bam` **and** `alignment/samtools_merge_bam/{sample}_{type}.bam.bai` as an input files, e.g. `qc_mosdepth_bed`. Don't forget the `--no-temp` parameter!

```bash
# Run snakemake command with the extra config parameter called sequenceid
snakemake --profile snakemakeprofile --configfile config.yaml --config sequenceid="230202-test" -s /path/to/marple/workflow/Snakefile --no-temp --until qc_mosdepth_bed

```
### :books: Input files 
Four different files need to be available in your runfolder and to be adapted to your compute-environment and sequence run; `samples.tsv`, `units_references.tsv`, `config_references.yaml` and `resources.yaml`.
#### Samples and Units
To run the references pipeline you can use the same `samples.tsv` format as the standard pipeline, but of course only with samples you want in your normalpool. The `units_references.tsv` looks a bit different than the standard format.

| Column id | Description | 
| --- | --- |
|sample| sample name that matches the `samples.tsv`|
|type| data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA  |
|bam| Path to bam file produced by Marple `alignment/samtools_merge_bam/{sample}_{type}.bam`|

#### Config
The bare-bone version of the config file can be found in the `config/config_references.yaml`. This need to be adapted to match the local paths to referencefiles, bedfiles, caches etc on your system. Remember that the `config['reference']['design_bed']` need to be the same bedfile used later in standard Marple (`config['exomedepth_call']['bedfile']`)!


/// details | Expand to view current reference config.yaml
```yaml
{% include "includes/config_references.yaml" %}
```
///

#### Resources
An `resources.yaml` file can also be found in the `config/`-folder. This is adapted to the Uppsala Clinical Genomics' compute cluster but can be used as an indication of resources needed for the different programs. 

### :rocket: Run command 
```bash
#Activate the virtual enviorment
source virtual/environment/bin/activate

# Run snakemake command
snakemake --profile snakemakeprofile -s /path/to/marple/workflow/Snakefile_references --configfile config_references.yaml

```

## Result files
The Marple reference workflow only produces a single result file, a RData file containing the normalpool for exomedepth.
This file is located at `references/exomedepth_reference/RefCount.Rdata` and need to be moved and renamed before added to the standard `config.yaml`.

```yaml
exomedepth_call:
  container: "docker://hydragenetics/exomedepth:1.1.15"
  ref_count: "exomedepth_reference/RefCount_male.Rdata"
...
```
