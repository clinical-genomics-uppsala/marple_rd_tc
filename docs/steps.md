# Steps in Marple :woman_detective:
To go into details of the pipeline we dived the pipeline into modules similar to Hydra-Genetics module system. Default hydra-genetics settings/resources are used if no configuration is specified.

---
## Prealignment
See the **Prealignment** hydra-genetics module documentation on [ReadTheDoc](https://hydra-genetics-prealignment.readthedocs.io/en/latest/) or [github](https://github.com/hydra-genetics/prealignment/tree/v1.4.0) documentation for more details on the softwares. 

![dag plot](includes/images/prealignment.png){: style="height:30%;width:30%"}

### Pipeline output files
Only temporary intermediate files are created.

### Trimming
Trimming of fastq files is performed by **[fastp](https://github.com/OpenGene/fastp)** v0.20.1.  

### Merging
Merging of fastq files belonging to the same sample are performed by simply concatenating the files with **cat** if applicable.

---
## Alignment
See **Alignment** hydra-genetics module on [ReadTheDocs](https://hydra-genetics-alignment.readthedocs.io/en/latest/) or [github](https://github.com/hydra-genetics/alignment/tree/v0.7.0) for more information and documentation on the softwares used in the alignment steps. 

![dag plot](includes/images/alignment.png){: style="height:100%;width:100%"}

### Pipeline output files:

* `Results/{sample}/{sample}.cram`
* `Results/{sample}/{sample}.cram.crai`

### Alignment with BWA-mem
Alignment of fastq files into bam files is performed by **[bwa-mem](https://github.com/lh3/bwa)** v0.7.17 using the trimmed fastq files. This make it possible to speed up alignment by utilizing parallelization and also make it possible to analyze qc for lanes separately. Bamfiles are then directly sorted by **[samtools sort](http://www.htslib.org/doc/samtools-sort.html)** v1.15.

#### Read groups
Bam file read groups are set according to sequencing information in the `units.tsv` file.
The @RG read tag is set using the following options defined in the hydra-genetics bwa rule:
```
-R '@RG\tID:{ID}\tSM:{SM}\tPL:{PL}\tPU:{PU}\tLB:{LB}' -v 1
```
where the individual read groups are defined below:

| **RG tag** | **Value** |
|-------------|-|
| ID | sample_type.lane.barcode |
| SM | sample_type |
| PL | platform |
| PU | flowcell.lane.barcode |
| LB | sample_type |

### Bam splitting
The bam files are split into chromosome files for faster performance in downstream analysis. Split files are used by markduplicates. Splitting is performed by **[samtools view](http://www.htslib.org/doc/samtools-view.html)** v1.15.

### Mark duplicates
Flagging duplicated reads are performed on individual chromosome bam files by **[picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)** v2.25.4.

### Merging
Merging of deduplicated bam files belonging to the same sample are performed by **[samtools merge](http://www.htslib.org/doc/samtools-merge.html)** v1.15.

### Sorting
Merged bamfile are sorted by **[samtools sort](http://www.htslib.org/doc/samtools-sort.html)** v1.15.

### Bam indexing
Bamfile indexing is performed by **[samtools index](http://www.htslib.org/doc/samtools-index.html)** v1.15.

---
## SNV indels
SNV and indels are called using the **SNV_indels** ([ReadTheDoc](https://hydra-genetics-snv-indels.readthedocs.io/en/latest/) or [github](https://github.com/hydra-genetics/snv_indels/tree/v1.3.0)) module. Annotation is then done with **Annotation** module ([ReadTheDocs](https://hydra-genetics-annotation.readthedocs.io/en/latest/) or [github](https://github.com/hydra-genetics/annotation/tree/v1.2.0)).

![dag plot](includes/images/snv_indels.png){: style="height:100%;width:100%"}


### Pipeline output files

* `Results/{sample}/{sample}.hard-filtered.vcf.gz`
* `Results/{sample}/{sample}.genome.vcf.gz`
* `Results/{sample}/mosaic/{sample}.deepmosaic.txt`
* `Results/{sample}/mosaic/{sample}.deepsomatic.vcf.gz`
* `Results/{sample}/mosaic/{sample}.mosaicforecast.phasing`
* `Results/{sample}/mosaic/{sample}.mosaicforecast.DEL.predictions`
* `Results/{sample}/mosaic/{sample}.mosaicforecast.INS.predictions`
* `Results/{sample}/mosaic/{sample}.mosaicforecast.SNP.predictions`

### SNV calling
Variants are called using [DeepVariant](https://github.com/google/deepvariant). Deepvariant is run per chromosome over the regions defined in `config["references"]["design_bed"]`. The model type is set to "WES" and `--output_gvcf` is used to ensure that both genome vcf as well as a standard vcf is produced. The AF field is added to the `INFO` column in the vcf:s using the `fix_af.py`. The vcf header in the standard vcf is also updated to include a reference line using the `add_ref_to_vcf.py` to ensure that programs such as Alissa acknowledge the use of Hg38.

### Normalizing
The standard vcf files is decomposed with [**vt decompose**](https://genome.sph.umich.edu/wiki/Vt#Decompose) followed by [**vt decompose_blocksub**](https://genome.sph.umich.edu/wiki/Vt#Decompose_biallelic_block_substitutions) v2015.11.10. The decomposed vcf files are then normalized by [**vt normalize**](https://genome.sph.umich.edu/wiki/Vt#Normalization) v2015.11.10.

### Annotation
The normalized standard VCF files are then annotated using **[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)** v109.3. Vep is run with the extra parameters `--assembly GRCh38 --check_existing --pick --variant_class --everything`.

See the [annotation hydra-genetics module](https://hydra-genetics-annotation.readthedocs.io/en/latest/) for additional information.

### Mosaic variant calling
Possible mosaic variants is called by [DeepSomatic](https://github.com/google/deepsomatic) which are then predicted to be real or not with [MosaicForecast](https://github.com/parklab/MosaicForecast) and [DeepMosaic](https://github.com/XiaoxuYangLab/DeepMosaic).


---
## CNV
CNVs are called using the Hydra-Genetics **CNV_SV** module ([ReadTheDocs](https://hydra-genetics-cnv-sv.readthedocs.io/en/latest/) or [github](https://github.com/hydra-genetics/cnv_sv/tree/78f270c)).
<br />

![dag plot](includes/images/cnv_sv.png){: style="height:60%;width:60%"}

### Pipeline output files

* `Results/{sample}_{sequenceid}/{sample}_{sequenceid}_exomedepth_SV.txt`
* `Results/{sample}_{sequenceid}/{sample}_{sequenceid}_exomedepth.aed`
* `Results/{sample}/{sample}.cnv.vcf.gz`
* `Results/{sample}/mobile_elements/{sample}.ALU.vcf.gz`
* `Results/{sample}/mobile_elements/{sample}.LINE1.vcf.gz`
* `Results/{sample}/mobile_elements/{sample}.HERVK.vcf.gz`
* `Results/{sample}/mobile_elements/{sample}.SVA.vcf.gz`

### Exomedepth
To call larger structural variants **[Exomedepth](https://cran.r-project.org/web/packages/ExomeDepth/index.html)** v1.1.15 is used. Exomedepth does **not** use a window approach but evaluates each row in the bedfile as a segment, therefor the bedfile need to be split into appropriate large windows (e.g. using `bedtools makewindows`). Exomedepth also need a `RData` file containing the normal pool, this can be created using the [Marple - references workflow](/running_ref). Lines with no-change calls (`reads.ratio == 1`) are removed from the output for Alissa compatibility. Since no sex-chromosome are included in the design Exomedepth is run with the same normalpool irregardless of the sample's biological sex. Marple is designed on HG38 therefor a genes and exonfile are also needed for annotation. 

### Mobile elements
To find mobile elements, like ALU, HERVK, LINE1 and SVA, [MELT](https://melt.igs.umaryland.edu/index.php) is used. It is free to use for academic purposes and you can buy a licence if used in other cases, like for this pipeline as part of the clinical workflow at the hospital. Since you should ask to use MELT, the singularity we use are local.

---
## QC
For quality control the **QC** module ([ReadTheDocs](https://hydra-genetics-qc.readthedocs.io/en/latest/) or [github](https://github.com/hydra-genetics/qc/tree/ca947b1)) is used and the results are summarized/aggregated into a MultiQC-report.

<br />

![dag plot](includes/images/qc.png){: style="height:70%;width:70%"}

### Pipeline output files
* `Results/{sequenceid}_MultiQC.html`

### MultiQC
A MultiQC html report is generated using **[MultiQC](https://github.com/ewels/MultiQC)** v1.11. The report starts with a general statistics table showing the most important QC-values followed by additional QC data and diagrams. See [result files page](/result_files) for more in-depth info on general stats columns. The qc data is collected and generated using Fastp, FastQC, samtools, picard and Mosdepth.

The samples in the MultiQC report is renamed and ordered based on the order of the `SampleSheet.csv` to enable easier integration with the wet lab follow-up.


The report is configured based on a MultiQC config file. 

/// details | Expand to view current MultiQC config.yaml
```yaml
{% include "includes/multiqc_config.yaml" %}
```
///


### FastQC
**[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** v0.11.9 is run on the raw fastq-files before trimming.

### Mosdepth
**[Mosdepth](https://github.com/brentp/mosdepth)** v0.3.2 is used together with a bedfile covering all coding exons (`config[reference][exon_bed]`) and thresholds (`10,20,50`) to calculate coverage.

### Samtools
**[Samtools stats](http://www.htslib.org/doc/samtools-stats.html)** v1.15 is run on BWA-mem aligned and merged bam files without any designfile.

### Picard
**[Picard](https://broadinstitute.github.io/picard/)** v2.25.4 is run on BWA-mem aligned and merged bam files collecting a number of metrics. The metrics calculated are listed below:

* **picard CollectAlignmentSummaryMetrics** - using fasta reference genome file
* **picard CollectDuplicationMetrics**
* **picard CollectGCBiasMetrics**
* **picard CollectHsMetrics** - using fasta reference genome file, the full capture design bedfile (`config[reference][design_bed]`), and with the option COVERAGE_CAP=5000
* **picard CollectInsertSizeMetrics**
* **picard CollectMultipleMetrics** - using fasta reference genome file, and the full bedfile (`config[reference][design_intervals]`)
---
