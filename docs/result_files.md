# Result files
Marple produces a lot of intermediate and result files but only files defined in [`output_files.yaml`](https://github.com/clinical-genomics-uppsala/marple_rd_tc/blob/main/config/output_files.yaml) are kept, the rest are temporary and will be deleted when not needed in the any consecutive rules. If other files than the predefined are wanted you need to edit `output_files.yaml` or add `--no-temp` to the running command.

## Files
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


## MultiQC report
Marple produces a **[MultiQC](https://github.com/ewels/MultiQC)**-report for the entire sequencing run to enable easier QC tracking. The report starts with a general statistics table showing the most important QC-values followed by additional QC data and diagrams. The entire MultiQC html-file is interactive and you can filter, highlight, hide or export data using the ToolBox at the right edge of the report.

<br />

The report is configured based on a MultiQC config file. 

/// details | Expand to view current MultiQC config.yaml
```yaml
{% include "includes/multiqc_config.yaml" %}
```
///

### General Statistics
The general statistics table are ordered based on the fastq-file  "S"-index, e.g. `sampleT_S1_R1_001.fastq.gz` will be before `sampleA_S2_R1_001.fastq.gz`. This is done by renaming the samples in two steps using the script `sample_order_multiqc.py`. To toggle between "Sample Order" and "Sample Name" use the buttons just above General Stats header.

<br />

| Column Name | Origin | Comment | 
| --- | --- | --- |
| K Reads | [Samtools stats](http://www.htslib.org/doc/samtools-stats.html)  | Total number of reads in inputfile (`alignment/samtools_merge_bam/{sample}_{type}.bam`) |
| % Mapped| [Samtools stats](http://www.htslib.org/doc/samtools-stats.html) | Percent reads mapped, anywhere in the reference (no design file used) |
| % Proper pairs| [Samtools stats](http://www.htslib.org/doc/samtools-stats.html) | Only reads on target (`config[reference][design_bed]`) |
| Average Quality | [Samtools stats](http://www.htslib.org/doc/samtools-stats.html) | Ratio between sum of base quality over total length. Only reads on target (`config[reference][design_bed]`) |
| Median | [Mosdepth](https://github.com/brentp/mosdepth) | Median Coverage over coding exon in design (`config[reference][exon_bed]`) |
| >= 30X | [Mosdepth](https://github.com/brentp/mosdepth) | Fraction of coding exons (`config[reference][exon_bed]`) with coverage over 30x |
| >=50X |[Mosdepth](https://github.com/brentp/mosdepth) | Fraction of coding exons (`config[reference][exon_bed]`) with coverage over 50x |
| Bases on Target | [Picard](https://broadinstitute.github.io/picard/) HSMetrics | Bases inside the capture design (`config[reference][design_intervals]`) |
| Fold80 |[Picard](https://broadinstitute.github.io/picard/) HSMetrics | The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets (`config[reference][design_intervals]`) |
| % Dups | [Picard](https://broadinstitute.github.io/picard/) DuplicationMetrics |  |
| Mean Insert Size | [Picard](https://broadinstitute.github.io/picard/) InsertSizeMetrics |  |
| Target Bases with zero coverage [%] | [Picard](https://broadinstitute.github.io/picard/) HSMetrics | Percent target (`config[reference][design_intervals]`) bases with 0 coverage |
| % Adapter | [fastp](https://github.com/OpenGene/fastp) | |
