# Result files
Marple produces a lot of intermediate and result files but only files defined in [`output_files.yaml`](https://github.com/clinical-genomics-uppsala/marple_rd_tc/blob/develop/config/output_files.yaml) are kept, the rest are temporary and will be deleted when not needed anymore. If other files are wanted you need to edit `output_files.yaml` or add `--no-temp` to the running command.

## Files
The following output files are located in `Results/`-folder:

| File | Format |Description |
|---|---|---|
| `{sequenceid}_MultiQC.html` | html | Aggregated QC values for entire sequence run, open in browser |
|`{sample}_{sequenceid}/{sample}_{sequenceid}.xlsx`| xlsx | Excel file with QC stats (primarily coverage) for each sample|
|`{sample}_{sequenceid}/{sample}_{sequenceid}.bam`| bam |Deduplicated alignment file |
|`{sample}_{sequenceid}/{sample}_{sequenceid}.bam.bai`| bai | Index for alignment file|
|`{sample}_{sequenceid}/{sample}_{sequenceid}.vcf.gz`| vcf.gz | Compressed VCF-file normalized and annotated with vep |
|`{sample}_{sequenceid}/{sample}_{sequenceid}.vcf.gz.tbi`| tbi | Index for variant file |
|`{sample}_{sequenceid}/{sample}_{sequenceid}.merged.genome.vcf.gz`| genome.vcf.gz | Compressed VCF-file for all positions in the design|
|`{sample}_{sequenceid}/{sample}_{sequenceid}.merged.genome.vcf.gz.tbi`| tbi |Index for genome VCF-file |
|`{sample}_{sequenceid}/{sample}_{sequenceid}_exomedepth_SV.txt`| txt | Nexus SV text file with structural variants from Exomedepth|
|`{sample}_{sequenceid}/{sample}_{sequenceid}_exomedepth.aed`| aed | aed text file with structural variants from Exomedepth |


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
The general statistics table are ordered based on the sample order in `SampleSheet.csv`, this is done by renaming the samples in two steps using the script `sample_order_multiqc.py`. To toggle between "Sample Order" and "Sample Name" use the buttons just above General Stats.

<br />

| Column Name | Origin | Comment | 
| --- | --- | --- |
| K Total Seq | [Samtools stats](http://www.htslib.org/doc/samtools-stats.html) | Only reads on target (`config[reference][design_bed]`) |
| K Reads Mapped| [Samtools stats](http://www.htslib.org/doc/samtools-stats.html) | Only reads on target (`config[reference][design_bed]`) |
| % Mapped| [Samtools stats](http://www.htslib.org/doc/samtools-stats.html) | Only reads on target (`config[reference][design_bed]`) |
| % Proper pairs| [Samtools stats](http://www.htslib.org/doc/samtools-stats.html) | Only reads on target (`config[reference][design_bed]`) |
| Average Quality | [Samtools stats](http://www.htslib.org/doc/samtools-stats.html) | Ratio between sum of base quality over total length. Only reads on target (`config[reference][design_bed]`) |
| Median | [Mosdepth](https://github.com/brentp/mosdepth) | Median Coverage over coding exon in design (`config[reference][exon_bed]`) |
| >= 30X | [Mosdepth](https://github.com/brentp/mosdepth) | Fraction of coding exons (`config[reference][exon_bed]`) with coverage over 30x |
| >=50X |[Mosdepth](https://github.com/brentp/mosdepth) | Fraction of coding exons (`config[reference][exon_bed]`) with coverage over 50x |
| Bases on Target | [Picard](https://broadinstitute.github.io/picard/) HSMetrics | Bases inside the coding exons in the design (`config[reference][exon_bed]`) |
| Fold80 |[Picard](https://broadinstitute.github.io/picard/) HSMetrics | The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets (`config[reference][exon_bed]`) |
| % Dups | [Picard](https://broadinstitute.github.io/picard/) DuplicationMetrics |  |
| Mean Insert Size | [Picard](https://broadinstitute.github.io/picard/) InsertSizeMetrics |  |
| Target Bases with zero coverage [%] | [Picard](https://broadinstitute.github.io/picard/) HSMetrics | Percent target (`config[reference][exon_bed]`) bases with 0 coverage |
| K Reads | [Picard](https://broadinstitute.github.io/picard/) HSMetrics | Total number of reads in inputfile (`alignment/samtools_merge_bam/{sample}_{type}.bam`) |
| % Adapter | [fastp](https://github.com/OpenGene/fastp) | |
