# Result files
Marple produces a lot of intermediate and result files but only files defined in [`output_files.yaml`](https://github.com/clinical-genomics-uppsala/marple/blob/develop/config/output_files.yaml) are kept, the rest are temporary and will be deleted when not needed anymore. If other files are wanted you need to edit `output_files.yaml` or add `--no-temp` to the running command.

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


### MultiQC report
General stats