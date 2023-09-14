# Pipeline specific rules/softwares used in Marple
Rules created specifically for Marple pipeline are listed here.

## add_ref_to_vcf.smk
A pythonscript to add a line to vcf-files with reference path for Alissa to know which genome build used. The `##reference=`-line need to contain either hg38 or GRCh38 for Alissa to understand that the reference is not hg19.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__add_ref_to_vcf__add_ref_to_vcf#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__add_ref_to_vcf__add_ref_to_vcf#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__add_ref_to_vcf#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__add_ref_to_vcf#

---

## exomedepth_export.smk
A Rscript to create output files from exomedepth results. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__exomedepth_export__exomedepth_export#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__exomedepth_export__exomedepth_export#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__exomedepth_export#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__exomedepth_export#

---

## export_qc.smk
Rules that creates a `.xlsx` file per sample with aggregated coverage information.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__export_qc__export_qc_bedtools_intersect#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__export_qc__export_qc_bedtools_intersect#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__export_qc_bedtools_intersect#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__export_qc_bedtools_intersect#


### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__export_qc__export_qc_xlsx_report#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__export_qc__export_qc_xlsx_report#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__export_qc_xlsx_report#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__export_qc_xlsx_report#

---

## sample_order_multiqc.smk
A python script to create sample_replacement and sample_order files to be used in MultiQC to order samples based on order in SampleSheet.csv 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__sample_order_multiqc__sample_order_multiqc#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__sample_order_multiqc__sample_order_multiqc#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__sample_order_multiqc#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__sample_order_multiqc#

---
