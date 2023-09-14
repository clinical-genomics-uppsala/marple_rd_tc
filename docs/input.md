### :school_satchel: Preparations


#### Samples and units
Input data should be added to a [`samples.tsv`](https://github.com/clinical-genomics-uppsala/marple_rd_tc/blob/develop/config/samples.tsv)
and an [`units.tsv`](https://github.com/clinical-genomics-uppsala/marple_rd_tc/blob/develop/config/units.tsv).

The following information need to be added to these files:

| Column Id | Description |
| --- | --- |
| **`samples.tsv`** |
| sample | unique sample/patient id, one per row |
|tumor_content| tumor cell content estimation|
| **`units.tsv`** |
| sample | same sample/patient id as in `samples.tsv` |
| type | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA |
| platform | type of sequencing platform, e.g. `NovaSeq` |
| machine | specific machine id, e.g. NovaSeq instruments have `@Axxxxx` |
| flowcell | identifier of flowcell used |
| lane | flowcell lane number |
| barcode | sequence library barcode/index, connect forward and reverse indices by `+`, e.g. `ATGC+ATGC` |
| fastq1/2 | absolute path to forward and reverse reads |
| adapter | adapter sequences to be trimmed, separated by comma |

#### Configuration
A [`resources.yaml`](https://github.com/clinical-genomics-uppsala/marple_rd_tc/blob/develop/config/resources.yaml) and [`config.yaml`](https://github.com/clinical-genomics-uppsala/marple_rd_tc/blob/develop/config/config.yaml) need to be adapted to work with your computing environment.
