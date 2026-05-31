# Pipeline harmonization notes

This repository is being aligned with the CGU Hydra pipeline layout described in
the bootstrap/harmonization plan.

## Added in this pass

- `profiles/slurm/config.yaml`
- `profiles/local/config.yaml`
- `containers/manifest.yaml`
- `workflow/manifest.yaml`

The existing `config/config.miarka.yaml`, `profiles/config.yaml`, and
`profiles/miarka/config.yaml` remain in place for compatibility while start
scripts and cluster validation catch up.

## Current workflow target

The main workflow entrypoint is `workflow/Snakefile`; reference generation uses
`workflow/Snakefile_references`. Local rules, hydra modules, and script
candidates are described in `workflow/manifest.yaml`.

## Remaining refactor candidates

- Move remaining absolute reference paths behind a common `cgu/library` layout
  so site config only needs the installation parent and runtime resources.
- Convert container references in `config/config.yaml` to `{{CONTAINER_DIR}}/*.sif`
  after cluster validation of the manifest.
- Rename `workflow/rules/export_qc.smk` toward the shared `export.smk` naming
  convention once callers and docs are updated.
- Compare `workflow/scripts/export_qc_xlsx_report.py`,
  `workflow/scripts/sample_order_multiqc.py`, and VCF/BED helper scripts with
  Hastings/Poirot for possible migration to `hydra-genetics/report`.
