# Pipeline harmonization notes

This repository is being aligned with the CGU Hydra pipeline layout described in
the bootstrap/harmonization plan.

## Added in this pass

- `config/site/example.local.yaml`
- `config/site/miarka.local.yaml`
- `config/site/marvin.local.yaml`
- `profiles/slurm/config.yaml`
- `profiles/local/config.yaml`
- `containers/manifest.yaml`

The existing `config/config.miarka.yaml`, `profiles/config.yaml`, and
`profiles/miarka/config.yaml` remain in place for compatibility while start
scripts and cluster validation catch up.

## Current standard invocation target

```bash
snakemake \
  --profile profiles/slurm \
  --configfile config/config.yaml \
  --configfile config/site/marvin.local.yaml
```

For Miarka:

```bash
snakemake \
  --profile profiles/miarka \
  --configfile config/config.yaml \
  --configfile config/site/miarka.local.yaml
```

## Remaining refactor candidates

- Move remaining absolute reference paths out of `config/config.yaml`.
- Convert container references in `config/config.yaml` to `{{CONTAINER_DIR}}/*.sif`
  after cluster validation of the manifest.
- Rename `workflow/rules/export_qc.smk` toward the shared `export.smk` naming
  convention once callers and docs are updated.
- Compare `workflow/scripts/export_qc_xlsx_report.py`,
  `workflow/scripts/sample_order_multiqc.py`, and VCF/BED helper scripts with
  Hastings/Poirot for possible migration to `hydra-genetics/report`.
- Remove `github(...)` usage from `workflow/Snakefile_references` before treating
  the references workflow as offline-ready.
