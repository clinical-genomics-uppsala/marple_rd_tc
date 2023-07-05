__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2023, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

import itertools
import numpy as np
import pathlib
import pandas as pd
import yaml
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *


min_version("7.13.0")

try:
    validate(config, schema="../schemas/config_references.schema.yaml")
except WorkflowError as we:
    # Probably a validation error, but the original exsception in lost in
    # snakemake. Pull out the most relevant information instead of a potentially
    # *very* long error message.
    if not we.args[0].lower().startswith("error validating config file"):
        raise
    error_msg = "\n".join(we.args[0].splitlines()[:2])
    parent_rule_ = we.args[0].splitlines()[3].split()[-1]
    if parent_rule_ == "schema:":
        sys.exit(error_msg)
    else:
        schema_hiearachy = parent_rule_.split()[-1]
        schema_section = ".".join(re.findall(r"\['([^']+)'\]", schema_hiearachy)[1::2])
        sys.exit(f"{error_msg} in {schema_section}")


config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema="../schemas/units_references.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def compile_output_list(wildcards: snakemake.io.Wildcards):
    output_files = [
        "references/exomedepth_reference/RefCount.Rdata",
    ]
    return output_files
