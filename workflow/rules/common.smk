__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2023, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

import json
import pathlib
import pandas
import yaml
from datetime import datetime

from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.misc import replace_dict_variables
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics.utils.misc import extract_chr
from snakemake.utils import min_version
from snakemake.utils import validate

from hydra_genetics import min_version as hydra_min_version

from hydra_genetics.utils.misc import export_config_as_file
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import use_container

min_version("7.13.0")

hydra_min_version("3.0.0")

## Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")

config = replace_dict_variables(config)

try:
    validate(config, schema="../schemas/config.schema.yaml")
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


## Load and validate resources
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pandas.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)

validate(units, schema="../schemas/units.schema.yaml")


### Read and validate output file

with open(config["output"]) as output:
    if config["output"].endswith("json"):
        output_spec = json.load(output)
    elif config["output"].endswith("yaml") or config["output"].endswith("yml"):
        output_spec = yaml.safe_load(output.read())

validate(output_spec, schema="../schemas/output_files.schema.yaml")


# Check that fastq files actually exist. If not, this might result in other
# errors that can be hard to interpret
for fq1, fq2 in zip(units["fastq1"].values, units["fastq2"].values):
    if not pathlib.Path(fq1).exists():
        sys.exit(f"fastq file not found: {fq1}\ncontrol the paths in {config['units']}")
    if not pathlib.Path(fq2).exists():
        sys.exit(f"fastq file not found: {fq2}\ncontrol the paths in {config['units']}")


### get bam input for compression_samtools_view
def get_bam_input(wildcards, use_sample_wildcard=True, use_type_wildcard=True, by_chr=False):
    if use_sample_wildcard and use_type_wildcard is True:
        sample_str = "{}_{}".format(wildcards.sample, wildcards.type)
    elif use_sample_wildcard and use_type_wildcard is not True:
        sample_str = "{}_{}".format(wildcards.sample, "N")
    else:
        sample_str = wildcards.file
 
    bam_input = "alignment/samtools_merge_bam/{}.bam".format(sample_str)
    bai_input = "{}.bai".format(bam_input)

    return (bam_input, bai_input)


### get gvcf output for deepvariant
def get_gvcf_output(wildcards, name):
    if config.get(name, {}).get("output_gvcf", False):
        return f" --output_gvcf snv_indels/deepvariant/{wildcards.sample}_{wildcards.type}_{wildcards.chr}.g.vcf.gz "
    else:
        return ""


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []

    for f in output_spec["files"]:
        # Please remember to add any additional values down below
        # that the output strings should be formatted with.
        outputpaths = set(
            [
                f["output"].format(sample=sample, type=unit_type, sequenceid=config["sequenceid"])
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
            ]
        )

        for op in outputpaths:
            output_files.append(outdir / Path(op))

    return output_files


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec.get("directory", "./"))
    rulestrings = []

    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "_copy_{}".format("_".join(re.split(r"\W{1,}", f["name"].strip().lower())))
        input_file = pathlib.Path(f["input"])
        output_file = output_directory / pathlib.Path(f["output"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code = "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")',
                "@workflow.run\n",
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
                '\tshell("(cp --preserve=timestamps {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )

        rulestrings.append(rule_code)

    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)


generate_copy_rules(output_spec)
