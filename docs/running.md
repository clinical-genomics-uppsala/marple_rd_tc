# Running the pipeline
## :exclamation: Requirements

### Recommended hardware 
blubb

### Software

- [python](https://www.python.org/), version 3.9 or newer
- [pip3](https://pypi.org/project/pip/)
- [virtuelenv](https://docs.python.org/3/library/venv.html)
- [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)

If this is available installation is easily executed with pip.

## :computer: Installation
A list of releases for Marple - Inherited Cancer pipeline can be found on [github](https://github.com/clinical-genomics-uppsala/marple/releases).

### Clone the Marple git repo
To clone the repository and fetch the pipeline.
```bash
# Set up a working directory path
WORKING_DIRECTORY="/path_working_to_directory"

# Set version
VERSION="v0.1.0"

# Clone selected version
git clone --branch ${VERSION} https://github.com/clinical-genomics-uppsala/marple.git ${WORKING_DIRECTORY}
```

### Create python environment
To run the Marple pipeline a python virtual environment is needed. Create a virtual environment and then install pipeline requirements specified in `requirements.txt`.
```bash
# Create a new virtual environment
python3 -m venv ${WORKING_DIRECTORY}/virtual/environment

# Enter working directory
cd ${WORKING_DIRECTORY}

# Activate python environment
source virtual/environment/bin/activate

# Install requirements
pip install -r requirements.txt
```
This will install all required softwares needed to run the pipeline in an virtual environment which you will have to remember to activate before running the pipeline each time. 

## :books: Input files 
Four different files need to be adapted to your compute environment and sequence data, `samples.tsv`, `units.tsv`, `config.yaml` and `resources.yaml`.
### Samples and units
A `samples.tsv` and a `units.tsv` file which contain sample information, sequencing meta information and location of fastq-files are needed. Specification of what columns are needed can be found at [Marple schemas](https://github.com/clinical-genomics-uppsala/marple/tree/develop/workflow/schemas). The `.tsv`-files can also be generated with the help of `hydra-genetics create-input-files`.

``` bash
hydra-genetics create-input-files -d /path/to/fastq/ --tc 1.0 -a 'adaptersequence1,adaptersequence2' --sample-type N  --sample-regex "^([0-9]{4})_S"

# If the fastq file is missing barcode in the header a default barcode can be set with -b NNNNNNNN
```


### Config
The bare-bone version of the config file can be found in the `config/config.yaml`. This need to be adapted to match the local paths to referencefiles, bedfiles, caches etc on your system.


/// details | Expand to view current config.yaml
```yaml
{% include "includes/config.yaml" %}
```
///

### Resources
An `resources.yaml` file can also be found in the `config/`-folder. This is adapted to the Uppsala Clinical Genomics' compute cluster but can be used as an indication of resources needed for the different programs. 

/// details | Expand to view current resources.yaml
```yaml
{% include "includes/resources.yaml" %}
```
///


## :rocket: Run command 
```bash
# Activate the virtual environment
source virtual/environment/bin/activate

# Run snakemake command with the extra config parameter called sequenceid
snakemake --profile snakemakeprofile --configfile config.yaml --config sequenceid="230202-test" -s /path/to/marple/workflow/Snakefile

```

