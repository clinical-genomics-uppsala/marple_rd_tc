#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

# !!! Make a mei_list.txt that works for the cluster where it is used (full_path needed) !!! #

# Clone git
git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO}
cd ${PIPELINE_NAME}

# Create and activate conda envrionmnet in the current directory, then install pipeline requirements
conda create --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env python=${PYTHON_VERSION} -y
conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
conda install -c conda-forge pip -y

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi

mkdir ${PIPELINE_NAME}_${TAG_OR_BRANCH}
git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO} ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}
./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env/bin/pip3 install -r ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/requirements.txt 
conda pack --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env -o ${PIPELINE_NAME}_${TAG_OR_BRANCH}/env.tar.gz

# Clone snakemake-wrappers and hydra-genetics
mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics

git clone https://github.com/snakemake/snakemake-wrappers.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/snakemake-wrappers

git clone https://github.com/hydra-genetics/alignment.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/annotation.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/cnv_sv.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/compression.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/compression
git clone https://github.com/hydra-genetics/filtering.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/prealignment.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/qc.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/snv_indels.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/snv_indels

# Download containers
conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
hydra-genetics prepare-environment create-singularity-files -c config/config.yaml -o apptainer_cache

cp /projects/wp3/Software/MELTv2.2.2/MELT_v2.2.2.sif apptainer_cache

#cp config/config.miarka.yaml config/config.yaml.copy
#hydra-genetics prepare-environment container-path-update -c config/config.yaml.copy -n config/config.miarka.yaml -p apptainer_cache

# Pack all cloned repositories
tar -zcvf ${PIPELINE_NAME}_${TAG_OR_BRANCH}.tar.gz ${PIPELINE_NAME}_${TAG_OR_BRANCH}

# Download references
for reference_config in "$@"
do
    hydra-genetics --debug references download -o design_and_ref_files -v $reference_config
done

conda deactivate

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
fi

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi

# Compress data
tar -czvf design_and_ref_files.tar.gz design_and_ref_files
