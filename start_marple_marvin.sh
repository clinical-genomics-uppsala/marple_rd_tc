#!/bin/bash

source /projects/wp3/nobackup/TwistCancer/Bin/venv_marple/bin/activate
pip install -r /projects/wp3/nobackup/TwistCancer/Bin/marple_rd_tc/requirements.txt
module load slurm-drmaa/1.1.3

## To match stackstorm outline
old_varible=$1
sequenceid=$2
outbox=$(pwd)/../../OUTBOX/

echo 'Creating outbox and scratch folders'
if [ ! -d "/scratch/wp3/TwistCancer/${sequenceid}/" ]
then
  mkdir /scratch/wp3/TwistCancer/${sequenceid}/
fi

if [ ! -d "${outbox}/${sequenceid}/" ]
then
  mkdir ${outbox}/${sequenceid}/
fi

## Collect statistics
echo "Statistics.."
export sequenceid
python3 -c 'import json; import os; from ductus.tools.utils import generate_elastic_statistics; writer=open("/projects/wp3/nobackup/Collect_statistics/"+os.environ["sequenceid"]+".json", "w"); [writer.write(json.dumps(d) + "\n") for d in generate_elastic_statistics(samplesheet="SampleSheet.csv", workpackage="wp3", tool="tc", analysis="Twist_Cancer", project="klinik", prep="Fresh")]; writer.close()'
## cp config and resources
echo "Creating sample.tsv and units.tsv and cp all .[tc]sv and .yaml to scratch" && \
hydra-genetics create-input-files -d /projects/wp3/nobackup/TwistCancer/Workarea/${sequenceid}/fastq/ \
  -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' -t N --tc 0  -b 'NNNNNNNNN+NNNNNNNNN' -f && \
## Need to remove 000 from miseq? sed -i 's/\t000000000-/\t/' units.tsv
cp *.[tc]sv /scratch/wp3/TwistCancer/${sequenceid}/ && \
cp /projects/wp3/nobackup/TwistCancer/Bin/marple_rd_tc/config/config.yaml . && \
cp /projects/wp3/nobackup/TwistCancer/Bin/marple_rd_tc/config/resources.yaml . && \
cp *.yaml /scratch/wp3/TwistCancer/${sequenceid}/ && \

## Run pipeline on scratch
echo "cd to scratch and run snakemake" && \
cd /scratch/wp3/TwistCancer/${sequenceid}/ && \
snakemake --profile /projects/wp3/nobackup/TwistCancer/Bin/marple_rd_tc/profiles/ --config sequenceid="${sequenceid}" && \

## Cp Results to OUTBOX and trigger stackstorm
echo 'Pipeline done, cp back data and trigger stackstorm' && \
rsync -ru /scratch/wp3/TwistCancer/${sequencerun}/Results/ ${outbox}/${sequencerun}/ && \
touch ${outbox}/${sequencerun}/Done.txt && \
echo 'All done!'