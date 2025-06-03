#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import pandas as pd
import re
import sys

if not os.path.exists("snv_indels/deepmosaic"):
    os.makedirs("snv_indels/deepmosaic")

bam_path = snakemake.input.bam
vcf_path = snakemake.input.vcf
sample = snakemake.params.name

path = "{PWD}"
if re.search("TE", path):
    depth = 100
elif re.search("TC", path):
    depth = 200
else:
    depth = 40

file = snakemake.output.txt
with open(file, "w") as output:
    output.write("#sample_name\tbam\tvcf\tdepth\tsex\n" +
                 str(sample) + "\t" + str(bam_path) +
                 "\t" + str(vcf_path) + "\t" + str(depth) +
                 "\t" + "female" + "\n")
