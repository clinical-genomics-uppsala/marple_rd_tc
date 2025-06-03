#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import re
import sys

if not os.path.exists("snv_indels/deepmosaic"):
    os.makedirs("snv_indels/deepmosaic")

bam_path = snakemake.input.bam
depth = 200
sample = snakemake.params.name
vcf_path = snakemake.input.vcf
 
print(sample)
print(vcf_path)
print(bam_path)

file = snakemake.output.txt
with open(file, "w") as output:
    output.write("#sample_name\tbam\tvcf\tdepth\tsex\n" +
                 str(sample) + "\t" + str(bam_path) +
                 "\t" + str(vcf_path) + "\t" + str(depth) +
                 "\t" + "female" + "\n")
