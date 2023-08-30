#!/bin/python3
import sys
from pysam import VariantFile

vcf_in = VariantFile(snakemake.input.vcf)  # dosen't matter if bgziped or not. Automatically recognizes

# Add reference_description descriptions to new header
new_header = vcf_in.header
# new_header.add_line("reference="+ sys.argv[2])
new_header.add_line("##reference=" + snakemake.input.ref)

# start new vcf with the new_header
vcf_out = VariantFile(snakemake.output.vcf, 'w', header=new_header)


for record in vcf_in.fetch():
    vcf_out.write(record)


# import pdb; pdb.set_trace()