#!/usr/bin/env python
# Fredrik Boulund
# 2014-08-01
# First test with ruffus

from sys import argv, exit
from subprocess import Popen, PIPE
import shlex
import os
import argparse
from ruffus import *



def parse_commandline(argv):
    """Parses commandline arguments"""
    parser = argparse.ArgumentParser(description="Proteotyping pipeline. Fredrik Boulund 2014.")
    parser.add_argument("-i", "--input", required=True,
            help="The peptide fragment FASTA file to search with")
    parser.add_argument("-r", "--referencepath", 
            help="Path to dir with FASTA files of reference genomes",
            default="/shared/genomes/NCBI/bacterial/20140228_panindex/fasta/")

    if len(argv)<2:
        parser.print_help()
        exit()
    return parser.parse_args()


options = parse_commandline(argv)

reference_genomes = []
for refdir,_,filenames in os.walk(options.referencepath):
    reference_genomes.extend([os.path.join(refdir,f) for f in filenames])
    break

print reference_genomes


# Stage 1: Map peptide fragments to reference genome
@transform(reference_genomes, 
           suffix(".fasta"), 
           ".psl",
           options.input)
def map_aa_fragments(input_file, output_file, query_file):
    """Maps aa fragments in FASTA format to a reference genome"""

    blat_call = "blat {database} {query} {output} -t=dnax -q=prot".format(database=input_file,
                                                          query=query_file,
                                                          output=output_file)
    blat = Popen(shlex.split(blat_call))



# Stage 2: Read mapping results
@transform(map_aa_fragments, 
           suffix(".psl"), 
           ".txt")
def read_mapping_results(input_file, output_file):
    """Read mapping results"""
    with open(output_file, "w") as f:
        f.write("Test\n")




pipeline_run()
