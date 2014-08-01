#!/usr/bin/env python
# Fredrik Boulund
# 2014-07-25
# Recursively read *.fna files in directories and parse
# FASTA headers. Create a dictionary object with the 
# information obtained and pickle it to disk for later use.

from sys import argv, exit
import fnmatch
import os

#if len(argv) <2:
#    print "usage: script.py something"
#    exit()


def parse_fasta(filename):
    """Reads information form FASTA headers."""

    def trim_header(h):
        return h[1:].strip().split(" ", 1)


    output = []
    with open(filename) as f:
        line = f.readline()
        if not line.startswith(">"):
            raise Exception("File does not appear to be valid FASTA\n  "+filename)
        else:
            output.append(trim_header(line))
        for line in f:
            if line.startswith(">"):
                output.append(trim_header(line))

    return output

def find_files(directory, pattern):
    """Generator that yields files by recursively searching a dir with a glob pattern."""
    for root, subfolders, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


pickled_file = "/shared/genomes/NCBI/bacterial/20140228_pydict.pickle"
genomesdir = "/shared/genomes/NCBI/bacterial/20140228"
bigdict = {}

with open(pickled_file, "w") as fout:
    for fasta_file in find_files(genomesdir, "*.fna"):
        try:
            id_no, other = parse_fasta(fasta_file)[0]
            for seqinfo in parse_fasta(fasta_file):
                id_no, other = seqinfo
                bigdict[id_no] = other
        except Exception, e:
            print e.message

import pickle
pickle.dump(bigdict, open(pickled_file, 'wb'))
