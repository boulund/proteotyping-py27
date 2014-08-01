#!/usr/bin/env python
# Fredrik Boulund
# 2014-08-01
# Parse blat output and find the most commonly occuring 
# hits

from sys import argv, exit
from collections import Counter
import os

if len(argv) <2:
    print "usage: script.py output.psl N"
    exit()
else:
    blat_output = argv[1]
    if len(argv)==3:
        most_common_count = int(argv[2])
    else:
        most_common_count = 10


def create_id_species_dict(filename):
    """Creates a dictionary of NC_001122 to specicies name mappings"""
    with open(filename) as f:
        d = {}
        for line in f:
            identifier, name = line.split("\t")
            d[identifier] = name.strip()
        return d

def parse_id(s):
    """Parses the NC_001122 id from a regular gi string"""
    return s.split("|")[-2][0:-2]

def parse_gi_hits(filename):
    """Generator that yields gi hits by parsing blat output"""
    with open(filename) as f:
        # Reads first five lines before actual mapping output
        line = f.readline()
        if not line.startswith("psLayout"):
            print "File not blat output format?"
            exit()
        else:
            f.readline()
            f.readline()
            f.readline()
            f.readline()

        for line in f:
            yield line.split()[13]

def generate_words(gis, id_species):
    """Produce a list of all bacterial words in the most common hits"""
    wordlist = []
    for gi in gis:
        gi = gi[0]
        wordlist.extend(id_species[parse_id(gi)].split())
    return wordlist

genomes_map = "/shared/genomes/NCBI/bacterial/20140228_panindex/bacterial_genomes_map.txt"
id_species = create_id_species_dict(genomes_map)

c = Counter(parse_gi_hits(blat_output))
words = Counter(generate_words(c.most_common(most_common_count), id_species))

for gi, count in c.most_common(most_common_count):
    print "{:>4}: {:<32} {}".format(count, gi, id_species[parse_id(gi)])


for word, count in words.most_common(most_common_count):
    print "{:>4}: {}".format(count, word)



