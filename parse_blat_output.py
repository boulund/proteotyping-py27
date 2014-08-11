#!/usr/bin/env python
# Fredrik Boulund
# 2014-08-01
# Parse blat output and find the most commonly occuring hits

from __future__ import division
from sys import argv, exit
from collections import Counter
import re
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
    """Creates a dictionary of ref id -> species name mappings"""
    with open(filename) as f:
        d = {}
        for line in f:
            identifier, name = line.split(" ", 1)
            identifier = parse_ref_id(identifier)
            d[identifier] = name.strip()
        return d


def parse_ref_id(s):
    """Parses the ref id (e.g. NC_001122.22) from a string"""

    # Regex finds e.g. NC_001122.1 in a FASTA header string
    # such as "gi|158333233|ref|NC_009925.1|"
    # hit.group(1) is whole ref id
    # hit.group(2) is base ref id
    ref_id_regex = re.compile(r"\|ref\|((\w\w_\w+)\.\d{1,2})\|")
    hit = re.search(ref_id_regex, s)
    if hit is not None:
        return hit.group(1)
    else:
        raise Exception("ERROR finding ref id string in {}".format(s))


def gi_generator_from_blat_output(filename):
    """Generator that yields gi hits by parsing blat output"""
    with open(filename) as f:
        # Reads first five lines (version info and table headers)
        # before actual mapping output
        line = f.readline()
        if not line.startswith("psLayout"):
            print "File '{}' not blat output format?".format(filename)
            exit()
        else:
            [f.readline() for x in range(0,4)] # Read past column headers

        for line in f:
            yield line.split()[13] # 13 is the target name


def parse_blat_output(filename):
    """Parses blat output and returns counts for all hits"""
    with open(filename) as f:
        line = f.readline()
        if not line.startswith("psLayout"):
            print "File '{}' not blat output format?".format(filename)
            exit()
        else:
            [f.readline() for x in range(0,4)] # Read past column headers

        hits = {}
        for line in f:
            split_line = line.split()
            fragment_id = split_line[9]
            target_id = parse_ref_id(split_line[13])
            try:
                hits[fragment_id].append(target_id)
            except KeyError:
                hits[fragment_id] = [target_id]
            
        return hits

def generate_words(gis, id_species):
    """Produce a list of all bacterial words in the most common hits"""
    wordlist = []
    for gi in gis:
        gi = gi[0]
        wordlist.extend(id_species[parse_ref_id(gi)].split())
    return wordlist


def count_most_common_occurences(id_species):
    """Simple implementation using Counter. Probably biased towards what is 
    most common in the reference database."""
    c = Counter(gi_generator_from_blat_output(blat_output))
    words = Counter(generate_words(c.most_common(most_common_count), id_species))
    
    for gi, count in c.most_common(most_common_count):
        print "{:>4}: {:<32} {}".format(count, gi, id_species[parse_ref_id(gi)])
    
    for word, count in words.most_common(most_common_count):
        print "{:>4}: {}".format(count, word)
    
def score_most_common_occurences(id_species):
    """Scores the most commonly occuring hits with the reciprocial of their
    occurence per fragment."""

    hits = parse_blat_output(blat_output)
    scores = {}
    for key in hits.keys():
        numhits = len(hits[key])
        for hit in hits[key]:
            try:
                scores[hit] += 1/numhits
            except KeyError:
                scores[hit] = 1/numhits

    scores_list = list(scores.iteritems())
    scores_list_reversed = [(pair[1], pair[0]) for pair in scores_list]
    sorted_scores_list = sorted(scores_list_reversed, 
                                key=lambda pair: pair[0], #sort by score
                                reverse=True)
    for i in xrange(0,most_common_count):
        print "{:>3f}: {:<}".format(sorted_scores_list[i][0],
                                    id_species[sorted_scores_list[i][1]])


if __name__ == "__main__":
    genomes_names = "/shared/genomes/NCBI/bacterial/20140228_panindex/bacterial_genomes_names.txt"
    id_species = create_id_species_dict(genomes_names)

    print "------------------------------ Count ranking"
    count_most_common_occurences(id_species)

    print "------------------------------ Score ranking"
    score_most_common_occurences(id_species)

