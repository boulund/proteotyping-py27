#!/usr/bin/env python
# encoding: utf-8
# Fredrik Boulund
# 2014-08-01
# Parse blat output and find the most commonly occuring hits

from __future__ import division
from sys import argv, exit
from collections import Counter
import argparse
import logging
import re
import os


def parse_commandline(argv):
    """Parse commandline arguments"""

    desc = """Proteotyping pipeline. (c) Fredrik Boulund 2014."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("PSLFILE", help="BLAT output file.")

    parser.add_argument("--counts", "--rankbycounts", action="store_true",
            default=True,
            help="Rank results by naïve counts. Not recommended as it is probably severely biased towards your specific database.")
    parser.add_argument("--frequencies", "--rankbyfrequencies", action="store_true",
            default=True,
            help="Rank results by best hits and relative scoring.")
    parser.add_argument("--display", type=int, default=10,
            help="Number of results to display.")

    devoptions = parser.add_argument_group("Developer options", "Voids warranty ;)")
    devoptions.add_argument("--loglevel", choices=["INFO", "DEBUG"],
            default="INFO", help="Set logging level")


    if len(argv) < 2:
        parser.print_help()
        exit()

    options = parser.parse_args(argv[1:])

    logging.basicConfig(level=options.loglevel)

    return options


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


def generate_words(gis, id_species):
    """Produce a list of all bacterial words in the most common hits"""
    wordlist = []
    for gi in gis:
        gi = gi[0]
        wordlist.extend(id_species[parse_ref_id(gi)].split())
    return wordlist

def choose_best_hit(hitlist, fragment_id):
    """Chooses the best hit(s) from a list of hits."""

    startlen = len(hitlist)
    hitlist.sort(key=lambda hit: hit[0], reverse=True)
    # First, filter out hits with maximum number of matches
    max_matches = max(hitlist)
    max_matches_positions = [i for i, hit in enumerate(hitlist) if hit[0] == max_matches[0]]
    hitlist = [hitlist[pos] for pos in max_matches_positions]

    # Second, filter out hits with more mismatches than minimum
    min_mismatches = min(hitlist, key=lambda hit: hit[1])
    min_mismatches_positions = [i for i, hit in enumerate(hitlist) if hit[1] == min_mismatches[1]]
    hitlist = [hitlist[pos] for pos in min_mismatches_positions]
    logging.debug("Removed {} hits for fragment {}".format(startlen-len(hitlist), fragment_id))
    
    return [hit[2] for hit in hitlist]


def parse_blat_output(filename):
    """Parses blat output and returns counts for all hits"""
    with open(filename) as f:
        line = f.readline()
        if not line.startswith("psLayout"):
            logging.debug("File '{}' not blat output format?".format(filename))
            exit()
        else:
            [f.readline() for x in range(0,4)] # Read past column headers

        hits = {}
        for line in f:
            split_line = line.split()
            matches, mismatches, repmatches = split_line[0:3]
            ncount, qinserts, qbaseinserts, tinserts, tbaseinserts = split_line[3:8]
            fragment_id = split_line[9]
            target_id = parse_ref_id(split_line[13])
            try:
                hits[fragment_id].append((matches, mismatches, target_id))
            except KeyError:
                hits[fragment_id] = [(matches, mismatches, target_id)]
        
        for fragment_id, hitlist in hits.iteritems():
            hits[fragment_id] = choose_best_hit(hitlist, fragment_id)

        return hits


def count_most_common_occurences(id_species, blat_output, most_common_count):
    """Simple implementation using Counter. Probably biased towards what is 
    most common in the reference database."""
    c = Counter(gi_generator_from_blat_output(blat_output))
    words = Counter(generate_words(c.most_common(most_common_count), id_species))
    
    for gi, count in c.most_common(most_common_count):
        print "{:>4}: {:<32} {}".format(count, gi, id_species[parse_ref_id(gi)])
    
    for word, count in words.most_common(most_common_count):
        print "{:>4}: {}".format(count, word)
    

def score_most_common_occurences(id_species, blat_output, most_common_count):
    """Scores the most commonly occuring hits with the reciprocial of their
    occurence per fragment."""

    hits = parse_blat_output(blat_output)
    scores = {}
    for key in hits.iterkeys():
        numhits = len(hits[key])
        for hit in hits[key]:
            try:
                scores[hit] += 1/numhits
            except KeyError:
                scores[hit] = 1/numhits

    scores_list = list(scores.items())
    scores_list_reversed = [(pair[1], pair[0]) for pair in scores_list]
    sorted_scores_list = sorted(scores_list_reversed, 
                                key=lambda pair: pair[0], #sort by score
                                reverse=True)
    for i in xrange(0,most_common_count):
        print "{:>3f}: {:<}".format(sorted_scores_list[i][0],
                                    id_species[sorted_scores_list[i][1]])



if __name__ == "__main__":
    
    options = parse_commandline(argv)

    genomes_names = "/shared/genomes/NCBI/bacterial/20140228_panindex/bacterial_genomes_names.txt"
    id_species = create_id_species_dict(genomes_names)

    print "------------------------------ Count ranking"
    count_most_common_occurences(id_species, blat_output=options.PSLFILE, most_common_count=options.display)

    print "------------------------------ Score ranking"
    score_most_common_occurences(id_species, blat_output=options.PSLFILE, most_common_count=options.display)

