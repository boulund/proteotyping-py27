#!/usr/bin/env python
# encoding: utf-8
# Fredrik Boulund
# 2014-08-01
# Parse blat output and find the most commonly occuring hits

from __future__ import division
from sys import argv, exit
from collections import namedtuple
import taxtree
import argparse
import logging
import re


def parse_commandline(argv):
    """Parse commandline arguments"""

    desc = """Proteotyping pipeline. (c) Fredrik Boulund 2014."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("PSLFILE", help="BLAT output file.")

    parser.add_argument("-d", "--display", dest="display", type=int, default=10,
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



def parse_accno(s, group=1):
    """Parses the accession numbers (e.g. NC_001122.22) from a FASTA header string.
    
    User can chose to return either:
        group=1: ACCESSION
        group=2: ACCESSION.VERSION
    """

    # Regex finds e.g. NC_001122.1 in a FASTA header string
    # such as "gi|158333233|ref|NC_009925.1|"
    # hit.group(1) is whole accession number
    # hit.group(2) is accession number without VERSION
    ref_id_regex = re.compile(r"\|ref\|((\w{1,2}_\w+)\.\d{1,2})\|")
    hit = re.search(ref_id_regex, s)
    if hit is not None:
        return hit.group(group)
    else:
        raise Exception("ERROR finding accession string in {}".format(s))



def parse_blat_output(filename):
    """Parses blat output. Returns all hits for each fragment.
       A hit is a named tuple with the following fields:
       matches, mismatches, accno, score.
    """
    with open(filename) as f:
        line = f.readline()
        if not line.startswith("psLayout"):
            logging.debug("File '{}' not blat output format?".format(filename))
            exit()
        else:
            [f.readline() for x in range(0,4)] # Read past column headers

        Hit = namedtuple("Hit", ["matches", "mismatches", "accno", "score"])
        hits = {}
        for line in f:
            split_line = line.split()
            matches, mismatches, repmatches = split_line[0:3]
            ncount, qinserts, qbaseinserts, tinserts, tbaseinserts = split_line[3:8]
            fragment_id = split_line[9]
            accno = parse_accno(split_line[13], group=2)
            hit = Hit(matches, mismatches, accno, 0)
            try:
                hits[fragment_id].append(hit)
            except KeyError:
                hits[fragment_id] = [hit]

        return hits



def choose_best_hits(hits):
    """Choose the best hit(s) from a list of hits."""

    for fragment_id, hitlist in hits.iteritems():
        startlen = len(hitlist)
        # First, filter out hits with maximum number of matches
        hitlist.sort(key=lambda hit: hit.matches, reverse=True)
        hit_max_matches = max(hitlist)
        max_matches_positions = [i for i, hit in enumerate(hitlist) if hit.matches == hit_max_matches.matches]
        hitlist = [hitlist[pos] for pos in max_matches_positions]

        # Second, filter out hits with more mismatches than minimum
        hit_min_mismatches = min(hitlist, key=lambda hit: hit.mismatches)
        min_mismatches_positions = [i for i, hit in enumerate(hitlist) if hit.mismatches == hit_min_mismatches.mismatches]
        hitlist = [hitlist[pos] for pos in min_mismatches_positions]
        logging.debug("Removed {} hits for fragment {}".format(startlen-len(hitlist), fragment_id))

        hits[fragment_id] = ([hit.accno for hit in hitlist], startlen)
    return hits 



def score_hits(hits):
    """Scores the most commonly occuring hits with the reciprocial of their
    number of hits per fragment."""

    scores = {}
    for key in hits.iterkeys():
        numhits = hits[key][1]
        for accno in hits[key][0]:
            try:
                scores[accno] += 1/numhits
            except KeyError:
                scores[accno] = 1/numhits

    scores_list = list(scores.items())
    scores_list_reversed = [(pair[1], pair[0]) for pair in scores_list]
    sorted_scores_list = sorted(scores_list_reversed, 
                                key=lambda pair: pair[0], #sort by score
                                reverse=True)
    return sorted_scores_list



def insert_scores_into_tree(tree, scored_hits):
    """Walk the tree and insert the score for nodes with hits."""
    scoredict = dict([(h[1], h[0]) for h in scored_hits])
    accnos = set(scoredict.keys())
    for node in tree.traverse():
        if node.accno and node.accno[0] in accnos:
            node.score = scoredict[node.accno[0]]



# TODO: Make this accept a list of nodes instead
def print_hits(sorted_scores_list, tree, n=10):
    """Prints the top 'n' hits."""
    for i in xrange(0,n):
        score, accno = sorted_scores_list[i]
        nodes = taxtree.search_for_accno(tree, accno)
        taxname = nodes[0].taxname
        print "{:>3f}: {:<11} {:<}".format(score, accno, taxname)



if __name__ == "__main__":
    
    options = parse_commandline(argv)

    # Load taxtree
    accno_name_file =  "/shared/genomes/NCBI/bacterial/20140228_panindex/bacterial_genomes_map.txt"
    dumpdir = "/shared/genomes/NCBI/taxonomy/taxdump/"
    tree = taxtree.load_ncbi_tree_from_dump(dumpdir, accno_name_file)

    print "------------------------------ Score ranking"
    hits = parse_blat_output(options.PSLFILE)
    best_hits = choose_best_hits(hits)
    scored_hits = score_hits(best_hits)
    insert_scores_into_tree(tree, scored_hits)
    #print_hits(scored_hits, tree, n=options.display)


