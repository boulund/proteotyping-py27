#!/usr/bin/env python
# encoding: utf-8
# Fredrik Boulund
# 2014-08-01
# Parse blat output and find the most commonly occuring hits

from __future__ import division
from sys import argv, exit
from collections import namedtuple
import cPickle
import os
import taxtree
import argparse
import logging
import re


def parse_commandline(argv):
    """Parse commandline arguments"""

    desc = """Proteotyping pipeline. (c) Fredrik Boulund 2014."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("FILE", nargs="*",
            help="BLAT output file.")

    parser.add_argument("-d", dest="display", type=int, metavar="N",
            default=10,
            help="Number of results to display [%(default)s].")
    parser.add_argument("--taxdumpdir", dest="taxdumpdir", metavar="DIR",
            default="/shared/genomes/NCBI/taxonomy/taxdump/",
            help="Path to NCBI Taxonomy dump folder (must contain 'names.dmp', 'nodes.dmp' [%(default)s].")
    parser.add_argument("--pickle", dest="id_gi_accno_pickle", metavar="FILE",
            default="id_gi_accno.pkl",
            help="Filename of ID_GI_ACCNO pickle [%(default)s].")
    parser.add_argument("--taxtree_pickle",  dest="taxtree_pickle", metavar="FILE",
            default="taxtree.pkl",
            help="Filename of pickled previously constructed taxtree to load instead of making it from scratch [%(default)s].")
    parser.add_argument("--rebase_tree", dest="rebase_tree", metavar="S", type=str,
            default="2", 
            help="Rebase the taxonomic tree to this node (taxid). Only applicable when creating tree from scratch [%(default)s].")
    parser.add_argument("--taxonomic_level", dest="taxonomic_level", metavar="LVL",
            choices=["subspecies", "species", "genus", "family"], #, "order", "class", "phylum", "superkingdom"],
            default="species",
            help="Set the taxonomic level on which hits are grouped [%(default)s].")
    parser.add_argument("-m", "--min_matches", dest="min_matches", metavar="m", type=int,
            default=15, #0
            help="Filter out hits with less than this  number of matches [%(default)s].")
    parser.add_argument("-M", "--max_mismatches", dest="max_mismatches", metavar="M", type=int,
            default=0, #100, #sys.maxint,
            help="Filter out hits with more than or equal to this number of mismatches [%(default)s].")
    parser.add_argument("-i", "--identity", dest="identity", metavar="I", type=float,
            default=90,
            help="Filter out hits with less than or equal to this percentage identity [%(default)s].")
    parser.add_argument("--remove_noninformative", dest="remove_noninformative", action="store_false",
            default=True, # TODO: This is unintuitive: reverse wording
            help="Remove fragments that match to more than one entity at the specific taxonomic level [%(default)s].")


    devoptions = parser.add_argument_group("Developer options", "Voids warranty ;)")
    devoptions.add_argument("--loglevel", choices=["INFO", "DEBUG"],
            default="INFO", 
            help="Set logging level [%(default)s].")
    devoptions.add_argument("--best_hits_only", dest="best_hits_only", action="store_true",
            default=False,
            help="For each fragment, remove hits that have less than the maximum number of matches and more than the minimum number of mismatches. [%(default)s].")
    devoptions.add_argument("--file_format", dest="file_format", 
            choices=["blast8","psl"],
            default="blast8",
            help="Specify format of mappings file from blat [%(default)s].")

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
    ref_id_regex = re.compile(r"\|ref\|((\w{1,2}_[\d\w]+)\.\d{1,2})\|")
    hit = re.search(ref_id_regex, s)
    if hit is not None:
        return hit.group(group)
    else:
        raise Exception("ERROR finding accession string in {}".format(s))



def parse_blat_output(filename, file_format="blast8"):
    """Parses blat output. Returns all hits for each fragment.

    A hit is a named tuple with the following fields:
    target_accno, identity, mathes, mismatches, score.
    """

    Hit = namedtuple("Hit", ["target_accno", "identity", "matches", "mismatches", "score"])
    hit_counter = 0
    hits = {}

    if file_format == "blast8":
        with open(filename) as blast8:
            counter = 0
            for line in blast8:
                blast8_line = line.split()
                fragment_id = blast8_line[0]
                target_accno = parse_accno(blast8_line[1])
                identity = float(blast8_line[2])
                matches = int(blast8_line[3])
                mismatches = int(blast8_line[4])
                hit = Hit(target_accno, identity, matches, mismatches, 0)
                try:
                    hits[fragment_id].append(hit)
                except KeyError:
                    hits[fragment_id] = [hit]
                hit_counter += 1
    elif file_format == "psl":
        with open(filename) as f:
            line = f.readline()
            if not line.startswith("psLayout"):
                logging.error("File '{}' not blat output format?".format(filename))
                exit()
            else:
                [f.readline() for x in range(0,4)] # Skip past column headers

            hits = {}
            for line in f:
                split_line = line.split()
                matches, mismatches, repmatches = [int(info) for info in split_line[0:3]]
                ncount, qinserts, qbaseinserts, tinserts, tbaseinserts = split_line[3:8]
                fragment_id = split_line[9]
                target_accno = parse_accno(split_line[13], group=1)
                hit = Hit(target_accno, -1, matches, mismatches, 0) # TODO: compute identity
                try:
                    hits[fragment_id].append(hit)
                except KeyError:
                    hits[fragment_id] = [hit]
    else:
        logging.error("Don't know what to do with file format '{}'.".format(file_format))
        exit()
    logging.info("Parsed {} hits for {} fragments.".format(hit_counter, len(hits)))
    return hits



def filter_hits(hits, remove_noninformative, identity, best_hits_only, matches, mismatches):
    """Filter hits based on user critera.
    """

    filtered_hits = {}
    for fragment_id, hitlist in hits.iteritems():
        startlen = len(hitlist)

        # Remove non-informative fragments
        if remove_noninformative:
            if startlen == 1:
                pass
            else:
                logging.debug("Fragment '{}' had {} hits and was removed.".format(fragment_id, startlen))
                continue

        # Filter hits based on user critera
        hitlist.sort(key=lambda hit: hit.matches, reverse=True)
        hitlist = [hit for hit in hitlist if hit.matches >= matches]
        filtered_hitlist = []
        for hit in hitlist:
            if hit.identity >= identity and hit.matches >= matches and hit.mismatches <= mismatches:
                filtered_hitlist.append(hit)
        if len(filtered_hitlist)>0:
            filtered_hits[fragment_id] = filtered_hitlist

    logging.info("Filtered {} fragments. {} fragments remain.".format(len(hits)-len(filtered_hits), len(filtered_hits)))

    return filtered_hits



def insert_hits_into_tree(tree, hits):
    """Updates counts on nodes in-place.
    """
    accnos_set = set()
    accno_fragment_id = {}
    for fragment_id, hitlist in hits.iteritems():
        for hit in hitlist:
            accnos_set.add(hit.target_accno)
            accno_fragment_id[hit.target_accno] = fragment_id
    for node in tree.traverse():
        for accno in node.accno:
            if accno in accnos_set:
                node.count += len(hits[accno_fragment_id[accno]])
                logging.debug("Updated count of node {} with accno {} to {}.".format(node.name, accno, node.count)) # Verbose


def sum_up_to_taxonomic_level(tree, taxonomic_level):
    """Sums the count of all counts below a specified level to nodes at the specified level.
    """
    for node in tree.traverse(is_leaf_fn=lambda n: n.rank == taxonomic_level):
        for child in node.iter_descendants():
            node.count += child.count
        logging.debug("Node {} now has a count of {}.".format(node.name, node.count)) # Verbose


def print_top_n_hits(tree, taxonomic_level, n=10):
    """Prints the top 'n' nodes with highest counts.
    """
    nodes_with_counts = []
    for node in tree.traverse():
        if node.rank == taxonomic_level and node.count > 0:
            nodes_with_counts.append(node)

    #print nodes_with_counts
    nodes_with_counts.sort(key=lambda n: n.count, reverse=True)
    #print nodes_with_counts

    num_nodes = len(nodes_with_counts)
    if num_nodes > 1:
        if num_nodes < n:
            n = num_nodes
        print "PERCENTAGE  #FRAGMENTS  ACCNO(S)    TAXNAME"  
        for i in xrange(0,n):
            n = nodes_with_counts[i]
            print "{:>10.3f}  {:>10}  {:<10}  {:<}".format(n.count/float(totalhits), n.count, ";".join(n.accno), n.taxname)
            #print n.taxname, n.accno, n.count, n.count/float(totalhits)
    else:
        print "Nothing filtered through at {} level.".format(taxonomic_level)


def load_taxtree(taxtree_pickle, taxdumpdir, id_gi_accno_pickle, rebase):
    """Determine if previous taxtree_pickle is available and load it, otherwise create new.
    """
    if os.path.isfile(taxtree_pickle):
        logging.debug("Found taxtree pickle '{}', loading...".format(taxtree_pickle))
        with open(taxtree_pickle, 'rb') as pickled_tree:
            tree = cPickle.load(pickled_tree)
        logging.debug("Taxtree pickle '{}' loaded.".format(taxtree_pickle))
    else:
        logging.debug("Found no pickled taxtree called '{}', creating one from scratch.".format(taxtree_pickle))
        tree = taxtree.load_ncbi_tree_from_taxdump(taxdumpdir, id_gi_accno_pickle, rebase=rebase)
        with open(taxtree_pickle, 'wb') as picklejar:
            cPickle.dump(tree, picklejar, -1)
        logging.debug("Pickled taxtree to '{}'.".format(taxtree_pickle))
    return tree




if __name__ == "__main__":
    
    options = parse_commandline(argv)

    # Load taxtree
    tree = load_taxtree(options.taxtree_pickle, 
            options.taxdumpdir, options.id_gi_accno_pickle, options.rebase_tree)

    for filename in options.FILE:
        hits = parse_blat_output(filename, options.file_format)

        filtered_hits = filter_hits(hits, 
                options.remove_noninformative,
                options.identity,
                options.min_matches,
                options.max_mismatches,
                options.best_hits_only)

        insert_hits_into_tree(tree, filtered_hits)
        sum_up_to_taxonomic_level(tree, options.taxonomic_level)

        totalhits = sum(map(len, [hits for hits in filtered_hits.itervalues()]))
        
        print "-"*68
        print "Results at {} level for file {}.".format(options.taxonomic_level, filename)
        print "Total number of hits: {}.".format(totalhits)
        print_top_n_hits(tree, options.taxonomic_level, n=options.display)

