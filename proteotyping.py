#!/usr/bin/env python
# encoding: utf-8
# Fredrik Boulund
# 2014-08-01
# Parse blat output and find the most commonly occuring hits

from __future__ import division
from sys import argv, exit
from collections import namedtuple
from numpy import cumsum
import operator
import cPickle
import os
import taxtree
import argparse
import logging
import re

from prepare_taxdump_refseq import find_files, Annotation

def parse_commandline(argv):
    """Parse commandline arguments"""

    desc = """Proteotyping pipeline. (c) Fredrik Boulund 2014."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("FILE", nargs="+",
            help="BLAT output file.")
    parser.add_argument("-d", dest="display", type=int, metavar="N",
            default=10,
            help="Number of results to display [%(default)s].")
    parser.add_argument("--taxdumpdir", dest="taxdumpdir", metavar="DIR",
            default="/shared/db/NCBI/taxonomy/taxdump",
            help="Path to NCBI Taxonomy dump folder (must contain 'names.dmp', 'nodes.dmp' [%(default)s].")
    parser.add_argument("--id_gi_accno_pickle", dest="id_gi_accno_pickle", metavar="FILE",
            default="id_gi_accno.pkl",
            help="Filename of ID_GI_ACCNO pickle [%(default)s].")
    parser.add_argument("--accno_annotation_pickle", dest="accno_annotation_pickle", metavar="FILE",
            default="accno_annotation.pkl",
            help="Filename of ID_GI_ACCNO pickle [%(default)s].")
    parser.add_argument("--taxtree_pickle",  dest="taxtree_pickle", metavar="FILE",
            default="taxtree.pkl",
            help="Filename of pickled previously constructed taxtree to load instead of making it from scratch [%(default)s].")
    parser.add_argument("--gene_info", dest="gene_info_file", metavar="FILE",
            default="/shared/db/NCBI/gene/gene_info",
            help="NCBI 'gene_info' file [%(default)s].")
    parser.add_argument("--rebase_tree", dest="rebase_tree", metavar="S", type=str,
            default="2", 
            help="Rebase the taxonomic tree to this node (taxid). Only applicable when creating tree from scratch [%(default)s].")
    parser.add_argument("--taxonomic_rank", dest="taxonomic_rank", metavar="LVL", type=str,
            choices=["no rank", "subspecies", "species", "genus", "family", "order", "class", "phylum", "superkingdom"],
            default="species",
            help="Set the taxonomic level on which hits are grouped [%(default)s].")
    parser.add_argument("--fragment_length", dest="fragment_length", metavar="L", type=int,
            default=6,
            help="Minimum fragment length [%(default)s].")
    parser.add_argument("-i", "--identity", dest="identity", metavar="I", type=float,
            default=90,
            help="Filter out hits with less than or equal to this percentage identity [%(default)s].")
    parser.add_argument("--fragment_coverage", dest="fragment_coverage", metavar="C", type=float,
            default=1,
            help="Amount of fragment covered in alignment [%(default)s].")
    parser.add_argument("-m", "--min_matches", dest="min_matches", metavar="m", type=int,
            default=15, #0
            help="Filter out hits with less than this  number of matches [%(default)s].")
    parser.add_argument("-M", "--max_mismatches", dest="max_mismatches", metavar="M", type=int,
            default=0, #100, #sys.maxint,
            help="Filter out hits with more than or equal to this number of mismatches [%(default)s].")
    parser.add_argument("--remove_noninformative", dest="remove_noninformative", action="store_false",
            default=True, # TODO: This is unintuitive: reverse wording
            help="Remove fragments that match to more than one entity at the specific taxonomic level [%(default)s].")
    parser.add_argument("--maxprint", dest="maxprint", metavar="N", type=int,
            default=50,
            help="Maximum number of hit annotations to print [%(default)s].")
    parser.add_argument("--print_all_hit_annotations", dest="print_all_hit_annotations", action="store_true",
            default=False,
            help="Print a listing of all hit annotations (and not only hits from filtered peptide fragments) [%(default)s].")


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
    devoptions.add_argument("--walk", dest="walk", action="store_true",
            default=False,
            help="Instead of only visiting 'leaf' nodes when printing, walk the entire distance from the root node down [%(default)s].")

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

    Hit = namedtuple("Hit", ["target_accno", 
        "identity", "matches", "mismatches", "score", "startpos", "endpos", "fragment_length", "fragment_coverage"])
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
                startpos = int(blast8_line[8])
                endpos = int(blast8_line[9])
                fragment_length = int(fragment_id.split("_")[1])
                fragment_coverage = matches/fragment_length
                hit = Hit(target_accno, identity, matches, mismatches, 0, startpos, endpos, fragment_length, fragment_coverage)
                try:
                    hits[fragment_id].append(hit)
                except KeyError:
                    hits[fragment_id] = [hit]
                hit_counter += 1
    elif file_format == "psl":
        logging.warning("BLAT psl format not completely supported. Identity not computed.") # TODO: compute identity
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
                hit = Hit(target_accno, -1, matches, mismatches, 0, 0, 0) # TODO: compute identity, startpos, endpos
                try:
                    hits[fragment_id].append(hit)
                except KeyError:
                    hits[fragment_id] = [hit]
    else:
        logging.error("Don't know what to do with file format '{}'.".format(file_format))
        exit()
    logging.info("Parsed {} hits for {} fragments.".format(hit_counter, len(hits)))
    return hits



def filter_hits(hits, options): #remove_noninformative, identity, best_hits_only, matches, mismatches):
    """Filter hits based on user critera.
    """

    filtered_hits = {}
    for fragment_id, hitlist in hits.iteritems():
        startlen = len(hitlist)

        # Remove non-informative fragments
        if options.remove_noninformative:
            if not startlen == 1:
                logging.debug("Fragment '{}' had {} hits and was removed.".format(fragment_id, startlen))
                continue

        # Filter hits based on user critera
        filtered_hitlist = []
        for hit in hitlist:
            logging.debug("About to filter {}.".format(hit))
            if hit.identity >= options.identity:
                logging.debug("  Passed identity.")
                if hit.fragment_length >= options.fragment_length:
                    logging.debug("  Passed fragment length.")
                    if hit.fragment_coverage >= options.fragment_coverage:
                        logging.debug("  Passed fragment coverage.")
                        filtered_hitlist.append(hit)
            #if hit.identity >= options.identity and hit.matches >= matches and hit.mismatches <= mismatches:
                #filtered_hitlist.append(hit)
        if len(filtered_hitlist)>0:
            filtered_hitlist.sort(key=lambda hit: hit.matches, reverse=True) # Just nice to have it sorted
            filtered_hits[fragment_id] = filtered_hitlist

    logging.info("Filtered {} fragments. {} fragments with {} hits remain.".format(len(hits)-len(filtered_hits), 
        len(filtered_hits), sum(map(len, [hits for hits in filtered_hits.itervalues()]))))


    return filtered_hits



def insert_hits_into_tree(tree, hits):
    """Updates counts on nodes in-place.
    """
    accnos_set = set()
    count_per_accno = {}
    for fragment_id, hitlist in hits.iteritems():
        for hit in hitlist:
            accnos_set.add(hit.target_accno)
            try:
                count_per_accno[hit.target_accno] += 1
            except KeyError:
                count_per_accno[hit.target_accno] = 1
    for node in tree.traverse():
        for accno in node.accno:
            if accno in accnos_set:
                node.count += count_per_accno[accno] 
                logging.debug("Updated count of node {} with accno {} to {}.".format(node.name, accno, node.count)) # Verbose


def sum_up_the_tree(tree):
    """Sums counts up the taxonomic tree, all the way to the root node.

    In the end, root node.count should be the same as the total number of
    counts in the tree when starting out.
    """
    for leaf in tree.get_leaves(): #(is_leaf_fn=lambda n: n.count > 0):
        lineage = leaf.get_ancestors()
        lineage.insert(0, leaf)
        counts_per_node = [n.count for n in lineage]
        cumulative_counts = list(cumsum(counts_per_node))
        for idx, parent in enumerate(lineage[1:]):
            parent.count += leaf.count


def print_top_n_hits(tree, taxonomic_rank, totalhits, n=10, walk=False):
    """Prints the top 'n' nodes with the highest counts (and thus percentages).
    """
    nodes_to_print = []
    if walk:
        for node in tree.traverse(is_leaf_fn=lambda n: n.rank == taxonomic_rank):
            if node.count > 0:
                nodes_to_print.append(node)
    else:
        for node in tree.iter_leaves(is_leaf_fn=lambda n: n.rank == taxonomic_rank and n.count > 0):
            nodes_to_print.append(node)
    nodes_to_print.sort(key=lambda n: n.count, reverse=True)
    num_nodes = len(nodes_to_print)
    if num_nodes > 0:
        if num_nodes < n:
            N = num_nodes
        else:
            N = n
        print "     %  #         TAXID       TAXNAME                         ACCNO(s)"  
        for i in xrange(0,N):
            n = nodes_to_print[i]
            print "{:>6.2f}  {:<8}  {:<10}  {:<30}  {:<}".format(100*n.count/float(totalhits), n.count, n.name, n.taxname, ";".join(n.accno))
    else:
        print "Nothing filtered through at {} level.".format(taxonomic_rank)


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


def load_gene_info(gene_info_file):
    """Reads NCBI gene_info a dictionary for geneID lookups.
    """
    gene_info = {}
    with open(gene_info_file) as f:
        if not f.readline().startswith("#Format:"):
            logging.error("'{}' does not appear to be an NCBI gene_info file.".format(gene_info_file))
            raise IOError("'{}' does not appear to be an NCBI gene_info file.".format(gene_info_file))
        for line in f:
            info = line.split("\t")
            taxid = info[0]
            geneID = info[1]
            symbol = info[2]
            desc = info[8]
            gene_info[geneID] = (symbol, desc, taxid)
    return gene_info



def count_annotation_hits(hits, annotation):
    """Annotates hits to determine what genes were hit.
    """

    missing_annotation = set()

    gene_counts = {}
    for hitlist in hits.itervalues():
        for hit in hitlist:
            accno = hit.target_accno
            try:
                for annot in annotation[accno]:
                    if hit.startpos >= annot.startpos and hit.startpos <= annot.endpos:
                        if hit.endpos <= annot.endpos and hit.endpos >= annot.startpos:
                            try:
                                gene_counts[annot.geneID] += 1
                            except KeyError:
                                gene_counts[annot.geneID] = 1
            except KeyError:
                missing_annotation.add(accno)
                #logging.warning("Found hits to {} but can't find annotation.".format(accno))

    if missing_annotation:
        logging.warning("Couldn't find annotation for:\n{}.".format("\n".join(missing_annotation)))
    return gene_counts


def print_gene_counts(gene_counts, gene_info, maxprint=50, sort=True):
    """Prints counts of annotated regions with some gene info.
    """
    if sort:
        counts = sorted(gene_counts.items(), key=operator.itemgetter(1), reverse=sort)
    printcounter = 0
    print "geneID    #     Symbol             Description" 
    for geneID, count in counts:
        if printcounter < maxprint:
            try:
                symbol, desc = gene_info[geneID][0:2]
            except KeyError:
                logging.warning("Couldn't find information for geneID {}.".format(geneID))
            print "{:<9} {:<5} {:<18} {:<}".format(geneID, count, symbol, desc)
            printcounter += 1
        else:
            break
    print "Printed {} out of {} hit annotated regions.".format(printcounter, len(counts))



def load_annotation(annotation_pickle):
    """Loads annotation from previously prepared gff info pickle.
    """

    with open(annotation_pickle, 'rb') as pkl:
        annotation = cPickle.load(pkl)
    return annotation



def main(filename, options):
    """Main function that runs the complete pipeline logic.
    """
    hits = parse_blat_output(filename, options.file_format)

    filtered_hits = filter_hits(hits, options)
    if logging.getLogger().getEffectiveLevel < 20:
        for fragment, hitlist in filtered_hits.iteritems():
            logging.debug("All filtered hits:\n{}".format(hitlist))
    totalhits = sum(map(len, [hits for hits in filtered_hits.itervalues()]))

    insert_hits_into_tree(tree, filtered_hits)
    sum_up_the_tree(tree)

    print "-"*68
    print "Results at {} level for file {}.".format(options.taxonomic_rank, filename)
    print_top_n_hits(tree, options.taxonomic_rank, totalhits, n=options.display, walk=options.walk)
    print " Total: {:<}".format(totalhits)

    print "-"*68
    if options.print_all_hit_annotations:
        gene_counts = count_annotation_hits(hits, options.accno_annotation_pickle)
        print "All hit annotated regions:"
        print_gene_counts(gene_counts, gene_info, options.maxprint)
    else:
        gene_counts = count_annotation_hits(filtered_hits, options.accno_annotation_pickle)
        print "Annotated regions hit by filtered fragments:"
        print_gene_counts(gene_counts, gene_info, options.maxprint)



if __name__ == "__main__":

    options = parse_commandline(argv)

    tree = load_taxtree(options.taxtree_pickle, options.taxdumpdir, 
            options.id_gi_accno_pickle, options.rebase_tree)
    gene_info = load_gene_info(options.gene_info_file)
    annotation = load_annotation(options.accno_annotation_pickle)

    for filename in options.FILE:
        main(filename, options)
