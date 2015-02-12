#!/usr/bin/env python
# encoding: utf-8
# Fredrik Boulund
# 2014-08-01
# Parse blat output and find the most commonly occuring hits

from __future__ import division
from sys import argv, exit
from collections import namedtuple
from multiprocessing import Pool
import time
import operator
import cPickle
import os
import argparse
import logging
import re

from prepare_taxdump_refseq import find_files, Annotation
import taxtree


# This needs to be defined at the top level for multiprocessing to work.
Hit = namedtuple("Hit", ["target_accno", "identity", "matches", "mismatches", 
    "score", "startpos", "endpos", "fragment_length", "fragment_coverage"])


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
            default="/shared/db/NCBI/taxonomy/taxdump/",
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
    parser.add_argument("--identity_difference", type=float, metavar="D",
            default=5.0,
            help="Maximum identity difference between highest and lowest hit for each fragment. Floating point between 0.0-100.0 [%(default)s].")
    parser.add_argument("--remove_nondiscriminative", dest="remove_nondiscriminative", action="store_false",
            default=True, # TODO: This is unintuitive: reverse wording
            help="Remove fragments that match to more than one entity at the specific taxonomic level [%(default)s].")
    parser.add_argument("--maxprint", dest="maxprint", metavar="N", type=int,
            default=50,
            help="Maximum number of hit annotations to print [%(default)s].")
    parser.add_argument("--print_all_hit_annotations", dest="print_all_hit_annotations", action="store_true",
            default=False,
            help="Print a listing of all hit annotations (and not only hits from filtered peptide fragments) [%(default)s].")
    parser.add_argument("--interactive", dest="interactive", action="store_true",
            default=False,
            help="Load all heavy stuff and then run interactively [%(default)s].")
    parser.add_argument("--output", dest="output",
            default="",
            help="Write results to this filename [results/FILE.results].")


    devoptions = parser.add_argument_group("Developer options", "Voids warranty ;)")
    devoptions.add_argument("--loglevel", choices=["INFO", "DEBUG", "VERBOSE"],
            default="INFO", 
            help="Set logging level [%(default)s].")
    devoptions.add_argument("--logfile", dest="logfile", 
            default="proteotyping.log",
            help="Filename for log output [%(default)s].")
    devoptions.add_argument("--walk", dest="walk", action="store_true",
            default=False,
            help="Instead of only visiting 'leaf' nodes when printing, walk the entire distance from the root node down [%(default)s].")
    devoptions.add_argument("--numCPUs", dest="numCPUs", type=int,
            default=16,
            help="Number of CPUs to utilize in parallel regions [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit()

    options = parser.parse_args(argv[1:])
    logger = logging.getLogger()
    if options.loglevel == "VERBOSE":
        logger.setLevel(0)
    else:
        logger.setLevel(options.loglevel)
    fh = logging.FileHandler(options.logfile)
    ch = logging.StreamHandler()
    file_formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    console_formatter = logging.Formatter("%(levelname)s: %(message)s")
    fh.setFormatter(file_formatter)
    ch.setFormatter(console_formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)

    logging.info("----------========= LOGGING STARTED {} ========---------".format(time.strftime("%c")))

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


def parse_blat_output(filename, options):
    """Parses blat output. Filters out hits based on user critera. 
    
    Critera used for filtering is:
     options.identity
     options.fragment_length
     options.fragment_coverage
    
    """

    fragment_ids = set()
    hit_counter = 0
    hits = {}

    logging.info("Parsing and performing absolute filtering of hits from {}".format(filename))
    with open(filename) as blast8:
        for line in blast8:
            blast8_line = line.split()
            fragment_id = blast8_line[0]

            # Fragment and hit counting logic
            if fragment_id not in fragment_ids:
                fragment_ids.add(fragment_id)
            hit_counter += 1

            target_accno = parse_accno(blast8_line[1])
            mapping = "{}::{}".format(fragment_id, target_accno)

            identity = float(blast8_line[2])
            if identity < options.identity:
                logging.log(0, "  {} passed identity.".format(mapping))
                continue
            matches = int(blast8_line[3])
            mismatches = int(blast8_line[4])
            startpos = int(blast8_line[8])
            endpos = int(blast8_line[9])
            fragment_length = int(fragment_id.split()[0].split("_")[-1])
            if fragment_length < options.fragment_length:
                logging.log(0, "  {} passed fragment length.".format(mapping))
                continue
            fragment_coverage = matches/fragment_length
            if fragment_coverage < options.fragment_coverage:
                logging.log(0, "  {} passed fragment coverage.".format(mapping))
                continue

            hit = Hit(target_accno, identity, matches, mismatches, 
                      0, startpos, endpos, fragment_length, fragment_coverage)
            try:
                hits[fragment_id].append(hit)
            except KeyError:
                hits[fragment_id] = [hit]

    num_hits_remaining = sum(map(len, [hitlist for hitlist in hits.itervalues()]))
    logging.info("Parsed {} hits for {} fragments.".format(hit_counter, len(fragment_ids)))
    logging.info("{} fragments with {} hits remain after absolute filtering ({} hits were removed).".format(len(hits), num_hits_remaining, hit_counter-num_hits_remaining))
    filter_relative(hits, options.identity_difference)
    return hits


def filter_relative(hits, identity_difference):
    """Perform relative filtering for each fragment.

    Removes hits from each fragment's hitlist that have identity less
    than max(hitlist)-identity_difference.
    Modifies hits in place!
    """

    logging.info("Performing relative filtering on hits for each fragment...")
    num_removals = 0
    for fragment, hitlist in hits.iteritems():
        highest_identity = max(hitlist, key=lambda hit: hit.identity).identity
        startlen = len(hitlist)
        hitlist = [hit for hit in hitlist if hit.identity >= highest_identity-identity_difference]
        postlen = len(hitlist)
        num_removals += startlen - postlen
        logging.log(0, "Relative filtering of hitlist for fragment {} removed {} out of {} hits. {} hits remain".format(fragment, startlen-postlen, startlen, postlen))
        hits[fragment] = hitlist
    logging.info("{} hits were removed by relative filtering criteria.".format(num_removals))


def informative_fragment(hitlist, tree, taxonomic_rank):
    """Return True if hitlist for a given fragment is considered informative.

    Searches the tree for the node corresponding to the first hit in the
    hitlist, then finds the ancestor at the correct taxonomic rank and creates
    a set of all nodes beneath that node to which all the remaining hits in the
    list are compared. If a hit is not to a node in the set, the fragment is
    not informative.

    Returns True if fragment is informative (if hitlist only contains hits that
    match beneath the same node of rank 'taxonomic_rank').
    Returns False otherwise.
    """
    
    if len(hitlist) > 0:
        hit_node = taxtree.search_for_accno(tree, hitlist[0].target_accno)
    else:
        logging,warning("Empty hitlist! (This should never happen)".format(hitlist))
        return False

    if hit_node:
        # TODO: Potential bug here, I assume only one hit_node in the list.
        logging.debug("Found hit to node {}.".format(hit_node[0].name))
        if hit_node[0].rank == taxonomic_rank:
            ancestor = hit_node[0]
        else:
            ancestors = hit_node[0].get_ancestors()
            try:
                ancestor = [n for n in ancestors if n.rank == taxonomic_rank][0]
            except:
                logging.warning("Found no ancestor of correct rank")
                return False
        logging.debug("Found ancestor {} at rank '{}'.".format(ancestor.name, ancestor.rank))
        child_nodes = []
        for descendant in ancestor.get_descendants():
            for accno in descendant.accno:
                child_nodes.append(accno)
        child_nodes = set(child_nodes)
    else:
        logging.warning("Found no node in the tree for {}.".format(hitlist[0].target_accno))
        return False

    for hit in hitlist:
        if hit.target_accno not in child_nodes:
            logging.debug("{} was not found in {}".format(hit.target_accno, child_nodes))
            return False
    logging.debug("All hits hit under the same taxonomic level, informative fragment!")
    logging.debug("Informative fragment hits to {}".format(hit_node[0].taxname))
    logging.debug("  Informative hits to {}".format(" ".join([hit.target_accno for hit in hitlist])))
    return True


def filter_parallel(fragment_hitlist):
    """Find discriminative hits.
    
    Designed to be used in multiprocessing: pool.map(filter_parallel, items).

    Uses global variables:
     options.taxonomic_rank
     tree

    Returns None if fragment is not informative.
    """

    fragment_id, hitlist = fragment_hitlist

    if options.remove_nondiscriminative:
        if not informative_fragment(hitlist, tree, options.taxonomic_rank):
            logging.debug("Fragment {} is not informative".format(fragment_id))
            return 
        if logging.getLogger().getEffectiveLevel() < 20:
            logging.debug("Fragment {} is informative with hits to:".format(fragment_id))
            for hit in hitlist:
                logging.debug("  {}".format(hit.target_accno))
    return (fragment_id, hitlist)


def filter_hits(hits, tree, options): 
    """Determine if hits are disciminative at options.taxonomic_rank.

    Runs in parallel using multiprocessing with options.numCPUs
    """

    logging.debug("Starting pool of {} workers to filter out discriminative fragments.".format(options.numCPUs))
    pool = Pool(options.numCPUs)
    result_list = pool.map(filter_parallel, hits.items())

    # Merge all results
    filtered_hits = {}
    for result in result_list:
        if result is not None:
            fragment_id, filtered_hitlist = result
            filtered_hits[fragment_id] = filtered_hitlist

    logging.info("Removed {} non-discriminative fragments. {} discriminative fragments with {} hits remain.".format(len(hits)-len(filtered_hits), 
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
        for idx, parent in enumerate(lineage[1:]):
            parent.count += leaf.count


def reset_tree(tree):
    """Resets counts on all nodes for use in interactive mode.
    """
    for node in tree.traverse():
        node.count = 0


def write_top_n_hits(outfilehandle, tree, taxonomic_rank, totalhits, n=10, walk=False):
    """Writes the top 'n' nodes with the highest counts (and thus percentages).

    outfilehandle is an already opened filehandle.
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
        outfilehandle.write("     %  #         TAXID       TAXNAME                         ACCNO(s)\n")
        for i in xrange(0,N):
            n = nodes_to_print[i]
            outfilehandle.write("{:>6.2f}  {:<8}  {:<10}  {:<30}  {:<}\n".format(100.0*n.count/float(totalhits), n.count, n.name, n.taxname, ";".join(n.accno)))
    else:
        outfilehandle.write("Nothing filtered through at {} level.".format(taxonomic_rank))


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
    logging.debug("Loading gene info from {}".format(gene_info_file))
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
    logging.debug("Gene info for {} genes loaded.".format(len(gene_info)))
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

    if missing_annotation:
        logging.warning("Couldn't find annotation for:\n{}.".format("\n".join(missing_annotation)))
    return gene_counts


def write_gene_counts(outfilehandle, gene_counts, gene_info, maxprint=50, sort=True):
    """Prints counts of annotated regions with some gene info.

    outfilehandle is an already opened filehandle.
    """
    if sort:
        counts = sorted(gene_counts.items(), key=operator.itemgetter(1), reverse=sort)
    printcounter = 0
    
    outfilehandle.write("geneID    #     Symbol             Description\n")
    for geneID, count in counts:
        if printcounter < maxprint:
            try:
                symbol, desc = gene_info[geneID][0:2]
            except KeyError:
                logging.warning("Couldn't find information for geneID {}.".format(geneID))
            outfilehandle.write("{:<9} {:<5} {:<18} {:<}\n".format(geneID, count, symbol, desc))
            printcounter += 1
        else:
            break
    outfilehandle.write("Printed {} out of {} hit annotated regions.\n".format(printcounter, len(counts)))


def load_annotation(annotation_pickle):
    """Loads annotation from previously prepared gff info pickle.
    """

    logging.debug("Loading annotations from {}".format(annotation_pickle))
    with open(annotation_pickle, 'rb') as pkl:
        annotation = cPickle.load(pkl)
        if len(annotation) == 0:
            errorstring = "{} contains no annotations!".format(annotation_pickle)
            logging.error(errorstring)
            raise IOError(errorstring)
    logging.debug("Annotations for {} genes loaded.".format(len(annotation)))
    return annotation


def wait_for_user_input(options):
    """Read input from terminal.
    """

    instructions = """
Ran with options: l: {}, c: {:2.2f}, i: {}, r: {}
Enter options to change filtering criteria and run again:
 fragment (l)ength, (i)dentity, fragment (c)overage, (r)emove noninformative.
 Some options can be combined with numbers, e.g. 'l 10' changes the filtering length to 10.
 Separate multiple options with ','.
 Type 'q' and press enter to quit.
Input: """.format(options.fragment_length, options.fragment_coverage, options.identity, options.remove_nondiscriminative)
    try:
        user_input = raw_input(instructions)
        for param in [a.strip() for a in user_input.split(",")]:
            if param.lower().startswith("l"):
                try:
                    options.fragment_length = int(param.split()[1])
                except ValueError:
                    print "ERROR: Cannot parse fragment length Try again."
                    return wait_for_user_input(options)
            elif param.lower().startswith("i"):
                try:
                    options.identity = float(param.split()[1])
                except ValueError:
                    print "ERROR: Cannot parse identity. Try again."
                    return wait_for_user_input(options)
            elif param.lower().startswith("c"):
                try:
                    options.fragment_coverage = float(param.split()[1])
                except ValueError:
                    print "ERROR: Cannot parse fragment coverage. Try again."
                    return wait_for_user_input(options)
            elif param.lower().startswith("r"):
                options.remove_nondiscriminative = not options.remove_nondiscriminative
            elif param.lower().startswith("q"):
                print "User exited."
                exit()
            else:
                print "ERROR: Option not recognized, try again."
                return wait_for_user_input(options)
        return options
    except KeyboardInterrupt:
        print "\nUser exited through Ctrl+C."
        exit()
    except EOFError:
        print "\nUser exited through Ctrl+D."
        exit()


def write_results(filename, tree, hits, filtered_hits, totalhits, annotation, options):
    """Write results to file.
    """

    if options.output:
        outfilename = options.output
    else:
        sourcename = os.path.basename(filename)
        if not (os.path.exists("results") and os.path.isdir("results")):
            os.mkdir("results")
        outfilename = "results/{}.results".format(os.path.splitext(os.path.basename(filename))[0])

    logging.info("Writing results to {}".format(outfilename))
    with open(outfilename, "w") as outfile:
        outfile.write("Results at rank '{}' for file {}\n".format(options.taxonomic_rank, filename))
        outfile.write("-"*70+"\n")
        write_top_n_hits(outfile, tree, options.taxonomic_rank, totalhits, n=options.display, walk=options.walk)
        outfile.write(" Total: {:<}\n".format(totalhits))

        outfile.write("-"*70+"\n")
        if options.print_all_hit_annotations:
            gene_counts = count_annotation_hits(hits, annotation)
            outfile.write("All hit annotated regions:\n")
            write_gene_counts(outfile, gene_counts, gene_info, options.maxprint)
        else:
            gene_counts = count_annotation_hits(filtered_hits, annotation)
            outfile.write("Annotated regions hit by filtered fragments:\n")
            write_gene_counts(outfile, gene_counts, gene_info, options.maxprint)


def main(filename, tree, options):
    """Main function that runs the complete pipeline logic.
    """

    hits = parse_blat_output(filename, options)
    discriminative_hits = filter_hits(hits, tree, options)
    totalhits = sum(map(len, [hits for hits in discriminative_hits.itervalues()]))

    if logging.getLogger().getEffectiveLevel() < 20:
        logging.log(0,"All filtered hits:")
        for fragment, hitlist in discriminative_hits.iteritems():
            for hit in hitlist:
                logging.log(0,"  {} {}".format(fragment, hit))

    insert_hits_into_tree(tree, discriminative_hits)
    sum_up_the_tree(tree)

    write_results(filename, tree, hits, discriminative_hits, totalhits, annotation, options)


if __name__ == "__main__":

    options = parse_commandline(argv)

    tree = load_taxtree(options.taxtree_pickle, options.taxdumpdir, 
            options.id_gi_accno_pickle, options.rebase_tree)
    gene_info = load_gene_info(options.gene_info_file)
    annotation = load_annotation(options.accno_annotation_pickle)

    if options.interactive:
        while True:
            reset_tree(tree)
            for filename in options.FILE:
                main(filename, tree, options)
            options = wait_for_user_input(options)
    else:
        for filename in options.FILE:
            main(filename, tree, options)
            reset_tree(tree)
