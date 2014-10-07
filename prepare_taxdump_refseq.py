#!/usr/bin/env python
# Fredrik Boulund
# 2014-10-07
# Prepare the files from NCBI taxdump and NCBI refseq
# to reduce loading times for the proteotyping pipeline

from sys import argv, exit
import re
import logging
import argparse
import cPickle
import fnmatch
import os


def parse_fasta(filename):
    """Reads information from FASTA headers.
    
    Can read multiple FASTA headers from a single file.
    Information from headers is split into two parts on the first space.
    Returns a list of information tuples for each sequence in the file.
    Raises IOError if file is not in FASTA format.
    """
    def trim_header(h):
        return tuple(h[1:].strip().split(" ", 1))

    output = []
    with open(filename) as f:
        line = f.readline()
        if not line.startswith(">"):
            raise IOError("File does not appear to be valid FASTA\n  "+filename)
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


def filter_GIs_from_dmp_to_pickle(gi_taxid_dmp, gis, pickle_file):
    """Creates a pickled dictionary with gi:taxid mappings for a set of GIs.
    
    Raises TypeError if 'gis' is not a set.
    """
    if type(gis) is not set:
        raise TypeError("gis needs to be a set")
    if not os.path.isfile(gi_taxid_dmp):
        raise IOError("Taxdump file {} not found.".format(gi_taxid_dmp))
    with open(gi_taxid_dmp) as f:
        gi_taxid = {}
        for line in f:
            gi, taxid = line.split()
            gi = int(gi)
            if gi in gis:
                gi_taxid[int(taxid)] = gi
    with open(pickle_file, 'wb') as pkl:
        cPickle.dump(gi_taxid, pkl, -1) # Highest protocol available


def create_set_of_GIs(refdir, pattern):
    """Creates a set of GIs identified in a directory structure of reference
    sequences in FASTA format.
    """
    gis = set()
    for fasta_file in find_files(refdir, pattern):
        for seqinfo in parse_fasta(fasta_file):
            gi, accno = parse_gi_accno_from_headers(seqinfo[0])
            gis.add(gi)
    return gis


def parse_gi_accno_from_headers(header):
    """Parses GI and accession number from strings (in regular FASTA header format.
    """
    # Regex identifies:
    #   NCBI gi numbers (digits)
    #   NCBI GenBank accession numbers (e.g. NC_001122.22)
    # from FASTA headers such as gi|158333233|ref|NC_009925.1|
    # hit.group(1) is gi
    # hit.group(2) is accession
    gi_accno_regex = r'gi\|(\d+)\|.*ref\|(\w{1,2}_[\d\w]+\.\d{1,2})\|'
    hit = re.search(gi_accno_regex, header)
    if hit:
        return (int(hit.group(1)), hit.group(2))
    else:
        raise Exception("Couldn't parse GI, accno from {}.".format(header))


def parse_commandline(argv):
    """Parse commandline arguments."""

    desc = """Prepares a pickled dict of relevant gi:taxid mappings for quick-loading."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("TAXDUMP", help="Path to gi_taxid dump (two column format).")
    parser.add_argument("REFDIR", help="Path to directory with reference sequences in FASTA format. Walks subfolders.")

    parser.add_argument("-o", dest="output", help="Output file [default %(default)s]", default="gi_taxid.pkl")
    parser.add_argument("-g", dest="globpattern", help="glob pattern for identifying FASTA files [%(default)s]", default="*.fna")

    devoptions = parser.add_argument_group("Developer options")
    devoptions.add_argument("--loglevel", choices=["INFO", "DEBUG"], default="DEBUG", help="Set logging level.")

    if len(argv) < 2:
        parser.print_help()
        exit()

    options = parser.parse_args()
    logging.basicConfig(level=options.loglevel)
    return options


if __name__ == "__main__":

    options = parse_commandline(argv)

    logging.debug("Creating a set of GIs to include in the pickle...")
    gis = create_set_of_GIs(options.REFDIR, options.globpattern)
    logging.debug("Found {} GIs.".format(len(gis)))

    logging.debug("Filtering GIs from taxdmp...")
    filter_GIs_from_dmp_to_pickle(options.TAXDUMP, gis, options.output)
    logging.debug("Wrote {}.".format(options.output))
