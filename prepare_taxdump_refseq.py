#!/usr/bin/env python
# Fredrik Boulund
# 2014-10-07
# Prepare the files from NCBI taxdump and NCBI refseq
# to reduce loading times for the proteotyping pipeline

from sys import argv, exit
from collections import namedtuple
import time
import re
import logging
import argparse
import cPickle
import fnmatch
import os


# This needs to be global for pickle to work
Annotation = namedtuple("Annotation", ["geneID", "startpos", "endpos", "strand"])

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


def parse_gff(filename):
    """Generator that parses gff files and yields annotations.
    """
    Annotation_info = namedtuple("Annotation", ["accno", "geneID", "startpos", "endpos", "strand"])
    geneID_regex = re.compile(r'GeneID:(\d+)')
    with open(filename) as gff_file:
        if not gff_file.readline().startswith("##gff-version 3"):
            logging.warning("'{}' doesn't appear to be a valid gff file, skipping.".format(filename))
            return
        [gff_file.readline() for x in xrange(0,4)] # skip header lines
        for line in gff_file:
            info = line.split("\t")
            if not info[0].startswith("###") and info[2] == "gene":
                accno = info[0]
                startpos = int(info[3])
                endpos = int(info[4])
                strand = info[6]
                infokeys = info[8]
                hit = re.search(geneID_regex, infokeys)
                if hit:
                    geneID = hit.group(1)
                    yield Annotation_info(accno, geneID, startpos, endpos, strand)


def create_dict_accno_gff(refdir, pattern="*.gff"):
    """Creates a dict with accno:gene_location mappings.
    """
    accno_gff = {}
    for gff_file in find_files(refdir, pattern):
        for info in parse_gff(gff_file):
            try:
                accno_gff[info.accno].append(Annotation(info.geneID, info.startpos, info.endpos, info.strand))
            except KeyError:
                accno_gff[info.accno] = [Annotation(info.geneID, info.startpos, info.endpos, info.strand)]
    return accno_gff


def filter_GIs_from_dmp_to_pickle(gi_taxid_dmp, gi_accno):
    """Creates a dictionary with gi:taxid mappings for a set of GIs.
    """
    if not os.path.isfile(gi_taxid_dmp):
        raise IOError("Taxdump file {} not found.".format(gi_taxid_dmp))
    with open(gi_taxid_dmp) as f:
        gis = set(gi_accno.keys())
        taxid_gi_accno = {}
        inserted_counter = 0
        for line in f:
            gi, taxid = line.split()
            gi = int(gi)
            if gi in gis:
                try:
                    taxid_gi_accno[taxid].append((gi, gi_accno[gi]))
                    logging.debug("Taxid {} is now associated with {} GIs: ({}).".format(taxid, 
                        len(taxid_gi_accno[taxid]), taxid_gi_accno[taxid]))
                except KeyError:
                    taxid_gi_accno[taxid] = [(gi, gi_accno[gi])]
                inserted_counter += 1
                gis.remove(gi)
    logging.debug("Inserted {} GIs".format(inserted_counter))
    if logging.getLogger().getEffectiveLevel < 20:
        logging.debug("The following GIs were not included:")
        for gi in gis:
            logging.debug("{:>10}".format(gi))
    return taxid_gi_accno





def create_dict_gi_accno(refdir, pattern):
    """Creates a set of GIs identified in a directory structure of reference
    sequences in FASTA format.
    """
    gi_accno = {}
    for fasta_file in find_files(refdir, pattern):
        for seqinfo in parse_fasta(fasta_file):
            gi, accno = parse_gi_accno_from_headers(seqinfo[0])
            gi_accno[int(gi)] = accno
    return gi_accno


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

    desc = """Prepares a pickled dict of relevant gi:taxid mappings for quick-loading. Fredrik Boulund (c) 2014"""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("TAXIDS", help="Path to gi_taxid_nucl.dmp (two column format).")
    parser.add_argument("REFDIR", help="Path to directory with reference sequences in FASTA format and sequence annotation in GFF format. Walks subfolders.")

    parser.add_argument("--id_gi_accno_pickle", dest="id_gi_accno_pickle", metavar="FILE",
            default="id_gi_accno.pkl",
            help="Taxid_GI_ACCNO output file [default %(default)s]")
    parser.add_argument("--accno_annotation_pickle", dest="accno_annotation_pickle", metavar="FILE",
            default="accno_annotation.pkl",
            help="ACCNO_Annotation output file [default %(default)s]")
    parser.add_argument("-f", dest="globpattern_fasta", metavar="'GLOB'",
            default="*.fna",
            help="Glob pattern for identifying FASTA files [%(default)s]")
    parser.add_argument("-g", dest="globpattern_gff", metavar="'GLOB'",
            default="*.gff",
            help="Glob pattern for identifying GFF files [%(default)s].")

    devoptions = parser.add_argument_group("Developer options")
    devoptions.add_argument("--loglevel", choices=["INFO", "DEBUG"], 
            default="INFO", 
            help="Set logging level [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit()

    options = parser.parse_args()
    logging.basicConfig(level=options.loglevel)
    return options


if __name__ == "__main__":

    options = parse_commandline(argv)

    tic = time.time()
    logging.debug("Creating a set of GIs (based on FASTA headers) to include in the pickle...")
    gi_accno = create_dict_gi_accno(options.REFDIR, options.globpattern_fasta)
    logging.debug("Found {} GIs in {} seconds.".format(len(gi_accno.keys()), time.time()-tic))

    tic = time.time()
    logging.debug("Filtering GIs from {}...".format(options.TAXIDS))
    id_gi_accno = filter_GIs_from_dmp_to_pickle(options.TAXIDS, gi_accno)
    logging.debug("{} entries remaining".format(len(id_gi_accno.keys())))
    logging.debug("Time to filter GIs: {}.".format(time.time()-tic))

    tic = time.time()
    logging.debug("Writing pickle to '{}'.".format(options.id_gi_accno_pickle))
    with open(options.id_gi_accno_pickle, 'wb') as pkl:
        cPickle.dump(id_gi_accno, pkl, -1) # Highest protocol available
    logging.debug("Time to write pickle: {}.".format(time.time()-tic))

    tic = time.time()
    logging.debug("Creating annotation dict (based on GFF files) to include in the pickle...")
    accno_gff = create_dict_accno_gff(options.REFDIR, options.globpattern_gff)
    logging.debug("Found {} annotations in {} seconds.".format(len(accno_gff.keys()), time.time()-tic))

    tic = time.time()
    logging.debug("Writing pickle to '{}'.".format(options.accno_annotation_pickle))
    with open(options.accno_annotation_pickle, 'wb') as pkl:
        cPickle.dump(accno_gff, pkl, -1) # Highest protocol available
    logging.debug("Time to write pickle: {}.".format(time.time()-tic))


