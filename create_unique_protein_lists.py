#!/usr/bin/env python2.7
# Fredrik Boulund 2014
# Create lists of unique proteins in X!tandem output report 

from sys import argv, exit
import argparse
import logging
from os import path, makedirs
from lxml import etree
import re
import cPickle


def parse_commandline():
    """Parse commandline.
    """

    desc = """Convert X!tandem XML output files to summary report. Fredrik Boulund 2015"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FILE", type=str, nargs="+",
        help="Filename of output report file(s) to summarize")
    parser.add_argument("-m", dest="mapping",
        help="Path to UniProt idmapping.dat.")
    parser.add_argument("-M", dest="mappings_pickle",
        help="Path to pickle with gi:uniprot mappings to use instead of reading idmapping.dat.")
    parser.add_argument("-o", dest="outdir", metavar="DIR", type=str,
        default="unique_protein_lists",
        help="Output directory [%(default)s].")
    parser.add_argument("-n", "--newline", dest="newline",
        choices=["windows", "unix"],
        default="windows", 
        help="Specify what line break type to use [%(default)s].")
    parser.add_argument("--loglevel", 
        choices=["INFO", "DEBUG"],
        default="INFO",
        help="Set logging level [%(default)s].")

    if len(argv)<2:
        parser.print_help()
        exit()

    options = parser.parse_args()

    logging.basicConfig(level=options.loglevel)
    return options


def read_uniprot_gi_mappings(idmappingfile):
    """ Read gi:uniprot mappings into dictionary.
    """

    gi_mappings = {}
    with open(idmappingfile) as f:
        for line in f:
            upid, kind, id = line.strip().split()
            if "GI" in kind:
                gi_mappings[id] = upid
    
    return gi_mappings


def find_unique_proteins(report):
    """ Find unique proteins in X!!tandem output report.
    """

    giregex = re.compile(r'gi\|(\d+)\|')
    
    unique_gis = set()

    with open(report) as f:
        for line in f:
            hit = re.search(giregex, line)
            if hit:
                unique_gis.add(hit.group(1))
    return unique_gis
            


def write_unique_protein_list(unique_gis, gi_mappings, newline, outdir, outputfile):
    """ Write list of unique proteins to file.
    """

    if not path.isdir(outdir):
        makedirs(outdir)

    outfile = path.join(outdir, outputfile)

    print "Writing to {}".format(outfile)
    if "windows" in newline:
        lb = "\r\n"
    elif "unix" in newline:
        lb = "\n"
    else:
        lb = "\n"

    nomaps = []
    with open(outfile, 'w') as f:
        f.write("UniProtID"+lb)
        for gi in unique_gis:
            try:
                f.write(gi_mappings[gi]+lb)
            except KeyError:
                logging.warning("Found no match for {}".format(gi))
                nomaps.append(gi)

        f.write("Found no mapping between gi and UniProt ID for:"+lb)
        for gi in nomaps:
            f.write(gi+lb)


    


def main(options):
    """Main.
    """

    if options.mappings_pickle and path.isfile(options.mappings_pickle):
        logging.debug("Found pickle, loading...")
        gi_mappings = cPickle.load(open(options.mappings_pickle, 'rb'))
    elif options.mappings_pickle and not path.isfile(options.mappings_pickle):
        logging.debug("Couldn't find pickle, creating...")
        gi_mappings = read_uniprot_gi_mappings(options.mapping)
        cPickle.dump(gi_mappings, open(options.mappings_pickle, "wb"))
    else:
        logging.warning("Loading mappings from file without creating pickle...")
        gi_mappings = read_uniprot_gi_mappings(options.mapping)

    for report in options.FILE:
        logging.debug("Finding unique proteins for {}".format(report))
        unique_gis = find_unique_proteins(report)
        logging.debug("Found {} unique proteins".format(len(unique_gis)))
        write_unique_protein_list(unique_gis, gi_mappings, options.newline, options.outdir, report+".unique_proteins.txt")
        
        


if __name__ == "__main__":
    options = parse_commandline()
    main(options)
