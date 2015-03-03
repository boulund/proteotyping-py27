#!/usr/bin/env python2.7
# Fredrik Boulund 2014
# Convert X!tandem output XML files to FASTA 

from sys import argv, exit
import argparse
import logging
from os import path, makedirs
from lxml import etree


def parse_commandline():
    """Parse commandline.
    """

    desc = """Convert X!tadenm XML output files to FASTA. Fredrik Boulund 2015"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FILE", type=str, nargs="+",
        help="Filename of output XML file to convert to FASTA")
    parser.add_argument("-d", dest="outdir", metavar="DIR", type=str,
        default="fasta",
        help="Output directory [%(default)s].")
    parser.add_argument("-O", dest="outfile", metavar="FILE", 
        default="",
        help="Output filename. If specified only one file is expected and outdir is disregarded.")
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



def generate_seqences_from_bioml_xml(xmlfile, first_seq_only=True):
    """Generates sequence entries from X!tandem BIOML XML file.

    first_seq_only  Default True. Returns only the first sequnece 
                    for each group (spectrum).
    """

    if first_seq_only:
        for _, element in etree.iterparse(xmlfile):
            if element.tag == "group":
                for child in element.iterdescendants("domain"):
                    yield element.attrib["label"], child.attrib["id"], child.attrib["expect"], child.attrib["seq"]
                    break
                continue
    else:
        # This generates ALL sequences from all groups (spectra). Not entirely useful.
        elements = (element for _, element in etree.iterparse(xmlfile))
        domains = (element for element in elements if element.tag == "domain")
        seqs = ((domain.attrib["id"], domain.attrib["expect"], domain.attrib["seq"]) for domain in domains)
        yield seqs




def convert_tandem_bioml_to_fasta(xmlfile, outdir, outfile):
    """Converts X!tandem output BIOML XML to FASTA, writes to file in outdir.
    """

    if outfile:
        outfilename = outfile
    else:
        outfilename = path.join(outdir, path.splitext(path.basename(xmlfile))[0]+".fasta")
        if not path.exists(outdir):
            makedirs(outdir)

    logging.debug("Writing FASTA to '{}'.".format(outfilename))
    sourceheaders = set()
    write_counter = 0
    with open(outfilename, 'w') as fastafile:
        for sourceheader, identity, expect, sequence in generate_seqences_from_bioml_xml(xmlfile):
            sourceheaders.add(sourceheader)
            logging.debug("Writing seq {} with length {}, expect {}.".format(identity, len(sequence), expect))
            header = ">{}:{}_{}".format(identity, expect, len(sequence))
            fastafile.write("{}\n{}\n".format(header, sequence))
            write_counter += 1
    logging.info("Wrote {} peptide fragments from {} unique protein sequences to {}.".format(write_counter, len(sourceheaders), outfilename))



def main(options):
    """Main.
    """
    for xmlfile in options.FILE:
        convert_tandem_bioml_to_fasta(xmlfile, options.outdir, options.outfile)


if __name__ == "__main__":
    options = parse_commandline()
    if options.outfile and len(options.FILE) > 1:
        logging.error("Cannot specify output filename with more than one file on command line")
        exit()
    main(options)
