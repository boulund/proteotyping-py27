#!/usr/bin/env python
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
    parser.add_argument("-o", dest="outdir", metavar="DIR", type=str,
            default="fasta",
            help="Output directory [%(default)s].")
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




def convert_tandem_bioml_to_fasta(xmlfile, outdir):
    """Converts X!tandem output BIOML XML to FASTA, writes to file in outdir.
    """

    if not path.exists(outdir):
        makedirs(outdir)

    outfilename = path.join(outdir, path.splitext(path.basename(xmlfile))[0]+".fasta")
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
    logging.info("Found {} unique source sequences in the XML.".format(len(sourceheaders)))
    logging.info("Wrote {} sequences to {}.".format(write_counter, outfilename))



def main(options):
    """Main.
    """
    for xmlfile in options.FILE:
        convert_tandem_bioml_to_fasta(xmlfile, options.outdir)


if __name__ == "__main__":
    options = parse_commandline()
    main(options)
