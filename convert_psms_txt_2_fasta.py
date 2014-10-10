#!/usr/bin/env python
# Fredrik Boulund 2014
# Convert psms.txt files from mascot to fasta

from sys import argv, exit
import argparse
import logging
from os import path, makedirs


def parse_commandline():
    """Parse commandline.
    """

    desc = """Convert psms.txt output files to FASTA. Fredrik Boulund 2014"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FILE", type=str, nargs="+",
            help="Filename of psms.txt file to convert to FASTA")
    parser.add_argument("-o", dest="outdir", metavar="DIR", type=str,
            default="fasta",
            help="Output directory [%(default)s].")
    parser.add_argument("--remove_duplicates", dest="remove_duplicates", action="store_true",
            default=False,
            help="Remove duplicates [%(default)s].")
    parser.add_argument("--loglevel", 
            choices=["INFO", "DEBUG"],
            default="DEBUG",
            help="Set logging level [%(default)s].")

    if len(argv)<2:
        parser.print_help()
        exit()

    options = parser.parse_args()

    logging.basicConfig(level=options.loglevel)
    return options


def convert_psms_to_fasta(psms, outdir, remove_duplicates=False):
    """Parses psms.txt file and converts to FASTA, writes to outdir.
    """

    if not path.exists(outdir):
        makedirs(outdir)

    outfilename = path.join(outdir, path.splitext(path.basename(psms))[0]+".fasta")
    logging.debug("Writing FASTA to '{}'.".format(outfilename))
    with open(outfilename, 'w') as fastafile:
        with open(psms) as psmsfile:
            psmsfile.readline() # Skip header information
            for num, line in enumerate(psmsfile):
                hitinfo = line.split()
                #hitinfo = [info.strip('"') for info in hitinfo]
                sequence = hitinfo[3].strip('"')
                fastafile.write("{}\n{}\n".format(">{}".format(num), sequence))


def main(options):
    """Main.
    """
    for filename in options.FILE:
        convert_psms_to_fasta(filename, options.outdir, options.remove_duplicates)



if __name__ == "__main__":
    options = parse_commandline()
    main(options)
