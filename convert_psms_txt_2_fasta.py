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
    parser.add_argument("--confidence_level", dest="confidence_level", metavar="LVL",
            choices=["Low", "Middle", "High"],
            default="High",
            help="Minimum confidence level of peptides to include in the FASTA output [%(default)s].")
    parser.add_argument("--psm_ambiguity", dest="psm_ambiguity", metavar="AMB", 
            choices=["Any", "Selected", "Unambiguous"],
            default="Unambiguous",
            help="Minimum PSM PSM ambiguity level of peptides to include in the FASTA output [%(default)s].")
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


def convert_psms_to_fasta(psms, outdir, remove_duplicates=False, confidence_level="High", psm_ambiguity="Unambiguous"):
    """Parses psms.txt file and converts to FASTA, writes to outdir.
    """

    if confidence_level == "Low":
        confidence_level = set(["Low", "Middle", "High"])
    elif confidence_level == "Middle":
        confidence_level = set(["Middle", "High"])
    elif confidence_level == "High":
        confidence_level = set(["High"])
    else:
        logging.warning("Confidence level not properly specified: {}. Assuming 'High'".format(confidence_level))
        confidence_level = set(["High"])
    logging.debug("Confidence level: {}.".format(confidence_level))
    if psm_ambiguity == "Any":
        psm_ambiguity = set(["Unconsidered", "Selected", "Rejected", "Unambiguous"])
    elif psm_ambiguity == "Selected":
        psm_ambiguity = set(["Selected", "Unambiguous"])
    elif psm_ambiguity == "Unambiguous":
        psm_ambiguity = set(["Unambiguous"])
    else:
        logging.warning("PSM ambiguity not properly specified: {}. Assuming 'Unambiguous'".format(psm_ambiguity))
        psm_ambiguity = set(["Unambiguous"])
    logging.debug("PSM ambiguity: {}.".format(psm_ambiguity))

    if not path.exists(outdir):
        makedirs(outdir)

    outfilename = path.join(outdir, path.splitext(path.basename(psms))[0]+".fasta")
    logging.debug("Writing FASTA to '{}'.".format(outfilename))
    write_counter = 0
    with open(outfilename, 'w') as fastafile:
        with open(psms) as psmsfile:
            psmsfile.readline() # Skip header information
            for num, line in enumerate(psmsfile):
                hitinfo = line.split()
                hitinfo = [info.strip('"') for info in hitinfo]
                confidence = hitinfo[0]
                sequence = hitinfo[3]
                ambiguity = hitinfo[5]
                if confidence in confidence_level and ambiguity in psm_ambiguity:
                    logging.debug("Writing {} with {} and {}.".format(sequence, confidence, ambiguity))
                    fastafile.write("{}\n{}\n".format(">{}".format(num), sequence))
                    write_counter += 1
    logging.info("Wrote {} sequences to {}.".format(write_counter, outfilename))

def main(options):
    """Main.
    """
    for filename in options.FILE:
        convert_psms_to_fasta(filename, options.outdir, 
                options.remove_duplicates, options.confidence_level, options.psm_ambiguity)



if __name__ == "__main__":
    options = parse_commandline()
    main(options)
