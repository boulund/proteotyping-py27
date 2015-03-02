#!/usr/bin/env python2.7
# Fredrik Boulund 2014
# Create a summary X!tandem output report 

from sys import argv, exit
import argparse
import logging
from os import path, makedirs
from lxml import etree


def parse_commandline():
    """Parse commandline.
    """

    desc = """Convert X!tandem XML output files to summary report. Fredrik Boulund 2015"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FILE", type=str, nargs="+",
            help="Filename of output XML file(s) to summarize")
    parser.add_argument("-o", dest="outdir", metavar="DIR", type=str,
            default="tandem_reports",
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



def generate_seqences_from_bioml_xml(xmlfile):
    """Generates sequence entries from X!tandem BIOML XML file.
    """

    for _, element in etree.iterparse(xmlfile):
        if element.tag == "group":
            if element.attrib["type"] == "model":
                group = element.attrib
            for domain in element.iterdescendants("domain"):
                yield group, domain.attrib



def create_tandem_output_report(xmlfile, outdir):
    """Converts X!tandem output BIOML XML to FASTA, writes to file in outdir.
    """

    if not path.exists(outdir):
        makedirs(outdir)

    outfilename = path.join(outdir, path.splitext(path.basename(xmlfile))[0]+"_report.txt")
    logging.debug("Writing report to '{}'.".format(outfilename))
    prev_group = ""
    prev_domain_string = ""
    totals = [0, 0]
    dom_counter = 0
    group_counter = 0
    first = True
    with open(outfilename, 'w') as reportfile:
        for group, domain in generate_seqences_from_bioml_xml(xmlfile):
            if group != prev_group:
                prev_group = group
                if not first:
                    reportfile.write("  x{}\n".format(dom_counter))
                    group_counter += dom_counter
                    reportfile.write("  Total: {}\n".format(group_counter))
                    group_counter = 0
                first = True
                dom_counter = 0
                reportfile.write("{id:<7} {mh:<12} {z} {sumI:<4} {maxI:<12} {fI:<10} {label}\n".format(**group))
                totals[0] += 1

            domain_string = "   {expect:<7} {hyperscore:<4} {seq}".format(**domain)
            dom_counter += 1
            totals[1] += 1
            if domain_string != prev_domain_string:
                prev_domain_string = domain_string
                if first:
                    first = False
                else:
                    reportfile.write("  x{}\n".format(dom_counter))
                    group_counter += dom_counter
                reportfile.write(domain_string)
        else:
            reportfile.write("  x{}\n".format(dom_counter))
            group_counter += dom_counter
            reportfile.write("  Total: {}\n".format(group_counter))

        reportfile.write("\nTotal number of spectra: {}\n".format(totals[0]))
    logging.info("Wrote {} peptides for {} spectra to {}.".format(totals[1], totals[0], outfilename))



def main(options):
    """Main.
    """
    for xmlfile in options.FILE:
        create_tandem_output_report(xmlfile, options.outdir)


if __name__ == "__main__":
    options = parse_commandline()
    main(options)
