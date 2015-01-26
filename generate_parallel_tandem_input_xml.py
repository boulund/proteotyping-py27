#!/usr/bin/env python
# Fredrik Boulund 2015

from sys import argv, exit
import argparse
import logging
from os import path, makedirs


XMLSTRING = """<?xml version="1.0"?>
<bioml>
	<note type="input" label="list path, default parameters">default_input.xml</note>
	<note type="input" label="list path, taxonomy information">taxonomy.xml</note>
	<note type="input" label="protein, taxon">bacteria</note>
	<note type="input" label="spectrum, path">{input_mzxml}</note>
	<note type="input" label="output, path">{output}.xml</note>
</bioml>
"""


def parse_commandline():
    """Parses commandline.
    """

    desc="""Generate input XML files for Parallel Tandem (X!!tandem). Fredrik Boulund 2015"""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-o", metavar="DIR",
            default=".",
            help="Directory where results are written. Default current dir.")
    parser.add_argument("mzXML", metavar="FILE", nargs="+",
            help="""mzXML file to search database with. Many files can be
            specified on the command line, but make sure to put them last in the
            command line!""")
    parser.add_argument("--loglevel", 
            choices=["INFO","DEBUG"],
            default="DEBUG",
            help="Set logging level [%(default)s]")
    
    if len(argv)<2:
        parser.print_help()
        exit()

    options = parser.parse_args()
    logging.basicConfig(level=options.loglevel)
    return options


def main(filename):
    """ Main function.
    """
    
    fn_abs = path.abspath(filename)
    fn_base = path.basename(filename)
    fn = path.splitext(fn_base)[0]

    with open("input_"+fn+".xml", 'w') as inputxml:
        inputxml.write(XMLSTRING.format(input_mzxml=fn_abs, output=fn))


if __name__ == "__main__":
    options = parse_commandline()

    for filename in options.mzXML:
        if not filename.lower().endswith(".mzxml"):
            logging.warning("{} doesn't have a mzxml file ending.".format(filename))
        main(filename)
