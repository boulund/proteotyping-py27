#!/usr/bin/env python
# Fredrik Boulund 2014
# Convenience script to run blat on multiple files

from sys import argv, exit
from subprocess import Popen, PIPE
from multiprocessing import Pool
import shlex
import argparse
import logging
from os import path, makedirs


def parse_commandline():
    """Parses commandline.
    """

    desc="""Convenience wrapper for running multiple BLAT instances in parallel. Fredrik Boulund 2014."""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("DB", metavar="DBFILE",
            help="Database to search against.")
    parser.add_argument("OUTDIR", metavar="DIR",
            help="Directory where results are written.")
    parser.add_argument("FASTA", metavar="FILE", nargs="+",
            help="""FASTA file to search database with. Many files can be
            specified on the command line, make sure to put them last in the
            command line!""")
    parser.add_argument("--blat_options", metavar="'S'", dest="blat_options",
            default="-out=blast8 -t=dnax -q=prot -tileSize=5 -stepSize=5 -minScore=10 -minIdentity=85",
            help="String of options to send to blat ['%(default)s'].")
    parser.add_argument("-N", dest="N", metavar="N", type=int,
            default=10,
            help="Number of processes to run at the same time [%(default)s].")
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


def parallel_run_blat(filename):
    """Runs blat on a single file.
    """
    database = options.DB
    blat_options = options.blat_options
    outdir = options.OUTDIR

    if not path.exists(outdir):
        makedirs(outdir)

    outputfile = path.join(outdir, path.splitext(path.basename(filename))[0]+".blast8")

    blat_call = shlex.split("blat {} {} {} {}".format(database, filename, blat_options, outputfile))
    blat = Popen(blat_call, stdout=PIPE, stderr=PIPE)
    blat_output = blat.communicate()
    if blat.returncode != 0:
        logging.error("BLAT: stdout: {}\n  stderr:{}".format(blat_output[0], blat_output[1]))


if __name__ == "__main__":
    options = parse_commandline()

    if options.N < len(options.FASTA):
        num_workers = options.N
    else:
        num_workers = len(options.FASTA)

    logging.debug("Creating pool of {} workers.".format(num_workers))
    pool = Pool(processes=num_workers)
    result = pool.map(parallel_run_blat, options.FASTA)
