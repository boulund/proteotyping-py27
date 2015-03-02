#!/usr/bin/env python2.7
# Fredrik Boulund 2015
# Fetch taxids from GenBank via ACCNO
# Input is list of lines with GI ACCNO separates by whitespace

from sys import argv, exit
from requests import get

if len(argv) < 2:
    print "usage: efetch_taxids.py FILE"
    print "Fredrik Boulund 2015"
    print "Fetch TAXIDs from Genbank using ACCNOs."
    print "Enter filename of file with one GI and ACCNO per line, separated by whitespace."
    exit()

with open(argv[1]) as f:
    for line in f:
        gi, accno = line.split()
        payload = {"db": "nuccore", 
                   "id": accno,
                   "rettype": "fasta",
                   "retmode": "xml"}
        xml = get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=payload)
        taxid = xml.text.split("taxid>")[1].split("<")[0]
        print taxid, gi, accno
