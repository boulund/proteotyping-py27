#!/usr/bin/env bash
# Fredrik Boulund 2014
# Download NCBI microbial FASTA (.fna) and annotations (.gff) current dir and
# uncompress files

read -p "About to download bacterial genomes from ftp.ncbi.nih.gov\nCtrl-C to abort, Press a key to continue."
mkdir NCBI_microbial
cd NCBI_microbial
wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz | tee -a download.log
wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gff.tar.gz | tee -a download.log
gunzip all.fna.tar.gz
gunzip all.gff.tar.gz
