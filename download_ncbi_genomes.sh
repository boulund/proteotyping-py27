#!/usr/bin/env bash
# Fredrik Boulund 2014
# Download NCBI microbial FASTA (.fna) and annotations (.gff) current dir and
# uncompress files

date=`date +%Y%m%d`
read -p "About to download bacterial genomes from ftp.ncbi.nih.gov to bacterial_$date"$'\nCtrl-C to abort, Press a key to continue.'
mkdir bacterial_$date
cd bacterial_$date
wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz | tee -a download.log
wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gff.tar.gz | tee -a download.log
tar -xf all.fna.tar.gz
tar -xf all.gff.tar.gz
