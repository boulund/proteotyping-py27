#!/usr/bin/env bash
# Fredrik Boulund 2014
# Download NCBI microbial FASTA (.fna) and annotations (.gff) to dir named
# after today's date and uncompress files

date=`date +%Y%m%d`
read -p "About to download bacterial genomes from ftp.ncbi.nih.gov to $date"$'\nCtrl-C to abort, Press a key to continue.'
mkdir $date
cd $date
wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz > >(tee -a download.log) 2> >(tee -a download.log >&2)
wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gff.tar.gz > >(tee -a download.log) 2> >(tee -a download.log >&2)
gunzip all.fna.tar.gz > >(tee -a download.log) 2> >(tee -a download.log >&2)
gunzip all.gff.tar.gz > >(tee -a download.log) 2> >(tee -a download.log >&2)
