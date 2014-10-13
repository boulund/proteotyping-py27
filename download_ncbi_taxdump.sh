#!/usr/bin/env bash
# Fredrik Boulund 2014
# Download NCBI taxonomy dump to current dir and gunzip files

date=`date +%Y%m%d`
read -p "About to download NCBI taxdump from ftp.ncbi.nih.gov to taxdump_$date"$'\nCtrl-C to abort, Press a key to continue.'
mkdir taxdump_$date
cd taxdump_$date
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt | tee -a download_taxdump.log
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz | tee -a download_taxdump.log
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid.readme | tee -a download_taxdump.log
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz | tee -a download_taxdump.log
tar -xf taxdump.tar.gz | tee -a download_taxdump.log
gunzip gi_taxid_nucl.dmp.gz | tee -a download_taxdump.log
