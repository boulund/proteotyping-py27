#!/usr/bin/env bash
# Fredrik Boulund 2014
# Download NCBI taxonomy dump to dir named after today's date and gunzip files

date=`date +%Y%m%d`
read -p "About to download NCBI taxonomy dump from ftp.ncbi.nih.gov to taxdump_$date"$'\nCtrl-C to abort, Press a key to continue.'
mkdir taxdump_$date
cd taxdump_$date
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt > >(tee -a download.log) 2> >(tee -a download.log >&2)
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz > >(tee -a download.log) 2> >(tee -a download.log >&2)
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid.readme > >(tee -a download.log) 2> >(tee -a download.log >&2)
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz > >(tee -a download.log) 2> >(tee -a download.log >&2)
tar -xf taxdump.tar.gz > >(tee -a download.log) 2> >(tee -a download.log >&2)
gunzip gi_taxid_nucl.dmp.gz > >(tee -a download.log) 2> >(tee -a download.log >&2)
