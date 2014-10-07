#!/usr/bin/env bash
# Fredrik Boulund 2014
# Download NCBI taxonomy dump to current dir and prepare files

mkdir taxdump_test
cd taxdump_test
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt | tee -a download_taxdump.log
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz | tee -a download_taxdump.log
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid.readme | tee -a download_taxump.log
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz | tee -a download_taxdump.log
gunzip taxdump.tar.gz | tee -a download_taxdump.log
gunzip gi_taxid_nucl.dmp | tee -a download_taxdump.log
