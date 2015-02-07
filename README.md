Proteotyping pipeline
=====================
Performs bactereotyping based on  amino acid fragments (peptides) produced by
mass spectrometry. It uses the placement of reference genomes in the taxonomic
tree to determine what peptide fragments are discriminative at a specified
taxonomic level and can thus be used to identify what species are present in
the sample. 

Overview of the workflow
------------------------

The files and programs in this suite of tools are the following:

    convert_psms_txt_2_fasta.py      (Legacy tool, deprecated)
    convert_tandem_xml_2_fasta.py
    create_tandem_search_report.py
    download_ncbi_genomes.sh
    download_ncbi_taxdump.sh
    prepare_taxdump_refseq.py
    proteotyping.py
    run_blat.py
    run_xtandem.py
    taxtree.py

To run the complete workflow, the following steps are required.

Preparations:

1. Run **download_ncbi_taxdump.sh** to download a recent copy of __NCBI Taxonomy__. 
   This is required to construct the taxonomic tree used in **proteotyping.py**.
2. Run **download_ncbi_genomes.sh** to download a recent copy of __NCBI bacterial__.
   This is required for placing reference genomes in the taxonomic tree, and 
   for finding annotations for discriminative hits in **proteotyping.py**.
2. Run **prepare_taxdump_refseq.py**. It requires the __taxdump__ from __NCBI Taxonomy__.
   You only need to run this once as the pickled files are reused.

The proteotyping workflow:

0. Run **ReAdW** to convert Thermo RAW files to mzXML or mzML. Could also use
   proteowizard tools. __NOTE: This must be done on Windows.__
1. Run **run_xtandem.py** to match mz spectra to peptides. We use __NCBI nr__
   and __HMP airways__ as references.
2. Run **create_tandem_search_report.py** to produce a simple overview of the 
   X!!tandem search results.
3. Run **convert_tandem_xml_2_fasta.py**.
4. Run **run_blat.py**. We use __NCBI bacterial__ as reference.
5. Run **proteotyping.py** to perform proteotyping based on the peptides in the
   BLAT output.


References and databases
------------------------

### Taxonomy
The taxonomic structure used in **proteotyping.py** is constructed from files
from __NCBI taxonomy__ and __NCBI bacterial__.  Some pre-processing of the
bacterial genomes is required. This preprocessing is done by
**prepare_taxdump_refseq.py** which produces a pickled file,
__id_gi_accno.pkl__, which is in turn used by **proteotyping.py**.

Files required by **prepare_taxdump_refseq.py**:

 * names.dmp
 * nodes.dmp
 * gi_taxid_nucl.dmp

The script **download_ncbi_taxdump.sh** gets the most recent ones from the NCBI
FTP server.

### Bacterial genomes
There are two main sources of reference sequences required in the workflow.
The first one is __NCBI bacterial__ which is used to produce the taxonomic reference
structure, as detailed above. This database is also used when mapping
the peptide fragments to reference genomes.

The second source of references is __NCBI nr__ and __HMP airways__ which are
used when identifying peptides from the mass spectra. This database should
contain as many variants of relevant peptides as possible, to increase the
chance of finding good peptide sequences. The next section summarizes what each
database is used for in the workflow.

### NCBI bacterial
__NCBI bacterial__ is used:
 * when constructing the taxonomic structure together with the taxdump from
   __NCBI taxonomy__. This step also requires the annotations for all genomes.
 * as a reference for BLAT when mapping peptides to reference genomes.

### NCBI nr and HMP Airways
__NCBI nr__  and __HMP airways__ are used:

 * as references for X!!tandem


Python 2.7 Virtualenv
=====================
To run the proteotyping workflow, prepare a virtual environment with Python 2.7
(I've been using Python 2.7.5 and virtualenv 1.11.6). Install the following
modules:

   * lxml (3.3.5)
   * MySQL-python (1.2.5)
   * numpy (1.8.1)
   * biopython (1.64)
   * bx-python (0.7.1)
   * master-sake (0.9.4.2)
   * pandas (0.15.2)
   * ete2 (2.2.1072)   -- Requires lxml, MySQL-python, numpy, and PyQt4 is optional.

NOTE: I never got PyQt4 installed properly. You get an annoying error message
everytime you run proteotyping.py, but unless you're going to print some pretty
trees you can do without it.  Make sure to install ete2 as the last package, as
it depends on some others and won't install them as dependencies.

    virtualenv venv
    source venv/bin/activate
    pip install lxml MySQL-python numpy biopython bx-python master-sake pandas 
    pip install ete2
