Proteotyping pipeline
=====================
Performs bactereotyping from amino acid fragments produced by mass spectrometry. 


References and database required
--------------------------------
The pipeline depends on files from NCBI taxonomy (see listing below), as well
as a set of reference genomes (e.g. NCBI bacterial). Some pre-processing of
these genomes is required to reduce runtimes (just parsing the gi:taxid
mappings take like 5 minutes(!). This preprocessing is done by
"prepare_taxdump_refseq.py" and produces a pickled file, "id_gi_accno.pkl",
that is used by "proteotyping.py". 


Files required (download_ncbi_taxdump.sh gets the most recent ones from NCBI
FTP server):

 * names.dmp
 * nodes.dmp
 * gi_taxid_nucl.dmp






Python 2.7 Virtualenv
=====================
Prepare a virtual environment with Python 2.7 (I've been using Python 2.7.5).
Make sure you use a virtualenv version greater than 1.10.1 because it contains
a bug that makes it impossible to install PyQt4 properly. I used version
1.11.6.  Initialize the virtualenv with the '--always-copy' command to ensure
it doesn't symlink to the globally installed libraries. Further, the
environment requires the following modules:

   * lxml (3.3.5)
   * MySQL-python (1.2.5)
   * numpy (1.8.1)
   * PyQt4             <-- See instructions below! Not really required.
   * ete2 (2.2.1072)   <-- Requires lxml, MySQL-python, numpy, PyQt4

NOTE: I never actually got PyQt4 installed properly. You get an annoying error
message but unless you're going to print some pretty trees you can do without
it.
Make sure to install ete2 as the last package, as it depends on all the others.

    virtualenv --always-copy venv
    source venv/bin/activate
    pip install lxml MySQL-python numpy
    [Follow instructions below to install SIP and PyQt4]
    pip install ete2


PyQt4
-----
To install PyQt4 you need to first make sure Qt4 is installed on the machine,
on my machine (Terra) it required the following packages from the RedHat
Server-6 repo:

  * qt (1:4.6.2-28)
  * qt-devel (1:4.6.2-28)
  * qt-sqlite (1:4.6.2-28)
  * qt-x11 (1:4.6.2-28)

They install the qmake binary as /usr/bin/qmake-qt4. Then follow instructions
below to install SIP. After SIP is installed, PyQt4 can be compiled and
installed as per instructions below. It takes a while to compile! 
I'm thankful to user Flyte on StackOverflow who provided the key information
to get this working in his post [1].

Links:

  1 http://stackoverflow.com/questions/19856927/install-sip-pyqt-on-virtualenv
  2 http://pyqt.sourceforge.net/Docs/sip4/
  3 http://pyqt.sourceforge.net/Docs/PyQt4/

SIP
---
    # Activate virtualenv! 
    cd venv
    mkdir build
    cd build
    wget http://sourceforge.net/projects/pyqt/files/sip/sip-4.16.2/sip-4.16.2.tar.gz
    tar -xf sip-4.16.2.tar.gz
    cd sip-4.16.2
    python configure.py --incdir=../../include/python2.7
    make
    make install

PyQt4
-----
    # Activate virtualenv!
    cd venv/build
    wget http://sourceforge.net/projects/pyqt/files/PyQt4/PyQt-4.11.1/PyQt-x11-gpl-4.11.1.tar.gz
    tar -xf PyQt-x11-gpl-4.11.1.tar.gz
    cd PyQt-x11-gpl-4.11.1
    python configure-ng.py --qmake /usr/bin/qmake-qt4 --static
    make
    make install


 
