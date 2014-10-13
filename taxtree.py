#!/usr/bin/env python
# Fredrik Boulund 
# 2014
# Taxonomic trees

from sys import argv, exit
import re
import cPickle
import time
import logging
import os
from ete2 import Tree



def read_gi_accno_from_fasta_headers(headers):
    """Parse gi and accno from fasta headers."""
    gi_accno = {}
    with open(headers) as f:
        # Regex identifies:
        #   NCBI gi numbers (digits) 
        #   NCBI GenBank accession numbers (e.g. NC_001122.22) 
        #   Leftover information in the header, separated from the rest with whitespace
        # from FASTA headers.
        # hit.group(1) is gi
        # hit.group(2) is accession
        # hit.group(3) is other information
        gi_ref_regex = re.compile(r'gi\|(\d+)\|ref\|(\w{1,2}_\d+\.\d{1,2})\| (.+)')
        for line in f:
            hit = re.search(gi_ref_regex, line)
            if hit:
                gi_accno[hit.group(1)] = (hit.group(2), hit.group(3))
    return gi_accno


def load_ncbi_tree_from_taxdump(dumpdir, taxid_gi_accno_pickle, rebase=0):
    """Loads NCBI taxonomy tree from NCBI taxdump files names.dmp, nodes.dmp.

    The rebase parameter can be used to rebase the tree to whatever node, e.g.
    specifying rebase="2" will return only the subtree from node 2 (bacteria).

    Uses mapping of GI to taxid (node name) from NCBI file gi_taxid_nucl.dmp.
    Credit for this code goes to jhcepas, large parts of code taken from:
      https://github.com/jhcepas/ncbi_taxonomy/blob/master/update_taxadb.py
    """

    namesdump = dumpdir+"names.dmp"
    nodesdump = dumpdir+"nodes.dmp"
    if not os.path.isfile(namesdump):
        raise Exception("Taxdump file {} not found.".format(namesdump))
    if not os.path.isfile(nodesdump):
        raise Exception("Taxdump file {} not found.".format(nodesdump))

    # Different name types can be included in the nodes if desired,
    # just add them to the set below.
    name_types = set([]) 
        #"acronym", "anamorph", 
        #"authority", "blast name", "common name", 
        #"equivalent name", "genbank acronym", 
        #"genbank anamorph", "genbank common name", 
        #"genbank synonym", "includes", "in-part", 
        #"misnomer", "misspelling", "synonym", 
        #"teleomorph", "type material"])

    logging.debug("Loading taxid:gi:accno mappings...")
    with open(taxid_gi_accno_pickle, 'rb') as pickled:
        taxid2gi_accno = cPickle.load(pickled)
    logging.debug("{:>10} taxid:gi:accno mappings loaded.".format(len(taxid2gi_accno)))

    parent2child = {}
    name2node = {}
    node2taxname = {}
    taxid2gi = {}
    synonyms = {}
    name2rank = {}
    logging.debug("Loading node names...")
    with open(namesdump) as f:
        for line in f:
            fields =  line.split("\t|\t")
            nodename = fields[0]
            name_type = fields[3].strip("\t|\n").lower() 
            taxname = fields[1]
            try:
                if name_type == "scientific name":
                    node2taxname[nodename] = taxname
                    synonyms[nodename].add(taxname)
                elif name_type in name_types:
                    synonyms[nodename].add(taxname)
            except KeyError:
                synonyms[nodename] = set([taxname])
    logging.debug("{:>10} names loaded.".format(len(node2taxname)))
    logging.debug("{:>10} synonyms loaded.".format(sum([len(c[1]) for c in synonyms.iteritems()])))

    logging.debug("Loading nodes...")
    counter_inserted = 0
    with open(nodesdump) as f:
        removals = 0
        for i, line in enumerate(f):
            if i%100000 == 0:
                logging.debug("Processed {:>7} nodes".format(i))
            fields = line.split("\t|\t")
            nodename = fields[0]
            parentname = fields[1]
            n = Tree()
            n.name = nodename # An int stored as str. Seems to work best that way with ETE2.
            n.taxname = node2taxname[nodename]
            n.rank = fields[2]
            n.synonyms = synonyms[nodename]
            n.gi = []
            n.accno = []
            n.score = -1
            n.count = 0
            try:
                n.gi = [gi_accno_tuple[0] for gi_accno_tuple in taxid2gi_accno[n.name]]
                n.accno = [gi_accno_tuple[1] for gi_accno_tuple in taxid2gi_accno[n.name]]
                #logging.debug("Inserted genome gi and accno into node {}".format(n.name))
                counter_inserted += 1
            except KeyError:
                pass
            parent2child[nodename] = parentname
            name2node[nodename] = n
    logging.debug("{:>10} nodes loaded.".format(len(name2node)))
    logging.debug("{:>10} genomes inserted into nodes.".format(counter_inserted))
    if counter_inserted < len(taxid2gi_accno):
        logging.warning("{:>10} gi and accno's were not inserted into the tree!".format(len(taxid2gi_accno)-counter_inserted))

    logging.debug("Linking nodes...")
    for node in name2node:
       if node == "1":
           t = name2node[node]
       else:
           parent = parent2child[node]
           parent_node = name2node[parent]
           parent_node.add_child(name2node[node])
    logging.debug("Linking completed. Tree loaded.")
    
    if rebase:
        return t.search_nodes(name=rebase)[0].detach()
    else:
        return t


def search_for_accno(tree, accno):
    """Searches a tree for a node with node.accno attribute == accno.
    
    The accno argument can also be a list of accno's to search for, then a list
    of nodes that match the critera are returned.
    The result is returned in a list, even if only one node was found, for
    consistency with when multiple nodes are searched for.  
    """
    if isinstance(accno, list):
        result = []
        for acession in accno:
            result.append(search_for_accno(tree, acession)[0])        
        return result
    elif isinstance(accno, str):
        return [node for node in tree.get_leaves(lambda n: bool(n.accno)) if accno in node.accno]
    else:
        raise Exception("Can't search for {} of type {}.".format(accno, type(accno)))


if __name__ == "__main__":
    # Some prototype and debugging code...
    logging.basicConfig(level=logging.DEBUG)
    taxid_gi_accno_pickle = "id_gi_accno.pkl"
    dumpdir = "/shared/db/NCBI/taxonomy/taxdump/"

    tic = time.time()
    logging.debug("Reading tree from taxdump {}".format(dumpdir))
    tree = load_ncbi_tree_from_taxdump(dumpdir, taxid_gi_accno_pickle)
    logging.debug("Time to read tree {}".format(time.time()-tic))

    #tic = time.time()
    #logging.debug("Pickling tree...")
    #with open("taxtree.pkl", 'wb') as jar:
    #    cPickle.dump(tree, jar, -1)
    #logging.debug("Time to pickle tree {}".format(time.time()-tic))

    #tic = time.time()
    #logging.debug("Reading back pickled tree...")
    #with open("taxtree.pkl", 'rb') as jar:
    #    tree = cPickle.load(jar)
    #logging.debug("Time to read pickled tree {}".format(time.time()-tic))

    
    print "Searching for 1280..."
    aureus = tree.search_nodes(name="1280")
    agnetis = tree.search_nodes(name="985762")
    print aureus, agnetis

    print "Searching for common ancestor to {} {}...".format(aureus[0].name, agnetis[0].name)
    ancestor = tree.get_common_ancestor([aureus[0].name, agnetis[0].name])
    print ancestor.name

    print "Searching for NC_022226.1"
    S_aureus = search_for_accno(tree, "NC_022226.1")
    print S_aureus
    print S_aureus[0].name
    print S_aureus[0].accno

    print "Searching for multiple accnos: NC_022226.1, NC_022222.1"
    multinodes = search_for_accno(tree, ["NC_022226.1", "NC_022222.1"])
    print multinodes
    for node in multinodes:
        print node.name, node.accno

