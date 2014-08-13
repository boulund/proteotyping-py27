#!/usr/bin/env python
# Fredrik Boulund 
# 2014
# Taxonomic trees

from sys import argv, exit
from string import strip
import logging
import os
import cPickle
from ete2 import Tree

def load_ncbi_tree_from_dump(dumpdir):
    # Credit goes to jhcepas, code taken from:
    #   https://github.com/jhcepas/ncbi_taxonomy/blob/master/update_taxadb.py

    namesdump = dumpdir+"names.dmp"
    nodesdump = dumpdir+"nodes.dmp"
    if not os.path.isfile(namesdump):
        raise Exception("Taxdump file {} not found.".format(namesdump))

    parent2child = {}
    name2node = {}
    node2taxname = {}
    synonyms = set()
    name2rank = {}
    logging.debug("Loading node names...")
    for line in open(namesdump):
        fields =  map(strip, line.split("|"))
        nodename = fields[0]
        name_type = fields[3].lower()
        taxname = fields[1]
        if name_type == "scientific name":
            node2taxname[nodename] = taxname
        elif name_type in set(["synonym", "equivalent name", "genbank equivalent name",
                               "anamorph", "genbank synonym", "genbank anamorph", "teleomorph"]):
            synonyms.add( (nodename, taxname) )
    logging.debug(len(node2taxname), "names loaded.")
    logging.debug(len(synonyms), "synonyms loaded.")

    logging.debug("Loading nodes...")
    for line in open(nodesdump):
        fields =  line.split("|")
        nodename = fields[0].strip()
        parentname = fields[1].strip()
        n = Tree()
        n.name = nodename
        n.taxname = node2taxname[nodename]
        n.rank = fields[2].strip()
        parent2child[nodename] = parentname
        name2node[nodename] = n
    logging.debug(len(name2node), "nodes loaded.")

    logging.debug("Linking nodes...")
    for node in name2node:
       if node == "1":
           t = name2node[node]
       else:
           parent = parent2child[node]
           parent_node = name2node[parent]
           parent_node.add_child(name2node[node])
    logging.debug("Tree loaded.")
    
    return t, synonyms


if __name__ == "__main__":

    if not os.path.exists("taxtree.pkl"):
        tree, synonyms = load_ncbi_tree_from_dump(dumpdir="/shared/genomes/NCBI/taxonomy/taxdump/")
        cPickle.dump((tree, synonyms), open("taxtree.pkl", "wb"))
    else:
        tree, synonyms = cPickle.load(open("taxtree.pkl"))

    print "Searching for 1280..."
    aureus = tree.search_nodes(name="1280")
    agnetis = tree.search_nodes(name="985762")
    print aureus, agnetis
    print "Searching for common ancestor to {} {}...".format(aureus[0].name, agnetis[0].name)
    ancestor = tree.get_common_ancestor([aureus[0].name, agnetis[0].name])
    print ancestor.name






