#!/usr/bin/env python
# Fredrik Boulund 
# 2014
# Taxonomic trees

from sys import argv, exit
from string import strip
# import difflib #TODO
import time
import logging
import os
import cPickle
import copy
from ete2 import Tree


def read_accession_names_map(mapfile):
    """Reads accession: name mappings"""
    accno_name = {}
    with open(mapfile) as f:
        for line in f:
            accno, name = line.strip().split("\t")
            last_token = name.split()[-1].lower()
            if last_token.startswith(("contig","ctg")):
                name = " ".join(name.split()[:-1])
            elif last_token.startswith("genome"):
                name = " ".join(name.split()[:-2])
            accno_name[name] = accno
    return accno_name
            

def load_ncbi_tree_from_dump(dumpdir, accno_name_file):
    # Credit goes to jhcepas, code inspired by:
    #   https://github.com/jhcepas/ncbi_taxonomy/blob/master/update_taxadb.py

    namesdump = dumpdir+"names.dmp"
    nodesdump = dumpdir+"nodes.dmp"
    if not os.path.isfile(namesdump):
        raise Exception("Taxdump file {} not found.".format(namesdump))

    accno_name = read_accession_names_map(accno_name_file)

    name_types = set([]) 
        #"acronym", "anamorph", 
        #"authority", "blast name", "common name", 
        #"equivalent name", "genbank acronym", 
        #"genbank anamorph", "genbank common name", 
        #"genbank synonym", "includes", "in-part", 
        #"misnomer", "misspelling", "synonym", 
        #"teleomorph", "type material"])

    parent2child = {}
    name2node = {}
    node2taxname = {}
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
    with open(nodesdump) as f:
        removals = 0
        for i, line in enumerate(f):
            if i%100000 == 0:
                logging.debug("Processed {:>7} nodes".format(i))
            fields = line.split("\t|\t")
            nodename = fields[0]
            parentname = fields[1]
            n = Tree()
            n.name = nodename
            n.taxname = node2taxname[nodename]
            n.rank = fields[2]
            n.synonyms = synonyms[nodename]
            n.accno = []
            n.score = -1
            n.count = 0
            for synonym in n.synonyms:
                try:
                    accno = accno_name[synonym]
                    n.accno.append(accno)
                    accno_name.pop(synonym)
                    removals += 1
                    if len(n.accno)>1:
                        logging.warning("Added second accno {} to node {}.".format(accno, n.name))
                except KeyError:
                    pass
            parent2child[nodename] = parentname
            name2node[nodename] = n
    logging.debug("Removed {} items".format(removals, len(accno_name)))
    logging.debug("The following {} items still lack assignment:\n{}".format(len(accno_name.keys()),accno_name.keys()))
    logging.debug("{:>10} nodes loaded.".format(len(name2node)))

    logging.debug("Linking nodes...")
    for node in name2node:
       if node == "1":
           t = name2node[node]
       else:
           parent = parent2child[node]
           parent_node = name2node[parent]
           parent_node.add_child(name2node[node])
    logging.debug("Linking completed. Tree loaded.")
    
    return t


def search_for_accno(tree, accno):
    """Searches a tree for a node with node.accno attribute == accno"""
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
    accno_name_file =  "/shared/genomes/NCBI/bacterial/20140228_panindex/bacterial_genomes_map.txt"
    dumpdir = "/shared/genomes/NCBI/taxonomy/taxdump/"

    tic = time.time()
    tree= load_ncbi_tree_from_dump(dumpdir, accno_name_file)
    logging.debug("Time to read tree {}".format(time.time()-tic))

    
    print "Searching for 1280..."
    aureus = tree.search_nodes(name="1280")
    agnetis = tree.search_nodes(name="985762")
    print aureus, agnetis

    print "Searching for common ancestor to {} {}...".format(aureus[0].name, agnetis[0].name)
    ancestor = tree.get_common_ancestor([aureus[0].name, agnetis[0].name])
    print ancestor.name

    print "Searching for NC_022226"
    S_aureus = search_for_accno(tree, "NC_022226")
    print S_aureus
    print S_aureus[0].name
    print S_aureus[0].accno

    print "Searching for multiple accnos: NC_022226, NC_022222"
    multinodes = search_for_accno(tree, ["NC_022226", "NC_022222"])
    print multinodes
    for node in multinodes:
        print node.name, node.accno

