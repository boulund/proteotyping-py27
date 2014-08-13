#!/usr/bin/env python
# Fredrik Boulund 
# 2014
# Taxonomic trees

from sys import argv, exit

#if len(argv)<3:
#    print "usage: script.py TAXDUMP OUTFILE"
#    exit()
#else:
#    taxdump = argv[1]
#    outfile = argv[2]


def traverse(tree, level=0):
    """Traverse tree."""
    pass #do something with tree at level
    for subtree in tree[1:]:
        traverse(subtree, level+1)


def print_tree(tree, level=0):
    """Pretty print a tree."""
    print "{char:<{level}}{node}".format(char="", level=level, node=tree[0])
    for subtree in tree[1:]:
        print_tree(subtree, level+1)


def parse_taxdump(taxdump):
    """Parse NCBI Taxonomy dumpfile (nodes.dmp)"""
    with open(taxdump) as f:
        parents = {}
        for line in f:
            node_id, parent_id =  line.split("\t|\t")[0:2]
            parents[node_id] = parent_id

    return parents

if __name__ == "__main__":
    tt = ['root', 
            ['A', 
                ['CC', 
                    ['CCTGATTACCG'], 
                    ['G']
                ], 
                ['TTACCG']
            ], 
            ['C', 
                ['C', 
                    ['CTGATTACCG'], 
                    ['TGATTACCG'], 
                    ['G']
                ], 
                ['TGATTACCG'], 
                ['G']
            ], 
            ['T', 
                ['GATTACCG'], 
                ['TACCG'], 
                ['ACCG']
            ], 
            ['GATTACCG']
         ]

    #print_tree(tt)
    parents = sorted(parse_taxdump("/shared/genomes/NCBI/taxonomy/taxdump/nodes.dmp"))

    print parents


