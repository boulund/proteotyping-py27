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


class Tree:
    """Tree node (and also entire trees)"""
    def __init__(self, payload=None, children=[]):
        self.payload = payload
        self.children = children

    def __str__(self):
        return str(self.payload)

    def traverse(self, level=0):
        """Traverse the tree."""
        pass #do something with tree at level
        for subtree in self.children:
            subtree.traverse()

    def _print(self, level=0):
        """Print tree recursively:"""
        print "{char:<{level}}{payload}".format(char="", level=level, payload=self.payload)
        for subtree in self.children:
            subtree._print(level+1)



def build_test_tree():
    """Construct a simple test tree, from leaf to root"""
    leaf1 = Tree("G")
    parent1 = Tree("CC", [leaf1])
    parent2 = Tree("TTACCG")
    parent3 = Tree("A", [parent1, parent2])
    leaf2 = Tree("CTGATTCCG")
    parent4 = Tree("C", [leaf2])
    parent5 = Tree("G")
    parent6 = Tree("C", [parent4, parent5])
    tree = Tree("root", [parent3, parent6])
    return tree


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
                    ['G']
                ], 
                ['TTACCG']
            ], 
            ['C', 
                ['C', 
                    ['CTGATTACCG'], 
                ], 
                ['G']
            ], 
         ]


    tree = build_test_tree()
    tree._print()


    parents = sorted(parse_taxdump("/shared/genomes/NCBI/taxonomy/taxdump/nodes.dmp"))



