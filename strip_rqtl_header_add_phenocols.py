#!/usr/bin/env python

'''given a source map and a new target map, emits only the data lines from target map, with column padding to match source map (padding up to ID column in source)
'''

import sys

source_map_f, new_map_f = sys.argv[1:]

smapfh = open(source_map_f)
nmapfh = open(new_map_f)

header = smapfh.readline().split(',')

for l in nmapfh.readlines()[3:]:
    print (',' * header.index('ID') ) + l.strip()
    

