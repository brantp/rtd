#!/usr/bin/env python

import os,sys
import overlap_preprocess

ifq,ofq = sys.argv[1:]

overlap_preprocess.convert_fastq(ifq,ofq)

print 'DONE'
