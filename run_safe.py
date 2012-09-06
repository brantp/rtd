#!/usr/bin/env python

import os,sys

cmd,done = sys.argv[1:]

if os.path.exists(done):
    print >> sys.stderr, 'completion flag exists: %s' % done
else:
    ret = os.system(cmd)
    if ret == 0:
        os.system('touch %s' % done)
