#!/usr/bin/env python

from subprocess import Popen,PIPE
import sys,re

def muscle(seqs,nohead=False):
    if nohead:
        label_seqs = list(enumerate(seqs))
    else:
        label_seqs = zip(seqs[:-1:2],seqs[1::2])

    mh = Popen(['muscle'],stdin=PIPE,stderr=PIPE,stdout=PIPE)
    mh.stdin.write('\n'.join(['>%s\n%s' % (l,s) for l,s in label_seqs]))
    alnstr = mh.communicate()[0]

    csep = re.sub(r'>(.+?)\n',r',\1,',alnstr).replace('\n','')[1:]

    if nohead:
        seqli = csep.split(',')
        return [seq for idx,seq in sorted(zip([int(i) for i in seqli[::2]],seqli[1::2]))]
    else:
        return csep

    #return ','.join(alnstr.strip().replace('\n','').split('>'))


if __name__ == '__main__':
    seqs = sys.argv[1].split(',')

    print muscle(seqs)
