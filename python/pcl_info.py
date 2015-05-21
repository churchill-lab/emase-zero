#! /usr/bin/env python
__author__ = 'asimons'

# A simple little program to dump the contents of a pickled, zipped,
# TranscriptHits file.

import sys, os, time
sys.path.append('/hpcdata/asimons/projects/churchill/emase')
from TranscriptHits import *

for fn in sys.argv[1:]:
    fn_name = os.path.split(fn)[1]
    fn_base = '.'.join(fn_name.split('.')[:-1])
    print 'Processing', fn_name
    th = TranscriptHits()
    print '  Reading the file'
    th.load(fn)

    # Output the strains
    ofn = fn_base + '.strains'
    print '  Writing', ofn
    with open(ofn, 'w') as f:
        for s in th.strains:    # sorted(th.strains):
            print >> f, s

    ofn = fn_base + '.reads'
    print '  Writing', ofn
    with open(ofn, 'w') as f:
        for r in th.read_names:   # sorted(th.read_names):
            print >> f, r

    ofn = fn_base + '.transcripts'
    print '  Writing', ofn
    with open(ofn, 'w') as f:
        for t in th.transcript_names:   # sorted(th.transcript_names):
            print >> f, t

    ofn = fn_base + '.map'
    print '  Writing', ofn
    th.dump_map_pretty(ofn)

    if th.gene_names:
        ofn = fn_base + '.genes'
        print '  Writing', ofn
        with open(ofn, 'w') as f:
            for g in th.gene_names:
                print >> f, g

        ofn = fn_base + '.tx_gene_mapping'
        print '  Writing', ofn
        with open(ofn, 'w') as f:
            for n in range(len(th.t_idx_to_g_idx)):
                print >> f, '{0}\t{1}\t{2}\t{3}'.format(
                    n, th.t_idx_to_g_idx[n],
                    th.transcript_names[n],
                    th.gene_names[th.t_idx_to_g_idx[n]])

    ofn = fn_base + '.apm'
    print '  Converting to AlignmentPropertyMatrix (this may take a ' \
          'while)'
    apm = th.to_AlignmentPropertyMatrix()
    print '  Writing', ofn
    apm.dump_KB_form(ofn)
    print
