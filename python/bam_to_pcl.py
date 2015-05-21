#! /usr/bin/env python
__author__ = 'asimons'

"""
A simple routine to convert a BAM file into the compressed pcl.bz2
format.
"""

import os

import TranscriptHits as th
import time
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser('BAM Size reducer v0.1')
    p.add_argument('-o', '--output', default='',
                   help='file to write the output to.\n'
                        'If not specified, will be written to the '
                        'current directory using the basename of the'
                        'input file, minus the .bam extension, with'
                        'the new extension .pcl.bz2')
    p.add_argument('-t', '--transcripts', default='',
                   help='file containing transcript names.  If not '
                        'specified, transcripts will be harvested '
                        'from the bam file.')
    p.add_argument('-g', '--gene-mappings', default='',
                   help='file containing gene to transcript mappings. '
                        'TSV; first column is gene name, remaining '
                        'are contained transcripts')
    p.add_argument('-s', '--strains', default='A,B,C,D,E,F,G,H',
                   help='a comma-delimited string containing'
                        'the strain names used in the bam file.'
                        '(default="A,B,C,D,E,F,G,H")')
    p.add_argument('bam', help='the bam file to reduce')
    args = p.parse_args()
    return args

def main():
    args = parse_args()

    print 'Processing {0}: {1}' \
            .format(args.bam, time.strftime("%I:%M:%S"))
    if not args.output:
        base = os.path.basename(args.bam)
        if base.endswith('.bam'):
            args.output = base[:-4]
        else:
            args.output = base

    strains = [x.strip() for x in args.strains.split(',')]

    t = th.TranscriptHits(strains=strains,
                         transcript_file=args.transcripts,
                         genes_file=args.gene_mappings)
    t.from_bam(args.bam)
    t.write_file(args.output)

    print 'Done: {0}'.format(time.strftime("%I:%M:%S"))

if __name__ == '__main__':
    main()

