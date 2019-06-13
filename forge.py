#!/usr/bin/python3

import argparse

from toolbox.fasta import Fasta
from toolbox.gff import Gff

## Command line stuff ##

parser = argparse.ArgumentParser(description='Something')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--target', required=True, type=str,
	metavar='<path>', help='path to target directory (%(type)s)')
parser.add_argument('--split', required=False, type=int,
	metavar='<int>', help='split data into <int> parts (%(type)s)')
parser.add_argument('--stub', required=False, type=str,
	metavar='<str>', help='filename stub for use with --split (%(type)s)')
parser.add_argument('--train', required=False, type=str,
	metavar='<type>', help='prok|euk|ime (%(type)s)')
args = parser.parse_args()

## Splitting ##

if args.split:
	ff = Fasta(args.fasta)
	gf = Gff(args.gff)

"""
for id in ff.ids:
	entry = ff.get(id)
	print(entry.id, entry.desc, entry.seq[0:50])


print(gf.types)
print(gf.chroms)
for f in gf.get(chrom='Chr1', beg=7500, end=8500):
	print(f.beg, f.end, f.type, f.strand)
"""
