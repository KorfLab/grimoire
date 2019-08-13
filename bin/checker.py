#!/usr/bin/env python3

import argparse
import re
import sys

import grimoire.io as io

extended_help = """
%(prog)s is used for checking annotation file formats to ensure that
they work with grimoire.
"""

parser = argparse.ArgumentParser(
	description='Annotation file sanitizer.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--file', required=True, type=str,
	metavar='<string>', help='path to file')
parser.add_argument('--type', required=True, type=str,
	metavar='<string>', help='gff3|gtf|bed12')
arg = parser.parse_args()

if arg.type == 'gff3':
	ann = io.GFF_file(filename=arg.file)

elif arg.type == 'gtf':
	ann = io.GTF_file(filename=arg.file)

elif arg.type == 'bed12':
	ann = io.BED12_file(filename=arg.file)

print(arg.file + ' looks okay')