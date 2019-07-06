#!/usr/bin/python3

import argparse
import base64
import io
import statistics
import sys

import matplotlib.pyplot as plt

import genome
import toolbox

## Command Line ##

parser = argparse.ArgumentParser(description='Genome annotation tool')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--title', required=True, type=str,
	metavar='<path>', help='title of report (%(type)s)')
parser.add_argument('--html', required=True, type=str,
	metavar='<path>', help='path to html output file (%(type)s)')
parser.add_argument('--json', required=False, type=str,
	metavar='<path>', help='path to json output file (%(type)s)')
parser.add_argument('--chrmap', required=False, type=str,
	metavar='<path>', help='path to file of chr name mappings (%(type)s)')
arg = parser.parse_args()

## needs a bunch of formatting options for the figures


####################
# Global Variables #
####################

FigureNumber = 0 # used for enumerating each figure in the document

ChrMap = None # used when GFF and FASTA chromosomes have different names
if arg.chrmap:
	ChrMap = {}
	with open(arg.chrmap, 'r') as file:
		for line in file:
			if line[0:1] == '#': continue
			line = line.rstrip()
			col = line.split('\t')
			if len(col) == 2:
				ChrMap[col[0]] = col[1]

########################
## Graphical Routines ##
########################

def stats(values):
	if values:
		return {
			'count' : str(len(values)),
			'min' : str(min(values)),
			'max' : str(max(values)),
			'med' : str(round(statistics.median(values))),
			'mean' : str(round(statistics.mean(values), 3)),
			'stdev' : str(round(statistics.stdev(values), 3)),
		}
	else:
		return {'count':0, 'min':'0', 'max':'0', 'med':'0', 'mean':'0', 'stdev':'0'}
		

def encode_fig(fig):
	buff = io.BytesIO()
	fig.savefig(buff, format='png', bbox_inches='tight')
	buff.seek(0)
	img = base64.b64encode(buff.read()).decode('utf-8')
	return '<br><img src="data:image/png;base64,{0}"><br>'.format(img)

def histogram(data, title, bins):
	global FigureNumber
	FigureNumber += 1
	
	fig, ax = plt.subplots()
	ax.hist(data, bins=bins)
	ax.set(title='Figure ' + str(FigureNumber) + ': ' + title)
		
	text = encode_fig(fig)
	text += '<br><b>Figure ' + str(FigureNumber) + '</b>: '
	text += str(stats(data))
	
	return text

def bar_chart(data, title):
	global FigureNumber
	FigureNumber += 1
	
	fig, ax = plt.subplots()
	names = list(data.keys())
	values = list(data.values())
	ax.barh(names, values)
	ax.set(title='Figure ' + str(FigureNumber) + ': ' + title)
	
	text = encode_fig(fig)
	text += '<b>Figure ' + str(FigureNumber) + '</b>: '
	text += 'Total: ' + str(sum(values)) + ': ' + str(data) + '.<p>'

	return text
	

###############################################################################
# Report Init
###############################################################################

doc = open(arg.html, 'w+')
doc.write('<html>')
doc.write('<head>')
doc.write('<meta charset="utf-8">')
doc.write('<title>' + arg.title + '</title>')
doc.write('</head>')
doc.write('<body>')
doc.write('<h1>' + arg.title + '</h1>')

###############################################################################
# Chromosome Report
###############################################################################

doc.write('<h2>Chromosome Report</h2>')
sizes = {}
fasta = toolbox.FASTA_stream(arg.fasta)
for entry in fasta:
	sizes[entry.id] = len(entry.seq)
doc.write('Origin: ' + arg.fasta + '<p>')
doc.write(bar_chart(sizes, 'Chromosome Sizes'))

###############################################################################
# Generic Annotation Report
###############################################################################

doc.write('<h2>Genome Annotation Report</h2>')
gff3 = toolbox.GFF_stream(arg.gff3)
types = {}
for entry in gff3:
	if entry.type not in types: types[entry.type] = 0
	types[entry.type] += 1
doc.write(bar_chart(types, 'Feature Types'))

###############################################################################
# Protein-Coding Genes Report
###############################################################################

doc.write('<h2>Protein-coding Genes Report</h2>')

tpg = [] # transcript per gene
ept = [] # exons per transcript
ipt = [] # introns per transcript
u5pt = [] # 5'utrs per transcript
u3pt = [] # 3'utrs per transcript
ilen = [] # intron lengths
elen = [] # exon lengths
u5len = [] # 5'utr lengths
u3len = [] # 3'utr lengths
tlen = [] # tx lengths
clen = [] # cds lengths

gen = genome.Genome(fasta=arg.fasta, gff3=arg.gff3, chr_map=ChrMap)
for chr in gen.chromosomes:
	for gene in chr.features:
		tpg.append(len(gene.children))
		for tx in gene.mRNAs():
				
			## count distributions
			ept.append(len(tx.exons))
			ipt.append(len(tx.introns))
			u5pt.append(len(tx.utr5s))
			u3pt.append(len(tx.utr3s))
								
			#  lengths
			for e in tx.exons: elen.append(e.end - e.beg + 1)
			for u in tx.utr5s: u5len.append(u.end - u.beg + 1)
			for u in tx.utr3s: u3len.append(u.end - u.beg + 1)
			for i in tx.introns: ilen.append(i.end - i.beg + 1)
			clen.append(len(tx.cds_str()))
			tlen.append(tx.end - tx.beg + 1)

doc.write(histogram(tpg, 'mRNAs per Gene', 10))
doc.write(histogram(ept, 'Exons per mRNA', 10))
doc.write(histogram(ipt, 'Introns per mRNA', 10))
doc.write(histogram(u5pt, 'UTR5s per mRNA', 10))
doc.write(histogram(u3pt, 'UTR3s per mRNA', 10))
doc.write(histogram(tlen, 'mRNA Lengths', 10))
doc.write(histogram(clen, 'CDS Lengths', 10))
doc.write(histogram(elen, 'Exon Lengths', 10))
doc.write(histogram(ilen, 'Intron Lengths', 10))
doc.write(histogram(u5len, 'UTR5 Lengths', 10))
doc.write(histogram(u3len, 'UTR3 Lengths', 10))

###############################################################################
# Protein-Coding Gene Errors : General
###############################################################################

doc.write('<h2>Protein-coding Gene Annotation Errors</h2>')

gene_ok, gene_err = 0, 0
tx_ok, tx_err = 0, 0
f_ok, f_err = 0, 0
gene_issues = {}
tx_issues = {}
f_issues = {}

for chr in gen.chromosomes:

	# gene level
	for gene in chr.features:
		if not gene.issues: 
			gene_ok += 1
		else:
			gene_err += 1
			for issue in gene.issues:
				if issue not in gene_issues: gene_issues[issue] = 0
				gene_issues[issue] += 1
		
		# transcript level
		for tx in gene.mRNAs():
			if not tx.issues:
				tx_ok += 1
			else:
				tx_err += 1
				for issue in tx.issues:
					if issue not in tx_issues: tx_issues[issue] = 0
					tx_issues[issue] += 1
			
			# feature level
			for f in tx.children:
				if not f.issues:
					f_ok += 1
				else:
					f_err += 1
					for issue in f.issues:
						if issue not in f.issues: f_issues[issue] = 0
						f_issues[issue] += 1

doc.write('Gene errors: {} / {} = {}'.format(gene_err, gene_ok,
	gene_err /(gene_err + gene_ok)))
doc.write(bar_chart(gene_issues, 'Gene Errors'))

doc.write('mRNA errors: {} / {} = {}'.format(tx_err, tx_ok,
	tx_err /(tx_err + tx_ok)))
doc.write(bar_chart(tx_issues, 'mRNA Errors'))

doc.write('Component errors: {} / {} = {}'.format(f_err, f_ok,
	f_err /(f_err + f_ok)))
doc.write(bar_chart(f_issues, 'Component Errors'))

###############################################################################
# Protein-Coding Gene Errors : Detailed
###############################################################################

doc.write('<h2>Detailed Error Report</h2>')
doc.write('<ul>')
for chr in gen.chromosomes:
	for gene in chr.features:
		if not gene.issues: continue
		errors = {}
		for eg in gene.issues.keys():
			if eg not in errors: errors[eg] = True 
		
		for tx in gene.mRNAs():
			if not tx.issues: continue

			for et in tx.issues.keys():
				if et not in errors: errors[et] = True
			
			for f in tx.children:
				if not f.issues: continue
				
				for ef in f.issues:
					if ef not in errors: errors[ef] = True
		doc.write('<li>' + gene.id + ': ')
		doc.write(', '.join(errors.keys()))
doc.write('</ul>')

###############################################################################
# Footer
###############################################################################

FOOTER = """
This document was prepared using <em>fathom</em>,
part of the <em>grimoire</em> tools from the KorfLab.
For more information, see
<a href="https://github.com/KorfLab/grimoire">grimoire</a> at GitHub.
"""

# probably add the date of the report to the footer

doc.write('</body>')
doc.write('<footer>' + FOOTER + '</footer>')
doc.write('</html>')
