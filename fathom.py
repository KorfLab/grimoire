#!/usr/bin/python3

import argparse
import base64
import io
import statistics

import matplotlib.pyplot as plt

import genome

## Command Line ##

parser = argparse.ArgumentParser(description='Genome annotation tool')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--title', required=True, type=str,
	metavar='<path>', help='title of report (%(type)s)')
parser.add_argument('--out', required=True, type=str,
	metavar='<path>', help='path to output file (%(type)s)')
arg = parser.parse_args()

## needs a bunch of formatting options for the figures


####################
# Global Variables #
####################

FigureNumber = 0 # used for enumerating each figure in the document

########################
## Graphical Routines ##
########################

def histogram(data, title, bins):
	global FigureNumber
	FigureNumber += 1
	
	fig, ax = plt.subplots()
	ax.hist(data, bins=bins)
	ax.set(title='Figure ' + str(FigureNumber) + ': ' + title)
	buff = io.BytesIO()
	fig.savefig(buff, format='png', bbox_inches='tight')
	buff.seek(0)
	img = base64.b64encode(buff.read()).decode('utf-8')
	text = '<img src="data:image/png;base64,{0}">'.format(img)
	text += '<br><b>Figure ' + str(FigureNumber) + '</b>: '
	text += 'min '    + str(min(data)) + ', '
	text += 'max '    + str(max(data)) + ', '
	text += 'mean '   + str(round(statistics.mean(data), 3)) + ', '
	text += 'stdev '  + str(round(statistics.stdev(data), 3)) + ', '
	text += 'median ' + str(statistics.median(data)) + '.<p>'
	return text

###############################################################################
# Header (might want to add some CSS)
###############################################################################

doc = open(arg.out, 'w+')
doc.write('<html>')
doc.write('<head>')
doc.write('<meta charset="utf-8">')
doc.write('<title>' + arg.title + '</title>')
doc.write('</head>')

###############################################################################
# Body
###############################################################################

doc.write('<body>')
doc.write('<h1>' + arg.title + '</h1>')

doc.write('<h2>Genome Characteristics</h2>')
doc.write('Number and sizes of chromosomes, percent N, filename, data, etc')


doc.write('<h2>Genome Annotation Report</h2>')

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

doc.write('Number and type of features')


doc.write('<h2>Protein-coding Genes Report</h2>')
gen = genome.Genome(fasta=arg.fasta, gff3=arg.gff3)
for chr in gen.chromosomes:
	for gene in chr.features:
		tpg.append(len(gene.children))
		for tx in gene.transcripts():
				
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


doc.write('<h2>Protein-coding Gene Annotation Errors</h2>')


doc.write('<h2>Detailed Error Report</h2>')


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
