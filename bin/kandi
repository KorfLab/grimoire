#!/usr/bin/env python3

import argparse
import sys
import os
import curses
import curses.textpad

import grimoire.genome as genome

## Command line stuff ##

extended_help = """
%(prog)s is a very simple terminal-based genome browser for examining
fasta and GFF at arbitrary resolution. It loads all the sequences and
annotations into memory and is therefore more suitable in a
test/debugging scenario than browsing complete genomes.
"""

parser = argparse.ArgumentParser(
	description='Terminal-based genome browser.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to input fasta file')
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='path to input GFF3 file')
arg = parser.parse_args()

CHR = []
GENE = []
STATE = {
	'id':None,
	'beg':None,
	'end':None,
	'pos':None,
}

def show_help(stdscr):
	curses.curs_set(0)

	while (1):
		stdscr.clear()
		H, W = stdscr.getmaxyx()
		
		# content
		stdscr.addstr(4, 1, 'KANDI help screen')		
		stdscr.addstr(6, 1, 'General')
		stdscr.addstr(7, 5, '+ The top of the screen shows current location')
		stdscr.addstr(8, 5, '+ The bottom of the screen shows menu items')
		stdscr.addstr(10, 1, 'Selector')
		stdscr.addstr(11, 5, '+ Use the arrow keys to move the cusor')
		stdscr.addstr(12, 5, '+ Use the return key to select chromosome')
		stdscr.addstr(14, 1, 'Browser')
		stdscr.addstr(15, 5, '+ Use the arrow keys to move the cusor')
		stdscr.addstr(16, 5, '+ Use the +/- keys to zoom')
		stdscr.addstr(17, 5, '+ Add the shift key to make bigger steps')
		stdscr.addstr(18, 5, '+ Use the return key to center the screen')
	
		# location & menu
		loc = '{}:{}-{} {}'.format(
			STATE['id'],STATE['beg'], STATE['end'], STATE['pos'])
		stdscr.addstr(0, 0, loc, curses.color_pair(1))
		stdscr.addstr(H-1, 0, '(s)elector (v)iewer (q)uit',
			curses.color_pair(1))

		# action
		k = stdscr.getch()
		if   k == ord('s'): show_selector(stdscr)
		elif k == ord('v'): show_viewer(stdscr)
		elif k == ord('q'): sys.exit(0)
		
		stdscr.refresh()

def show_selector(stdscr):
	curses.curs_set(2)
	top = 4
	bot = top + len(CHR) -1
	cx = 3
	cy = top
	stdscr.move(cy, cx)
	k = None

	while (1):
		stdscr.clear()
		H, W = stdscr.getmaxyx()
		
		# content
		stdscr.addstr(2, 1, 'KANDI chromosome selector')	
		for i in range(len(CHR)):
			stdscr.addstr(4 + i, 4, CHR[i].name)
	
		# location & menu
		loc = '{}:{}-{} {}'.format(
			STATE['id'],STATE['beg'], STATE['end'], STATE['pos'])
		stdscr.addstr(0, 0, loc, curses.color_pair(1))
		stdscr.addstr(H-1, 0, '(b)rowser (h)elp (q)uit',
			curses.color_pair(1))

		# action
		if   k == ord('h'): show_help(stdscr)
		elif k == ord('v'): show_viewer(stdscr)
		elif k == ord('q'): sys.exit(0)
		elif k == 258: cy = cy + 1
		elif k == 259: cy = cy - 1
		elif k == 10:
			reset_chrom(cy - top)
			show_viewer(stdscr)
		
		if cy < top: cy = top
		if cy > bot: cy = bot
		stdscr.move(cy, cx)
		stdscr.refresh()
		k = stdscr.getch()
		
def organize_genes(W):
	scale = len(STATE['seq']) / W 



def show_viewer(stdscr):
	curses.curs_set(2)
	H, W = stdscr.getmaxyx()
	cx = int(W/2)
	cy = int(H/2)
	stdscr.move(cy, cx)
		
	k = None
	while (1):
		stdscr.clear()
		H, W = stdscr.getmaxyx()
		scale = len(STATE['seq']) / W # scale
		
		mat = []
		max_rows = H-2
		for y in range(max_rows):
			mat.append([])
			for x in range(W):
				mat[y].append([])
				mat[y][x] = None

		for gene in STATE['genes']:
			c = 5 # color
			if gene.strand == '-': c = 7
			beg, end = None, None
			if gene.beg > STATE['beg'] and gene.beg < STATE['end']:
				beg = int(gene.beg / scale)
			if gene.end > STATE['beg'] and gene.end < STATE['end']:
				end = int(gene.end / scale)
			if not (beg or end): continue

			for y in range(max_rows):
				if beg == end:
					if mat[y][beg] == None:
						mat[y][beg] = gene
						stdscr.addstr(y+2, beg, '|', curses.color_pair(c))
						break
				elif beg and end:
					clear = True
					for x in range (beg, end+1):
						if mat[y][x] != None:
							clear = False
							break
					if not clear: continue
					for x in range (beg, end+1): mat[y][x] = gene
					stdscr.addstr(y+2, beg, '[', curses.color_pair(c))
					for i in range (beg+1, end):
						stdscr.addstr(y+2, i, '=', curses.color_pair(c))
					stdscr.addstr(y+2, end, ']', curses.color_pair(c))
					break
				elif beg:
					pass
				elif end:
					pass
			
		# location & menu
		loc = '{}:{}-{} {}'.format(
			STATE['id'],STATE['beg'], STATE['end'], STATE['pos'])
		stdscr.addstr(0, 0, loc, curses.color_pair(1))
		stdscr.addstr(H-1, 0, '(s)elector (h)elp (q)uit',
			curses.color_pair(1))

		# action
		if   k == ord('h'): show_help(stdscr)
		elif k == ord('s'): show_selector(stdscr)
		elif k == ord('q'): sys.exit(0)
		elif k == 258: cy = cy + 1    # arrow down
		elif k == 259: cy = cy - 1    # arrow up
		elif k == 261: cx = cx + 1    # arrow right
		elif k == 402: cx += int(W/3) # shift-arrow right
		elif k == 260: cx = cx - 1    # arrow left
		elif k == 393: cx -= int(W/3) # shift-arrow left
		elif k == 10:  cx = int(W/2)  # re-center cursor
		elif k == ord('-'): zoom(-1)
		elif k == ord('_'): zoom(-5)
		elif k == ord('='): zoom(+1)
		elif k == ord('+'): zoom(+5)
		
		if cx < 0:   cx = 0
		if cx > W-1: cx = W-2
		if cy < 1:   cy = 1
		if cy > H-2: cy = H-2
		
		stdscr.move(cy, cx)
		stdscr.refresh()
		k = stdscr.getch()



def browser(stdscr):

	# initialize curses	and begin on the help page
	curses.start_color()
	curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)
	curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
	curses.init_pair(3, curses.COLOR_YELLOW, curses.COLOR_BLACK)
	curses.init_pair(4, curses.COLOR_GREEN, curses.COLOR_BLACK)
	curses.init_pair(5, curses.COLOR_CYAN, curses.COLOR_BLACK)
	curses.init_pair(6, curses.COLOR_BLUE, curses.COLOR_BLACK)
	curses.init_pair(7, curses.COLOR_MAGENTA, curses.COLOR_BLACK)
	show_help(stdscr)

def reset_chrom(idx):

	# when a new chromosome is chosen, reset the zoom and cursor
	STATE['id'] = CHR[idx].name
	STATE['beg'] = 1
	STATE['end'] = len(CHR[idx].seq)
	STATE['pos'] = int(len(CHR[idx].seq) / 2)
	STATE['seq'] = CHR[idx].seq
	STATE['features'] = CHR[idx].ftable.features
	STATE['genes'] = GENE[idx]
	
if __name__ == '__main__':

	# read sequence and annotation into memory
	gen = genome.Reader(fasta=arg.fasta, gff=arg.gff3)
	for chr in gen:
		CHR.append(chr)
		genes = chr.ftable.build_genes()
		GENE.append(genes)
		
	# set the default chromosome to the first in the file
	reset_chrom(0)

	# start the curses session
	curses.wrapper(browser)

