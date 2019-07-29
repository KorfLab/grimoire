'''
Report

Takes in .eval files and find differences in predictions
'''

import argparse
import grimoire.toolbox as toolbox
import grimoire.sequence as sequence

parser = argparse.ArgumentParser(description='Quantify Evaluation')
parser.add_argument('--eval', required=True, type=str,
	metavar='<path>', help='path to input eval file (%(type)s)')
parser.add_argument('--report', required=False, type=str,
	metavar='<path>', help='path to output report file for prediction differences(%(type)s)')
parser.add_argument('--infasta', required=False, type=str,
	metavar='<path>', help='path to input FASTA file for protein translation(%(type)s)')
parser.add_argument('--outfasta', required=False, type=str,
	metavar='<path>', help='path to output FASTA file for protein translation(%(type)s)')
arg = parser.parse_args()

class ReportError(Exception):
    pass

class Eval_entry:
    def __init__(self, id, lorf, mRNA1, mRNA2):
        self.id = id
        self.lorf = lorf
        self.mRNA1 = mRNA1
        self.mRNA2 = mRNA2
        self.diff = False

#input file
lorf_mRNA1_diff = 0
lorf_mRNA2_diff = 0
mRNA1_mRNA2_diff = 0
all_diff = 0
total = 0
diff_entries = []

with open(arg.eval,'r') as eval:
    for line in eval:
        line = line.replace('\r','')
        line = line.replace('\n','')
        line = line.split('\t')
        if len(line)!=4: continue
        entry = Eval_entry(id=line[0],lorf=line[1],mRNA1=line[2],mRNA2=line[3])
        #Eval difference in predicton
        if (entry.lorf != entry.mRNA1):
            entry.diff = True
            lorf_mRNA1_diff += 1
        if (entry.lorf != entry.mRNA2):
            entry.diff = True
            lorf_mRNA2_diff += 1
        if entry.mRNA1 != entry.mRNA2:
            entry.diff = True
            mRNA1_mRNA2_diff += 1
        if entry.diff == True:
            all_diff += 1
            if arg.report:
                diff_entries.append(entry)
        total += 1

outfasta = []
if entry.diff == True:
    if arg.report:
        with open(arg.report,'w+') as report:
            for entry in diff_entries:
                report.write(entry.id+'\t'+entry.lorf+'\t'+entry.mRNA1+'\t'+entry.mRNA2+'\n')
        if not arg.infasta or not arg.outfasta:
            raise ReportError('Did not input all necessary arguments')
        fasta = toolbox.FASTA_file(arg.infasta)
        for entry in diff_entries:
            #ISSUE: THIS PART BREAKS SUPER EASILY IF ID DOES NOT 100% MATCH
            fasta_entry = fasta.get(entry.id)
            #at this point, all seq are nt
            length = len(fasta_entry.seq)
            lorf_seq = fasta_entry.seq[int(entry.lorf)-1:int(length)]
            mRNA1_seq = fasta_entry.seq[int(entry.mRNA1)-1:int(length)]
            mRNA2_seq = fasta_entry.seq[int(entry.mRNA2)-1:int(length)]
            #translate nt to pro
            lorf_seq = sequence.translate_str(lorf_seq)
            lorf_stop = lorf_seq.find('*')
            lorf_seq = lorf_seq[0:lorf_stop]
            mRNA1_seq = sequence.translate_str(mRNA1_seq)
            mRNA1_stop = mRNA1_seq.find('*')
            mRNA1_seq = mRNA1_seq[0:mRNA1_stop]
            mRNA2_seq = sequence.translate_str(mRNA2_seq)
            mRNA2_stop = mRNA2_seq.find('*')
            mRNA2_seq = mRNA2_seq[0:mRNA2_stop]

            lorf = toolbox.FASTA_entry(id=entry.id+'-lorf', desc=str(len(lorf_seq)), seq=lorf_seq)
            mRNA1 = toolbox.FASTA_entry(id=entry.id+'-m1', desc=str(len(lorf_seq)), seq=mRNA1_seq)
            mRNA2 = toolbox.FASTA_entry(id=entry.id+'-m2', desc=str(len(lorf_seq)), seq=mRNA2_seq)
            outfasta.append(lorf)
            outfasta.append(mRNA1)
            outfasta.append(mRNA2)

        with open(arg.outfasta,'w+') as out:
            for fasta_entry in outfasta:
                out.write('>'+fasta_entry.id+'\n')
                out.write(fasta_entry.seq+'\n')

print('LORF mRNA1 Diff: '+str(lorf_mRNA1_diff)+' out of '+str(total)+' predictions')
print('LORF mRNA2 Diff: '+str(lorf_mRNA2_diff)+' out of '+str(total)+' predictions')
print('mRNA1 mRNA2 Diff: '+str(mRNA1_mRNA2_diff)+' out of '+str(total)+' predictions')
print('All Diff: '+str(all_diff)+' out of '+str(total)+' predictions')
