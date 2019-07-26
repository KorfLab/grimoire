'''
Report

Takes in .eval files and find differences in predictions
'''

import argparse

parser = argparse.ArgumentParser(description='Quantify Evaluation')
parser.add_argument('--eval', required=True, type=str,
	metavar='<path>', help='path to input eval file (%(type)s)')
parser.add_argument('--report', required=False, type=str,
	metavar='<path>', help='path to output report file for prediction differences(%(type)s)')
arg = parser.parse_args()

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
report_entries = []

with open(arg.eval,'r') as eval:
    for line in eval:
        line = line.replace('\r','')
        line = line.replace('\n','')
        line = line.split('\t')
        if len(line)<4: continue
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
                report_entries.append(entry)
        total += 1

if entry.diff == True:
    with open(arg.report,'w+') as report:
        for entry in report_entries:
            report.write(entry.id+'\t'+entry.lorf+'\t'+entry.mRNA1+'\t'+entry.mRNA2+'\n\r')

print('LORF mRNA1 Diff: '+str(lorf_mRNA1_diff)+' out of '+str(total)+' predictions')
print('LORF mRNA2 Diff: '+str(lorf_mRNA2_diff)+' out of '+str(total)+' predictions')
print('mRNA1 mRNA2 Diff: '+str(mRNA1_mRNA2_diff)+' out of '+str(total)+' predictions')
print('All Diff: '+str(all_diff)+' out of '+str(total)+' predictions')
