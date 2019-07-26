'''
Report

Takes in .eval files and find differences in predictions
'''

import argparse

parser = argparse.ArgumentParser(description='Quantify Evaluation')
parser.add_argument('--eval', required=True, type=str,
	metavar='<path>', help='path to eval file (%(type)s)')
parser.add_argument('--report_diff', required=False, action='store_true',
	help='give detailed report on prediction differences(%(type)s)')
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

with open(arg.eval) as eval:
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
            if arg.report_diff:
                print(entry.id+'\t'+entry.lorf+'\t'+entry.mRNA1+'\t'+entry.mRNA2)
        total += 1

print('LORF mRNA1 Diff: '+str(lorf_mRNA1_diff)+' out of '+str(total)+' predictions')
print('LORF mRNA2 Diff: '+str(lorf_mRNA2_diff)+' out of '+str(total)+' predictions')
print('mRNA1 mRNA2 Diff: '+str(mRNA1_mRNA2_diff)+' out of '+str(total)+' predictions')
print('All Diff: '+str(all_diff)+' out of '+str(total)+' predictions')
