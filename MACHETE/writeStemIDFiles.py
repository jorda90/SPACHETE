# -*- coding: utf-8 -*-
'''
Created on Thu Feb  4 10:59:42 2016

@author: Gillian
@editor: Rob
'''

import argparse
import glob
import sys
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--origDir', required=True, help='orig directory')
parser.add_argument('-f', '--FJDir', required=True, help='Far Junctions Directory')
args = parser.parse_args()

#RB: 3/22/17 changed it to look in unaligned for denovo
FileList = glob.glob(os.path.join(args.origDir,'unaligned/forDenovoIndex/*'))
FileList = [os.path.basename(name) for name in FileList]
FileList = [name for name in FileList if 'sorted' not in name]
FileList = sorted(FileList)

#Matches anything that has the following suffixes
#(should only break on really treacherous names)
#   _1.
#   _2.
#   _r1.
#   _r2.
#   _R1.
#   _R2.
fastq_pattern = r'_[r|R]?[1|2]+\.'
FileList = sorted(FileList)
with open(os.path.join(args.FJDir,'StemList.txt'),'w') as FoutStems:
    for ind in range(len(FileList)-1):
        prefix_1,_ = re.split(fastq_pattern,FileList[ind])
        prefix_2,_ = re.split(fastq_pattern,FileList[ind+1])
        if prefix_1 == prefix_2:
            Stem = prefix_1.replace('unaligned_','')
            FoutStems.write(Stem+'\n')
        
            

