#!/bin/bash
#############################
# File Name : get_badfjs.py
#
# Purpose : Replace AppendNaive rept script w/ getting the badfjs here
#
# Creation Date : 12-06-2017
#
# Last Modified : Tue 13 Jun 2017 03:20:27 PM PDT
#
# Created By : Rob Bierman
#
##############################
from collections import defaultdict
import pandas as pd
import numpy as np
import glob
import ast
import os
import re

def count_indels(indel_file):
    noindel_counts = defaultdict(int)
    indel_counts = defaultdict(int)
    with open(indel_file) as h_indel:
        for line in h_indel:
            if line[0] == '@':
                continue

            jct_id,str_values = line.split('\t')
            values = ast.literal_eval(str_values) #<-- convert from str to list
            noindel_counts[jct_id] += values[len(values)/2] #<-- the mid ind is noindel
            indel_counts[jct_id] += sum(values)-noindel_counts[jct_id]
            
    return noindel_counts,indel_counts
        

def get_badfjs(sam_files):
    badfjs = defaultdict(int)
    for sam_file in sam_files:
        with open(sam_file) as h_sam:
            for line in h_sam:
                if line[0] == '@':
                    continue
                
                jct_id = line.split('\t')[0]
                badfjs[jct_id] = 1
            
    return badfjs


def post_process(OUTPUT_DIR,STEM,KNIFE_GLM):
    #Get the appended reports info
    naive_report = os.path.join(OUTPUT_DIR,'reports',STEM+'_naive_report.txt')
    glm_file = os.path.join(OUTPUT_DIR,'reports','glmReports',STEM+'_FUSION_W_ANOM_AND_INDEL_JUNCPOUT')

    #It is possible not to have a glm file produced if the input fq is too small
    glm_exists = os.path.exists(glm_file)
    naive = pd.read_table(naive_report,sep='\t')

    naive.rename(columns={'@Junction':'junction'}, inplace=True) #Rename the first naive column
    merged = naive

    if glm_exists:
        glm = pd.read_table(glm_file,sep='\t')
        merged = pd.merge(naive,glm) #Merge the two on their only shared column ('junction')

                  
    #Renaming might break if no glm
    rename_dict = {'genome-anomaly':'genome.anomaly',
                   'genome-pval':'genome.pval',
                   'reg-anomaly':'reg.anomaly',
                   'reg-pval':'reg.pval',
                   'junc-anom':'junc.anom',
                   'junc-pval':'junc.pval',
                   'FarJunc-anom':'FarJunc.anom',
                   'FarJunc-pval':'FarJunc.pval',
                   'logsum_2.x':'logsum.2.x',
                   'productPhat_lower.x':'productPhat.lower.x',
                   'junction_cdf.x':'junction.cdf.x',
                   'junction_cdf_lower.x':'junction.cdf.lower.x',
                   'logsum_2.y':'logsum.2.y',
                   'productPhat_lower.y':'productPhat.lower.y',
                   'junction_cdf.y':'junction.cdf.y',
                   'junction_cdf_lower.y':'junction.cdf.lower.y',
                   'junction_cdf_windel_diff':'junction.cdf.windel.diff'}
 
    merged.rename(columns=rename_dict, inplace=True) #Rename the merged columns

    #Split the junction column into the different components
    parse_locs = {0:'chr1',1:'gene1',2:'pos1',3:'strand1',
                  4:'chr2',5:'gene2',6:'pos2',7:'strand2',
                  8:'kind',10:'num',12:'score',14:'gap',
                  16:'don.dist',18:'acc.dist',20:'mapq',
                  22:'badfj.3',24:'jct.ind'}

    parse_vals = defaultdict(list)

    for jct in merged['junction']:
        parsed = re.split('=|:|,|\|',jct)
        for i,v in enumerate(parsed):
            if i in parse_locs:
                parse_vals[parse_locs[i]].append(v)
        
    colname_order = [x[1] for x in sorted(parse_locs.iteritems())]
    for k in colname_order:
        merged[k] = parse_vals[k]
        
    merged.replace('False',0,inplace=True)
    merged.replace('True',1,inplace=True)


    #Get the indel counts
    indels_path = os.path.join(OUTPUT_DIR,'IndelsHistogram')
    indel_file1,indel_file2 = sorted(glob.glob(os.path.join(indels_path,'*.txt')))
    noindel_counts1,indel_counts1 = count_indels(indel_file1)
    noindel_counts2,indel_counts2 = count_indels(indel_file2)
    
    merged['num.noindel.1'] = merged['junction'].map(noindel_counts1)
    merged['num.indel.1'] = merged['junction'].map(indel_counts1)
    merged['num.noindel.2'] = merged['junction'].map(noindel_counts2)
    merged['num.indel.2'] = merged['junction'].map(indel_counts2)


    #Get BadFJ1 and BadFJ2 info
    BadFJ1_sams = glob.glob(os.path.join(OUTPUT_DIR,'BadFJ',STEM,'*.sam'))
    BadFJ2_sams = glob.glob(os.path.join(OUTPUT_DIR,'BadFJ_ver2',STEM,'*.sam'))
    badfj1s = get_badfjs(BadFJ1_sams)
    badfj2s = get_badfjs(BadFJ2_sams)
    
    merged['badfj.1'] = merged['junction'].map(badfj1s)
    merged['badfj.2'] = merged['junction'].map(badfj2s)
    merged['badfj.1'].fillna(0,inplace=True)
    merged['badfj.2'].fillna(0,inplace=True)

    #Get exonL and exonR info
    #Want the maximum expressed linear juncs from the KNIFE GLM
    if KNIFE_GLM:
        gene_counts = defaultdict(int)
        with open(KNIFE_GLM) as h:
            h.readline()
            for line in h:
                jct,count,posterior = line.split('\t')[:3]
                if float(posterior) > 0.8:
                    chrom,gene1,pos1,gene2,pos2,kind,strand = re.split(':|\|',jct)
                    gene_counts[gene1] = max(gene_counts[gene1],int(count))
                    gene_counts[gene2] = max(gene_counts[gene2],int(count))

        
        merged['exon.1'] = merged['gene1'].map(gene_counts)
        merged['exon.2'] = merged['gene2'].map(gene_counts)
        merged['exon.1'].fillna(0,inplace=True)
        merged['exon.2'].fillna(0,inplace=True)

    else:
        merged['exon.1'] = -1               
        merged['exon.2'] = -1               

    #Add sequence data as an extra column to the appended report
    fasta_path = os.path.join(OUTPUT_DIR,"spork_out","novel_junctions.fasta")
    id_to_seq = {}
    with open(fasta_path,'r') as fasta_file:
        while True:
            header = fasta_file.readline()
            if not header:
                break

            jct_id = header[1:].strip()
            seq = fasta_file.readline().strip()
            id_to_seq[jct_id] = seq

    merged['seq'] = merged['junction'].map(id_to_seq)

    #Add emp_p as an extra column
    if glm_exists:
        badfj1s = merged[merged['badfj.1'] > 0]
        total_badfj1s = len(badfj1s)
        emp_ps = []

        for jct_cdf in merged['junction.cdf.y']:
            num_greater = sum(badfj1s['junction.cdf.y'] > jct_cdf)
            emp_ps.append(float(num_greater)/total_badfj1s)

        merged['emp.p'] = emp_ps


    #Write out the merged, appended, and parsed report
    appended_report = os.path.join(OUTPUT_DIR,"reports","AppendedReports",
                                   STEM+"_naive_report_Appended.txt")

    merged.to_csv(appended_report,sep='\t',index=False)

#--------------------------#
#           Main           #
#--------------------------#
if __name__ == '__main__':
    #Will be params, just testing a certain run
    OUTPUT_DIR = '/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/eng_fm_mixed_06_06_17_grch38/engs-mixed_trimmed'
    STEM = 'engs-mixed_trimmed'
    KNIFE_GLM = '/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/eng_fm_mixed_grch38_04_27_17/circReads/glmReports/engs-mixed_trimmed_R1__linearJuncProbs.txt'


