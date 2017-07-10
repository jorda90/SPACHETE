#!/bin/bash
#############################
# File Name : append_seq_emp.py
#
# Purpose : Append sequence info as a column of the appended report
#
# Creation Date : 26-04-2017
#
# Last Modified : Wed 05 Jul 2017 05:27:12 PM PDT
#
# Created By : Rob Bierman
#
##############################
import pandas as pd
import os
import re


def merge_for_appended(app_rep_path):
    #Kind of hack out the STEM and OUTPUT_DIR
    STEM = os.path.basename(app_rep_path).split('_naive_report_Appended')[0]
    OUTPUT_DIR = app_rep_path.split('/reports/')[0]

    naive_path = os.path.join(OUTPUT_DIR,'reports',STEM+'_naive_report.txt')
    glm_path = os.path.join(OUTPUT_DIR,'reports','glmReports',STEM+'_FUSION_W_ANOM_AND_INDEL_JUNCPOUT')

    appended_dir = os.path.join(OUTPUT_DIR,'reports','AppendedReports')
    if not os.path.exists(appended_dir):
        os.mkdir(appended_dir)

    appended_path = os.path.join(appended_dir,STEM+'_naive_report_Appended.txt')
    naive = pd.read_table(naive_path,sep='\t')
    glm = pd.read_table(glm_path,sep='\t')
    appended = pd.read_table(appended_path,sep='\t')

    #Rename the first naive column to match the first glm column
    naive.rename(columns={'@Junction':'junction'}, inplace=True)

    #Merge the two on their only shared column (junction)
    merged = pd.merge(naive,glm)

    #Build a dict of jct_ind to sequence
    fasta_path = os.path.join(OUTPUT_DIR,"spork_out","novel_junctions.fasta")
    ind_to_seq = {}
    jct_ind_p = r'jct_ind=(.*)'
    with open(fasta_path,'r') as fasta_file:
        while True:
            header = fasta_file.readline()
            if not header:
                break
            header = header.strip()
            seq = fasta_file.readline().strip()
            jct_ind = int(re.findall(jct_ind_p,header)[0])
            ind_to_seq[jct_ind] = seq

    #Add sequence data as an extra column to the appended report
    ordered_seqs = []
    for jct_name in merged['junction']:
        jct_ind = int(re.findall(jct_ind_p,jct_name)[0])
        ordered_seqs.append(ind_to_seq[jct_ind])
        
    merged['seq'] = ordered_seqs

    merged.to_csv(STEM+'_seq_appended_report.txt',sep='\t',index=False)


##########################################
#               Main                     #
##########################################
if __name__ == '__main__':

    #Make a txt file, where each line is the full path to the appended report output file, eg:
    #/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/CML_uconn_04_24_17/SRR3192409/reports/AppendedReports/SRR3192409_naive_report_Appended_seq_emp_p.txt
    #/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/CML_uconn_04_24_17/SRR3192410/reports/AppendedReports/SRR3192410_naive_report_Appended_seq_emp_p.txt
    app_rep_log = 'woo.txt'

    #Append seq and emp
    with open(app_rep_log,'r') as app_rep_f:
        for line in app_rep_f:
            app_rep_path = line.strip()
            new_app_rep_path = merge_for_appended(app_rep_path)

