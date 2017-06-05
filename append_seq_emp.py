#!/bin/bash
#############################
# File Name : append_seq_emp.py
#
# Purpose : Append sequence info as a column of the appended report
#
# Creation Date : 26-04-2017
#
# Last Modified : Fri 19 May 2017 01:09:15 PM PDT
#
# Created By : Rob Bierman
#
##############################
import pandas as pd
import subprocess
import sys
import os
import re


#Just make it a function to call it many times
def add_seq_emp_columns(appended_report_path,output_path=None):
    #RB 04/26/17 add a few extra columns to the appended reports file
    #A bit messy to get the STEM
    #1) The junction consensus sequence
    #2) The emp_pvalue
    #If not given then just go ahead an overwrite the file
    if not output_path:
        output_path = appended_report_path

    OUTPUT_DIR = appended_report_path.split('/reports')[0]
    NUM_FILES = 1
    STEM_FILE = os.path.join(OUTPUT_DIR,"StemList.txt")
    jct_ind_p = r'jct_ind=(.*)'
    for index in range(NUM_FILES):
        #Get the appended report location
        STEM = None
        with open(STEM_FILE,'r') as f:
            STEM = f.readlines()[index].strip()
        appended_report = os.path.join(OUTPUT_DIR,"reports","AppendedReports",str(STEM)+"_naive_report_Appended.txt")

        empty_appended_report = True
        with open(appended_report,'r') as f:
            if len(f.readlines()) > 0:
                empty_appended_report = False
        if empty_appended_report:
            return False

        #Build a dict of jct_ind to sequence
        fasta_path = os.path.join(OUTPUT_DIR,"spork_out","novel_junctions.fasta")
        ind_to_seq = {}
        with open(fasta_path,'r') as fasta_file:
            while True:
                header = fasta_file.readline()
                if not header:
                    break
                header = header.strip()
                seq = fasta_file.readline().strip()
                jct_ind = int(re.findall(jct_ind_p,header)[0])
                ind_to_seq[jct_ind] = seq
        sys.stdout.write('\tRead in and processed sequences\n')
        sys.stdout.flush()

        #Add sequence data as an extra column to the appended report
        df = pd.read_table(appended_report,sep='\t')
        ordered_seqs = []
        for jct_name in df['@Junction']:
            jct_ind = int(re.findall(jct_ind_p,jct_name)[0])
            ordered_seqs.append(ind_to_seq[jct_ind])
        sys.stdout.write('\tRead in and processed appended report\n')
        sys.stdout.flush()

            
        df['seq'] = ordered_seqs

        #Add emp_p as an extra column
        print df.columns
        badfj1s = df[df['BadFJ=1'] == 1]
        total_num = len(badfj1s)
        emp_ps = []
        for jct_cdf in df['junction_cdf.y']:
            num_greater = sum(badfj1s['junction_cdf.y'] > jct_cdf)
            emp_ps.append(float(num_greater)/total_num)
        sys.stdout.write('\tRead in and processed emp_p\n')
        sys.stdout.flush()

        
        df['emp_p'] = emp_ps

        #Write out the new appended report (overwrite the old one)
        df.to_csv(out_path,sep='\t',index=False)
        return True


##########################################
#               Main                     #
##########################################
if __name__ == '__main__':
    app_rep_log = '/scratch/PI/horence/rob/SPACHETE_dirs/current_finished_reports.csv'

    #Append seq and emp
    if False:
        empty_app_reps = []
        non_empty_reps = []
        with open(app_rep_log,'r') as app_rep_f:
            for line in app_rep_f:
                app_rep_path, tag = line.split(',')
                out_path = app_rep_path.split('.txt')[0]+'_seq_emp_p.txt'
                sys.stdout.write(app_rep_path+'\n')
                sys.stdout.write(out_path+'\n')
                sys.stdout.flush()
                res = add_seq_emp_columns(app_rep_path,out_path)
                if not res:
                    empty_app_reps.append(app_rep_path)
                else:
                    non_empty_reps.append(app_rep_path)


        print 'Number empty:',len(empty_app_reps)
        print 'Number non-empty:',len(non_empty_reps)

    #Clean out the empty files and make a new log file
    if True:
        new_app_rep_log = '/scratch/PI/horence/rob/SPACHETE_dirs/current_finished_reports_appended.csv'
        f_out = open(new_app_rep_log,'w')

        with open(app_rep_log,'r') as f_in:
            for line in f_in:
                app_rep_path,tag = line.split(',')
                out_path = app_rep_path.split('.txt')[0]+'_seq_emp_p.txt'
                if os.path.exists(out_path):
                    f_out.write(out_path+','+tag)

        f_out.close()                    

