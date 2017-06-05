#!/bin/bash
#############################
# File Name : make_mm10_indels.sbatch
#
# Purpose : Make mm10 indel indices for SPACHETE
#
# Creation Date : 20-03-2017
#
# Last Modified : Mon 20 Mar 2017 11:26:09 AM PDT
#
# Created By : Rob Bierman
#
##############################

#########################################################
#                                                       #
#      Parameters that should be set every time         #
#                                                       #
#########################################################
#Single commented out SBATCH are active, more than that is inactive
#SBATCH -p horence                          # which node to run on
#SBATCH --time=48:00:00                     # how much time to run
#SBATCH --mem=40000                         # how much mem in MBs
#SBATCH --job-name=mm10_indels


#########################################################
#                                                       #
#               More optional parameters                #
#                                                       #
#########################################################
##SBATCH --output=slurm.out                 # slurm.out file
##SBATCH --error=slurm.err                  # slurm.err file
##SBATCH --mail-type=END,FAIL               # notifications for job done & fail
##SBATCH --mail-user=rbierman@stanford.edu  # email to send it to

########################################
#  Directions from Gillian over email  #
########################################
#The program is 
#/scratch/PI/horence/gillian/MACHETE/AddIndelsToFasta.py

#The fields are 
#parser.add_argument("-i", "--infile", required=True, help = "input file")
#parser.add_argument("-o", "--outDir", required=True, help = "output directory")
#parser.add_argument("-s", "--stem", required = True, help = "stem name" )
#parser.add_argument("-n", "--MaxInDel", required=True, help = "# of indels on each side to test")

#You would use the following
#-i /share/PI/horence/circularRNApipeline_Cluster/index/mm10_junctions_reg.fa -o <whatever you want your output dir to be called> -s <"stem" for me means the name of the experiment> -n 5 (I assume you still want 5 indels)

#The output is going to be 
#<outputdir>/<stem>_FJ_indels_1.fa, ..._2.fa, ..._3.fa, ..._4.fa, and ..._5.fa

#########################################
#  Calling the directions from Gillian  #
#########################################
exec_dir="/scratch/PI/horence/gillian/MACHETE/"
mm10_jcts="/share/PI/horence/circularRNApipeline_Cluster/index/mm10_junctions_reg.fa"
out_dir="/scratch/PI/horence/rob/SPACHETE_dirs/indel_indices/mm10_indel_inds/"
stem="mm10"
indels="5"

cd $exec_dir
python AddIndelsToFasta.py -i $mm10_jcts -o $out_dir -s $stem -n $indels


