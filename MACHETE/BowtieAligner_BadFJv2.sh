#!/bin/bash -eu

#  BowtieAligner.batch.sh
#  
#
#  Created by Gillian Hsieh on 11/10/15.
#

#RB 4/4/17 I have to get rid of this shell script clutter

BOWTIEPARAM=${1}

#module load bowtie/2.2.4


bowtie2 ${1}
