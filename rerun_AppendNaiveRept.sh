#!/bin/sh

#  AppendNaiveRept.sh
#  
#
#  Created by Gillian Hsieh on 2/1/16.
#
### The AppendNaiveRept.sh shell calls the AppendNaiveRept.py script.  This reads in the IndelsHistogram, BadFJ and BadFJ_ver2 files, and GLM report results and outputs all the results into a single file in /FJDir/reports/AppendedReports/<STEM>_naive_report_Appended.txt
FarJuncDir=${1}
GLMReportDir=${2}
INSTALLDIR=${3}
FJGLMReportsDir=${4}
if [ $# -ge 5 ]
then
OUTPUTDIR="-o ${5}"
fi

STEMFILE=${1}StemList.txt
STEM=`awk 'FNR == '${SLURM_ARRAY_TASK_ID}' {print $1}' ${STEMFILE}`
echo ""
echo "In rerun_AppendNaiveRept.sh: FARJUNCDIR=$FarJuncDir"
echo "In rerun_AppendNaiveRept.sh: GLMReportDir=$GLMReportDir"
echo "In rerun_AppendNaiveRept.sh: INSTALLDIR=$INSTALLDIR"
echo "In rerun_AppendNaiveRept.sh: FJGLMReportsDir=$FJGLMReportsDir"
echo ""
ml load python/2.7.5

python ${INSTALLDIR}rerun_AppendNaiveRept.py -f ${1} -g ${2} -s ${STEM} -G ${4} ${OUTPUTDIR}

echo "rerun_AppendNaiveReports.sh completed for ${STEM} -- check ${1}reports/AppendedReports/${STEM}_naive_report_Appended.txt.    This is the last step of the MACHETE."  >> ${1}MasterError.txt
