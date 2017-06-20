#!/bin/bash
#SPACHETE call wrapper
#
#First please set the EXONS_DIR, INDEL_INDICES, and CIRC_REF which will be the same as MACHETE
#
# To run spachete you need to specify:
#(1) The path to the KNIFE named directory of interest
#(2) The path to where you want to save the output.
#    will be built if it doesn't exist
#
#Can optionally provide a list of keys to specify which files in the circpipe
#directory get run. Only samples that contain at least one of the keys get run.
#If this is not specified that every file in the circpipe directory will get run
#(3) This is the STEM_INCLUDE_LIST=("key") optional variable

#####################################
#    Example Input and outputs      #
#####################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/peter/A011_and_A012/"               #<-- Input path to the KNIFE output
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/peter_ovarian"                  #<-- Output path, will get made if does not exist
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("A012-CTTGTA_S3_L006") #<-- Optional, if not set all files will be run
#Choices for KNIFE dir
#A002-CGATGT_S1_L006/ A011_and_A012/       A016_and_A018/       HA08_and_HA09/
#A003_and_A007/       A014_and_A015/       A019_and_HA07/

#############################
#    CML UConn samples      #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/CML_UConn/circpipe_K562"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/CML_uconn_05_05_17"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("SRR3192409")

#############################
#     CML test samples      #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/CML_test/aligned/CML"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/cml_test_05_07_17"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("ENCFF000HOC2")

##########################################
#             Kami Test Set              #
##########################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/unaligned_T4_sandbox"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/kami_T4_3_22_17"
#MODE="hg19"

##########################################
#         Mouse Brain Test               #
##########################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/mouse_PE_4_3_17_2M"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/all_CDR1as/mouse_test_4_5_17"
#MODE="mm10"

##########################################
#          Human Brain Test              #
##########################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/hg19_brain_encode_04_06_17_p3"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/all_CDR1as/hg19_brain_4_7_17"
#MODE="hg19"

#############################
#       Ewing samples       #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/Ewing/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/ewing_05_08_17"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("SRR1594020"  "SRR1594021"  "SRR1594022"  "SRR1594023"  "SRR1594024"  "SRR1594025")
#STEM_INCLUDE_ONLY_LIST=("SRR1594020" "SRR1594021")

#################################
#    Normal breast samples      #
#################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_breast/circpipe"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/normal_breast_04_21_17"
#STEM_INCLUDE_ONLY_LIST=("SRR1027188" "SRR1027189" "SRR1027190")
#MODE="hg19"

################################
#    Normal fetal samples      #
################################
#KNIFE_DIR="/scratch/PI/horence/gillian/normal_fetal/circpipe_fetal"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/normal_04_21_17"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("Fetal_Adrenal_360_CTTGTA_L006" "Fetal_Intestine_392_GTGAAA_L005" "Fetal_Lung_361_CAGATC_L005" "Fetal_Adrenal_403b_GTCCGC_L008" "Fetal_Intestine_395_CGATGT_L005" "Fetal_Lung_384_CTTGTA_L007" "Fetal_Adrenal_405_GTGAAA_L008" "Fetal_Intestine_397_TGACCA_L005" "Fetal_Lung_388_AGTCAA_L007")

##########################################
#             CDR1as Sandbox             #
##########################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/CDR1as_sandbox_untrim/"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/CDR1as_sb_new_indels"
#MODE="hg19"

##########################################
#             CDR1as Full                #
##########################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/CDR1as_full"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/CDR1as_full_no_mq_3_13_17_25bp"
#MODE="hg19"

##########################################
#    FusionMap Engstrom Mixed Test       #
##########################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/eng_fm_mixed_04_25_17"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/eng_fm_mixed_06_13_17"
#MODE="hg19"
#MODE="grch38"

##########################################
#       FusionMap Positive Test          #
##########################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/hg19_fm_positive_05_04_17"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/fm_pos_06_13_17"
#MODE="hg19"
#MODE="grch38"

#####################################
#             Fetal Test            #
#####################################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/fetal_lung/fetal_lung"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/fetal_lung_test"
#MODE="hg19"

#####################################
#          CML Test                 #
#####################################
#KNIFE_DIR="/scratch/PI/horence/gillian/CML_UConn/circpipe_K562"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/CML_uconn_test_04_20_17"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("SRR3192409")
#FLANK_LEN="30"

#####################################
#          Engstrom                 #
#####################################
#KNIFE_DIR="/scratch/PI/horence/gillian/Engstrom/circpipe_engstrom"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/engstrom_4_28_17"
#MODE="hg19"

#############################
#      RNaseR samples       #
#############################
#KNIFE_DIR="/scratch/PI/horence/gillian/ov_RNaseR_Qatar/ovCircPipe"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/ovcar_05_09_17"
#MODE="hg19"
#STEM_INCLUDE_ONLY_LIST=("SRR1772257" "SRR1772957" "SRR1777309" "SRR1777310")

#############################
#    Ovarian2014 samples    #
#############################
#KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/ovarian2014_subset5"
#OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/ovarian2014_subset5_05_17_17"
#MODE="hg19"

#############################
#    Natalie ES Cells       #
#############################
KNIFE_DIR="/scratch/PI/horence/rob/KNIFE_dirs/knife_outputs/hg19_es_cells_06_14_17"
OUT_DIR="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/natalie_es_cells_06_14_17"
MODE="grch38"
STEM_INCLUDE_ONLY_LIST=("E14_S"" shZ_S")

#####################################
#                                   #
#       References Paths            #
# (should be the same across runs)  #
#                                   #
#####################################
#NOTE!!
#Have to change these paths to to same ones that MACHETE points to
#INDEL_INDICES="/scratch/PI/horence/gillian/HG19_reg_indels/toyIndelIndices/"
#INDEL_INDICES="/scratch/PI/horence/indices/indel_indices/hg19" #<-- make this more of a param
INDEL_INDICES="/scratch/PI/horence/indices/indel_indices" #<-- make this more of a param
#CIRC_REF="/share/PI/horence/circularRNApipeline_Cluster/index"
CIRC_REF="/scratch/PI/horence/indices/genomes"


#####################################
#      Dont have to change          #
#####################################
#Get the absolute path to this wrapper script
ABS_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #<-- line gets loc of this dir
TIME_STAMP=`date +"%Y_%m_%d-%H_%M"`
#OUT_DIR=$OUT_DIR"_"$TIME_STAMP #<-- uncomment to avoid write-over
if [ -z "$FLANK_LEN" ]; then
    echo "FLANK_LEN not set, setting to 25"
    FLANK_LEN="25"
fi

#Build up the OPTIONS string
OPTIONS=""
OPTIONS="$OPTIONS --root-dir $ABS_PATH"
OPTIONS="$OPTIONS --circpipe-dir $KNIFE_DIR"
OPTIONS="$OPTIONS --mode $MODE"
OPTIONS="$OPTIONS --output-dir $OUT_DIR"
#OPTIONS="$OPTIONS --hg19Exons $EXONS_DIR" #RB 3/13/17 not being used
OPTIONS="$OPTIONS --reg-indel-indices $INDEL_INDICES"
OPTIONS="$OPTIONS --circref-dir $CIRC_REF" #<-- adjust the flank len
OPTIONS="$OPTIONS --flank-len $FLANK_LEN"
OPTIONS="$OPTIONS --SLURM True" #<-- comment this line out if not on SLURM
if [ "$STEM_INCLUDE_ONLY_LIST" -a ${#STEM_INCLUDE_ONLY_LIST[@]} -gt 0 ]
then
    INCLUDE_LIST=$(printf ",%s" "${STEM_INCLUDE_ONLY_LIST[@]}")
    INCLUDE_LIST=${INCLUDE_LIST:1}
    OPTIONS="$OPTIONS --stem-include-list $INCLUDE_LIST"
fi

LOG_STEM="$ABS_PATH/logs/spachete_feeder_$TIME_STAMP" 
python "$ABS_PATH/wrappers/spachete_feeder.py" $OPTIONS 1> "$LOG_STEM.out" 2> "$LOG_STEM.err"
#python "$ABS_PATH/wrappers/spachete_feeder.py" $OPTIONS


