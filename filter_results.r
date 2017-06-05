#!/bin/bash
library(tidyr) #<-- used to split the jct into chr1, pos1, etc
require(data.table)

report_stem <- '/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/ovarian2014_05_15_17'

#############################################
#     Change the file to actually run on    #
#############################################
#report_suffix <- 'ovarian2014_ACAGTG/reports/AppendedReports/ovarian2014_ACAGTG_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_ACTTGA/reports/AppendedReports/ovarian2014_ACTTGA_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_ATCACG/reports/AppendedReports/ovarian2014_ATCACG_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_CAGATC/reports/AppendedReports/ovarian2014_CAGATC_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_CGATGT/reports/AppendedReports/ovarian2014_CGATGT_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_CTTGTA/reports/AppendedReports/ovarian2014_CTTGTA_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_GATCAG/reports/AppendedReports/ovarian2014_GATCAG_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_GCCAAT/reports/AppendedReports/ovarian2014_GCCAAT_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_GGCTAC/reports/AppendedReports/ovarian2014_GGCTAC_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_TAGCTT/reports/AppendedReports/ovarian2014_TAGCTT_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_TGACCA/reports/AppendedReports/ovarian2014_TGACCA_naive_report_Appended.txt'
#report_suffix <- 'ovarian2014_TTAGGC/reports/AppendedReports/ovarian2014_TTAGGC_naive_report_Appended.txt'
report_suffix <- 'ovarian2014_unmatched/reports/AppendedReports/ovarian2014_unmatched_naive_report_Appended.txt'
stem <- strsplit(report_suffix,split='/')[[1]][1]
out_stem <- paste(report_stem,stem,sep='/')
report_path <- paste(report_stem,report_suffix,sep='/')

##########################################
#                                        #
#       Read in the spachete jcts        #
#                                        #
##########################################
sj <- read.table(report_path,sep='\t',header=TRUE,fill=TRUE)
sj <- na.omit(sj) #<-- remove rows that have na's in them (not sure why they do)

#Parse the spacehte jct header
#chr5:GPR98:90358447:+|chr5:GPR98:90373534:+|no_fusion-niether_local-intrachrom_plus_reg,num=2,score=0.0,gap=-10,don-dist:-9974,acc-dist:5270,mapq=0,badfj3:False,jct_ind=43072
sj['X.Junction'] <- lapply(sj['X.Junction'],function(x) { gsub('\\|',':',x)})
sj['X.Junction'] <- lapply(sj['X.Junction'],function(x) { gsub('=',':',x)})
sj['X.Junction'] <- lapply(sj['X.Junction'],function(x) { gsub(',',':',x)})
groups <- c('chr1','gene1','pos1','strand1','chr2','gene2','pos2','strand2','fusion_type','num=','num','score=','score','gap=','gap','don_dist=','don_dist','acc_dist=','acc_dist','mapq=','mapq','badfj3=','badfj3','jct_ind=','jct_ind')
drop_groups <- c('num=','score=','gap=','don_dist=','acc_dist=','mapq=','badfj3=','jct_ind=')
sj <- separate(data=sj,col='X.Junction',into=groups, sep=':')
sj <- sj[,!(colnames(sj) %in% drop_groups)]
sj <- data.table(sj)

#Convert all these columns to numeric (probably a better way to do this)
sj$pos1 <- as.numeric(as.character(sj$pos1))
sj$pos2 <- as.numeric(as.character(sj$pos2))
sj$num <- as.numeric(as.character(sj$num))
sj$score <- as.numeric(as.character(sj$score))
sj$don_dist <- as.numeric(as.character(sj$don_dist))
sj$acc_dist <- as.numeric(as.character(sj$acc_dist))
sj$mapq <- as.numeric(as.character(sj$mapq))
sj$jct_ind <- as.numeric(as.character(sj$jct_ind))
sj$NetPValue <- as.numeric(as.character(sj$NetPValue))
sj$junction_cdf.y <- as.numeric(as.character(sj$junction_cdf.y))
sj$numReads.y <- as.numeric(as.character(sj$numReads.y))
sj$junction_cdf_lower.y <- as.numeric(as.character(sj$junction_cdf_lower.y))
sj$productPhat.y <- as.numeric(as.character(sj$productPhat.y))

##########################################
#                                        #
#       Filter out bad spachete jcts     #
#                                        #
##########################################
print('Initial spach_jcts len:')
print(nrow(sj))

##Replace dash with 0 and get rid of rows with missing data
#sj = sj.replace(r'^-$','0',regex=True)

#Compute the empirical p values
badfj1s <- sj[BadFJ.1 == 1,]
badfjs <- badfj1s[chr1 != chr2 | abs(pos2-pos1) > 1e6 | strand1 != strand2]
badfj_cdfs <- badfjs$junction_cdf.y
print(length(badfj_cdfs))
get_emp_p <- function(x){ return(sum(badfj_cdfs > x)/length(badfj_cdfs)) }
sj$emp_p <- sapply(sj$junction_cdf.y,get_emp_p)

#Filter based on statistical score
p.thresh <- 0.1
min.reads <- 1
min.cdf.y <- 0.2
min_lower_cdf <- 0.5
min.phat <- 0.5

#Filter
sj <- sj[((emp_p<p.thresh & junction_cdf.y>min.cdf.y & numReads.y>min.reads) |
          (numReads.y==min.reads & junction_cdf_lower.y>min.phat & productPhat.y>min.phat))]
         #(emp_p<p.thresh & numReads.y>min.reads)]

print('Len after filtering on stats:')
print(nrow(sj))

#Check how many fusions there were
fusions <- sj[chr1 != chr2]
print('Num fusions:')
print(nrow(fusions))

#Filter on badfj1, badfj2, and badfj3
len_prior <- nrow(sj)

###Different filtering options, worked great w/ badfj3 == "True"
###Not having it increases the FPs a lot
#sj <- sj[BadFJ.1 == 0 & BadFJv2.1 == 1 & badfj3 == "True"]
#sj <- sj[BadFJ.1 == 0 & badfj3 == "True"]
#sj <- sj[BadFJv2.1 == 1] #<-- not good enough
sj <- sj[badfj3 == "True"] #<-- works best right now

print('There were BadFJ that made it through stat filtering')
print(len_prior-nrow(sj))

#Sort the rows
sj$diff_genes <- sj$gene1 != sj$gene2
sj <- sj[with(sj,order(-diff_genes,emp_p,junction_cdf.y)), ]

#Rearrange the columns
sj <- as.data.frame(sj) #<-- reordering columns is different for data table
front_cols <- c('diff_genes','emp_p','junction_cdf.y','junction_cdf_lower.y')
potential_front_cols <- c('chr1','gene1','pos1','strand1','chr2','gene2','pos2','strand2')
col_order <- c(front_cols,setdiff(names(sj),front_cols))
sj <- sj[col_order]

write.table(sj,paste(out_stem,'jcts_passed.tab',sep='_'),sep='\t',row.names=FALSE)

