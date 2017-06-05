#/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/engstrom_4_4_17/Human_simulated_reads1

library(binom)

## simulated data anlysis
par(mfrow=c(2,4))
rm(list=ls())
i=0

require(data.table)
## eventually want one of these loops
#for (name in c("fm_engstrom","normal_fetal","ov_RNaseR_Qatar","Engstrom","Ewing","normal_breast" ,"CML_test","CML_UConn")){


for (name in "Engstrom"){ #Unfiltered

for (spork in c(0,1)){# should be 0,1

## ROB, pls check these paths/confirm?

if (spork==1){

if (name %like% "CML_U"){ name ="CML_uconn_500K_chromcheck_11_18"}
if (name %like% "CML_test"){ name ="CML_test_9_21"}
if (name %like% "Ewing"){ name ="ewing_9_14"}
if (name %like% "fetal"){ name ="normal_fetal_500K_11_19"}
if (name %like% "breast"){ name ="normal_breast_9_21"}
if (name %like% "Engstrom"){name="engstrom_9_9" ## rob pls change this name/path so its consistent w/ your current nomenclature
name="/engstrom_4_4_17/"#Human_simulated_reads1/"
#/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/engstrom_4_4_17/Human_simulated_reads1
}



if (name %like% "ov_RNaseR_Qatar"){name="ovcar3_9_15"}
}

if (!(name %like% "fm_engstrom")){
dir=paste("/scratch/PI/horence/gillian/All_AppendedReports_Jun20/",name,"/",sep="")
}
if ((name %like% "fm_engstrom")){
dir=paste("/scratch/PI/horence/rob/MACHETE_dirs/machete_outputs/fm_engstrom/reports/AppendedReports/")
}
# for spachete
if (spork==1){

dir=paste("/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/",name,"/",sep="")
}
print ("dir")
print(dir)
require(data.table)

for (tmyfile in list.files(dir)){
myfile=tmyfile

if (spork==1){
myfiledir1=paste(dir,tmyfile,"/reports/AppendedReports/",sep="")

myfile=c(list.files(myfiledir1))
}

print(paste("myfile is ",myfile))

print ("ROB, I HARCODED THE NEXT LINE- let's pls make it programatic")
myfile="/scratch/PI/horence/rob/SPACHETE_dirs/spachete_outputs/engstrom_4_4_17/Human_simulated_reads1/reports/AppendedReports/Human_simulated_reads1_naive_report_Appended.txt"


if (length(myfile)>0){

if ((myfile %like% "Appended" & myfile %like% "naive" )){

print (paste(myfile))

## diff directory structure for spork
if (!(spork==1)){
if (!is.null(tryCatch(read.delim(paste(dir,myfile,sep=""),sep="\t"), error=function(e) NULL))){
print (myfile)
m=data.table(read.delim(paste(dir,myfile,sep=""),sep="\t"))
}
}
## diff directory structure for spork
if (spork==1){
if (!is.null(tryCatch(read.delim(paste(myfiledir1,myfile,sep=""),sep="\t"), error=function(e) NULL))){
print (myfile)
m=data.table(read.delim(paste(myfiledir1,myfile,sep=""),sep="\t"))
}
}

m[,sample:=paste("spork_",spork,"_",myfile,sep="")]
m[,sampleType:=name ]
g=data.frame(m)
names(m)[1]="junction"
if (dim(m)[2]>17){
names(m)[19]="badfj1is1"
names(m)[20]="badfj2is1"

if ((is.null(match("productPhat.y",names(m))))){
setnames(m,"p_predicted.x","productPhat.x")
setnames(m,"p_predicted.y","productPhat.y")
setnames(m,"p_value.x","junction_cdf.x")
setnames(m,"p_value.y","junction_cdf.y")
}
g=data.frame(m[ productPhat.y !="-" & productPhat.x!="-",])
print (head(g))
g$junction_cdf.y=as.numeric(as.vector(g$junction_cdf.y))
g$productPhat.y=as.numeric(as.vector(g$productPhat.y))
g$numReads.y=as.numeric(as.vector(g$numReads.y))
g$junction_cdf.x=as.numeric(as.vector(g$junction_cdf.x))
g$productPhat.x=as.numeric(as.vector(g$productPhat.x))
}

print ("finished ")
print (name)
if (i==0){allg=g}
if (i>0){allg=rbind(g,allg)}
i=i+1
print (head(allg))
}
}
}
}
}
print("finishedloop")
allg=data.table(allg)

allg[,junction:=gsub("([|])", ":", paste(junction))]
# example:     chr19:DAZAP1:1432689:+:chr22:SEPT5:19708072:+:fusion      0

allg[,numF:= as.character(lapply(strsplit(paste(allg$junction), split=","), "[", 2))]
allg[,numIN:= as.character(lapply(strsplit(paste(allg$numF), split="="), "[", 2))]
allg[,scoreIN:= as.character(lapply(strsplit(paste(allg$junction), split=","), "[", 3))]
allg[,score:= as.character(lapply(strsplit(paste(allg$scoreIN), split="="), "[", 2))]

allg[,chr1:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 1))]
allg[,gene1:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 2))]
allg[,pos1:= as.numeric(as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 3)))]
allg[,strand1:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 4))]

allg[,chr2:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 5))]
allg[,gene2:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 6))]
allg[,pos2:= as.numeric(as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 7)))]
allg[,strand2:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 8))]

allg[,type:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 9))]

## star used:
#g= g[SpliceType %like% "ONLY" & (chr1!=chr2 | abs(pos1-pos2)>1000000 | strand1!=strand2 )]
allg[,anomSum:=genome.anomaly+FarJunc.anom+reg.anomaly+junc.anom+NoPartner+unaligned]
epsilon=.1
allg[,AnomPhat:=((epsilon+sum(anomSum))/(epsilon+sum(numReads.y)+sum(anomSum))), by=list(sample,junction)]
## make sure that the AnomPhat (in this case) is the same across the listed variables
allg[,AnomSD:=sqrt((1/(sum(numReads.y+anomSum)))*AnomPhat*(1-AnomPhat)), by=list(sample,junction)]
allg[,upperCI:=AnomPhat + 2*AnomSD]

p_pred_thresh=0
p_val_thresh=.2

##MOD

## informative plots and calling thresholds
i=0

par(mfrow=c(2,2))

for (my.sample in unique (allg$sample)){

train.good=allg [ paste(sample)==my.sample & (badfj1is1+badfj2is1)==0]

# these should be junctions that are bad whereas fj2is bad are more likley to just be linear artifacts

if (spork==1){
# cryptic exons may be detected by spork and also badfj1
train.bad=allg[(sample %like% spork) &(!(junction %like% "no_fusion")) & paste(sample)==my.sample & badfj1is1>0 & badfj2is1==0 ]
train.bad=allg[(sample %like% "spork_0") & paste(sample)==my.sample & badfj1is1>0 ]
}
if (spork==0){
train.bad=allg[(sample %like% spork) & paste(sample)==my.sample & badfj1is1>0 ]
}

if (  dim(train.bad) [1]>0 ){
plot(train.bad[sample==my.sample]$junction_cdf.y, train.bad[sample==my.sample]$numReads.y, main=paste(unique(train.bad[sample==my.sample,sampleType]),"badFJs"))
}
if (  dim(train.bad) [1] == 0 ){
plot(c(0,0), main=paste("NODATA", unique(allg[sample==my.sample,sampleType]),"badFJs"))
}

if (  dim(train.good) [1] == 0 ){
plot(c(0,0), main=paste("NODATA", unique(allg[sample==my.sample,sampleType]),"goodFJs"))
}
if (  dim(train.good) [1] > 0 ){

plot(train.good[sample==my.sample]$junction_cdf.y, train.good[sample==my.sample]$numReads.y, main=paste(unique(train.good[sample==my.sample,sampleType]),"goodFJs", my.sample))
}
INC=100

train.good[, discrete_p:=round(INC*junction_cdf.y)/INC]
for (th in c(0:INC)/INC){
## assign a threshold
tot.bad=length(train.bad$junction_cdf.y)
tot.good=length(train.good$junction_cdf.y)
prop=tot.bad/tot.good
print ("adding proportion of good and bad-- CHANGE")
current.p=sum(train.bad$junction_cdf.y>th)/tot.bad

print(paste("theshold is",th,"current ",current.p))			      


train.good[discrete_p==th,emp_p:=current.p]

print(train.good[discrete_p==th,])
print (paste("current thresh and fdr",th,current.p))
}

if (i>0 & dim(train.good)[1]>0){all.passed=rbind(all.passed,train.good)}
if (i==0){
all.passed=train.good
}
i=i+1
}

print ("COMPLETED")
#

## falsely called over all clled:
## find the good threshold, then report posterior and fdr
all.passed[emp_p=="NaN", emp_p:=0]

# misnomer!

all.passed[,sinfo:=paste(sample,sampleType)]
all.passed[,fusionInfo:=junction]

min.reads=1

p.thresh=.1


## sqrt motivation by scaling of standard error for a junction cdf measuement as 1/sqrt(n)

## statistical cut-offs

test= all.passed[(( sinfo %like% "spork_0" ) & ((emp_p<p.thresh &  junction_cdf.y >.2 & numReads.y>min.reads) | ( numReads.y==1 & junction_cdf_lower.y >.5 & productPhat.y>.5))) | (emp_p<p.thresh & sinfo %like% "spork_1" &   numReads.y>min.reads & (fusionInfo %like% ":fusion-both"))]

sam=unique(test$sample)
test[,sampleid:=match(sample,sam)]

test[,sinfo:=paste(sample,sampleType)]
test[,fusionInfo:=junction]




write.table(file="Jan_2017_V2_supp_machet_processed.tab",test,sep="\t",quote=F)
#1:    chr9:TLE1:83956441:-:chr7:TMEM178B:141137408:+:no_fusion-acceptor_interchrom_inversion_invert,num=10,score=0.0,gap=0,don-dist:-242157,acc-dist:0,jct_ind=35803




