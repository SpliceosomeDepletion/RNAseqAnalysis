#!/bin/bash

# Iterate over the the raw files and run the pipeline
# Create a directory for each run
RUNSCRIPT=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Internship/PRPF8_rnaseq2019/src/RNAseqAnalysis/rnaseq_pipeline_fromalignment.sh;
PARAMFILE=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Internship/PRPF8_rnaseq2019/src/RNAseqAnalysis/rnaseq_pipeline_config.sh;

DATADIR=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Internship/PRPF8_rnaseq2019/data/RNAseqRaw/; # Use same raw data as in previous project
CURRDIR=$(pwd -LP);
OUTDIR=/nfs/nas21.ethz.ch/nas/fs2102/biol_ibt_usr_s1/mfrank/Internship/PRPF8_rnaseq2019/data/RNAseq/Alignments/;
cd $DATADIR;
for i in *r_1.fq;
do
  R1=$i;
  R2=$(echo $R1 | sed 's/r_1\.fq/r_2\.fq/');
  NAME=$(echo $R1 | sed 's/\.811\.s_1.*//');
  echo $NAME
  #echo $R2;
#  mkdir ${OUTDIR}/$NAME
  $RUNSCRIPT 1 3 $NAME ${OUTDIR}/$NAME $PARAMFILE ${DATADIR}/$R1 ${DATADIR}/$R2;
  #cd ..;
done

cd $CURRDIR;
