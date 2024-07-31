#!/bin/bash
# usage: qsub -t 1-367 _run-rnaelem-h1.sh <train.fq> <test.fq> <motifs> <outdir> <dataset ID>

#$ -cwd
#$ -e /home/iyak/.ugeerr
#$ -o /home/iyak/.ugeout
#$ -S /bin/bash
##$ -m es -M h.m@tuta.io
#$ -V
#$ -l s_vmem=3G,mem_req=3G
#$ -r no
#$ -pe def_slot 4

function realpath { echo $(cd $(dirname $1); pwd)/$(basename $1); }
trainfq=`realpath $1`
testfq=`realpath $2`
motifs=`realpath $3`
outdir=`realpath $4`
dID=$5
bench=$HOME/proj/RNAelem/bench/18-perf-compe
date=171102

set -eu

cd $outdir
    mkdir -p pattern-$SGE_TASK_ID
    cd pattern-$SGE_TASK_ID
        motif=`sed "${SGE_TASK_ID}q;d" $motifs`
        #$HOME/proj/RNAelem-dev/build/bin/RNAelem \
        $bench/$date/RNAelem \
            -f ../train.fq \
            -m $motif \
            -w 50 \
            -c 30 \
            --out1 a.model \
            --out2 a.svg \
            --out3 train.raw \
            -t 4 \
            -i 300 \
            --font /usr/share/fonts/vlgothic/VL-Gothic-Regular.ttf \
            2>log
            #--theta-softmax \
            #--lik-ratio \
            #--param-set 4-43 \
            #--rho-s 0 \
            #--rho-s 1e-3 \
        #$HOME/proj/RNAelem-dev/build/bin/RNAelem scan\
        $bench/$date/RNAelem scan \
            -f ../test.fq \
            -q a.model \
            --out1 raw \
            -t 4 \
            2>>log
        python $bench/$date/main2.py $dID $SGE_TASK_ID 2>>log
