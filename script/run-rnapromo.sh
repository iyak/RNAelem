#!/bin/bash
# usage: run-rnapromo.sh <positive.fa> <negative.fa> <test.fa> <outdir>


function realpath { echo $(cd $(dirname $1); pwd)/$(basename $1); }
P_fa=`realpath $1`
N_fa=`realpath $2`
T_fa=`realpath $3`
OUTDIR=`realpath $4`
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


set -eu
mkdir -p $OUTDIR
cd $OUTDIR


docker run \
    --volume $P_fa:/p.fa \
    --volume $N_fa:/n.fa \
    --volume $T_fa:/t.fa \
    --volume $OUTDIR:/out \
    --env TIMEFORMAT='ELEM_TIME %3R' \
    --entrypoint /bin/bash \
    rnapromo:dev -c " \
    time ( \
        perl rnamotifs08_motif_finder.pl \
            -n 1 -min 3 -max 11 \
            -positive_seq /p.fa \
            -negative_seq /n.fa; \
    ); \
    time ( \
        cat /t.fa \
        | bin/ViennaRNA/ViennaRNA-1.6/Progs/RNAfold -noPS \
        | awk '0==NR%3{\$0=\$1}1' \
        | paste - - - \
        > /out/scan.in; \
        perl rnamotifs08_motif_match.pl \
            /out/scan.in \
            -cm cm_1.tab -b \
            > /out/cm_1.pos; \
    )" \
2> >(tee stderr.txt >&2) | tee stdout.txt


$DIR/_run-rnapromo-h0.py cm_1.pos $T_fa >$OUTDIR/res0 2>>$OUTDIR/stderr.txt 
