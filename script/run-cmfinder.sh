#!/bin/bash
# usage: run-cmfinder.sh <positive.fa> <negative.fa> <test.fa> <outdir>
# note: it does not use negative.fa


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
    --volume $T_fa:/t.fa \
    --volume $OUTDIR:/out \
    --env TIMEFORMAT='ELEM_TIME %3R' \
    --entrypoint /bin/bash \
    cmfinder:dev -c " \
    awk 'NR%2==1{print \">\"int(NR/2)}NR%2==0{print \$1}' /p.fa > /out/p.fa
    set -x;
    time (\
        perl bin/cmfinder.pl\
            -v\
            -n 3\
            -m 10\
            -M 20\
            -s 1\
            /out/p.fa;\
    );\
    time (\
        infernal/cmsearch\
            /out/p.fa.cm.h1.1\
            /t.fa\
            >/out/t.fa.match;
    );\
    "\
2> >(tee stderr.txt >&2) | tee stdout.txt


python $DIR/_run-cmfinder-h0.py $OUTDIR/t.fa.match $OUTDIR/t.fa >res0 2>>log

