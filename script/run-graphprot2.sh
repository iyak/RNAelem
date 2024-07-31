#!/bin/bash
# usage: run-graphprot2.sh <positive.fa> <negative.fa> <test.fa> <outdir>


set -eu
function realpath { echo $(cd $(dirname $1); pwd)/$(basename $1); }


P_fa=`realpath $1`
N_fa=`realpath $2`
T_fa=`realpath $3`
OUTDIR=`realpath $4`
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


mkdir -p $OUTDIR
cd $OUTDIR
set -x


docker run --gpus all \
    --rm \
    --volume ${P_fa}:/p.fa \
    --volume ${N_fa}:/n.fa \
    --volume ${T_fa}:/t.fa \
    --volume ${OUTDIR}:/out \
    --env TIMEFORMAT='ELEM_TIME %3R' \
    --entrypoint /bin/bash \
    graphprot2:dev -c " \
    time ( \
        graphprot2 gt --in /p.fa --neg-in /n.fa --out /out/gt_out; \
        graphprot2 gt --in /p.fa --neg-in /n.fa --out /out/gt_out; \
        graphprot2 train --in /out/gt_out --out /out/train_out; \
    ); \
    time ( \
        graphprot2 gp --in /t.fa --train-in /out/train_out --out /out/gp_out; \
        graphprot2 predict --in /out/gp_out --train-in /out/train_out --out /out/pred_out --mode 1; \
        graphprot2 predict --in /out/gp_out --train-in /out/train_out --out /out/pred_out --thr -1 --mode 2; \
    ) " \
2> >(tee stderr.txt >&2) | tee stdout.txt


python $DIR/_run-graphprot2-h0.py pred_out >res0
