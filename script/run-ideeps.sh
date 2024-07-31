#!/bin/bash
# usage: run-ideeps.sh <positive.fa> <negative.fa> <test.fa> <outdir>


set -eu
function realpath { echo $(cd $(dirname $1); pwd)/$(basename $1); }


P_fa=`realpath $1`
N_fa=`realpath $2`
T_fa=`realpath $3`
OUTDIR=`realpath $4`
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


mkdir -p $OUTDIR
cd $OUTDIR


bash $DIR/_run-ideeps.h0.sh\
    $P_fa\
    $N_fa\
    $T_fa\
    $OUTDIR/train.fa.gz\
    $OUTDIR/test.fa.gz
rm -f $OUTDIR/structure.gz


docker run\
    --volume $OUTDIR:/out\
    --env TIMEFORMAT='ELEM_TIME %3R'\
    --entrypoint /bin/bash\
    --rm\
    ideeps:dev -c "\
    time (\
        python -W ignore ideeps.py\
            --train=True\
            --data_file=/out/train.fa.gz\
            --model_dir=/out/model;\
        mv -f /out/structure.gz /out/train.structure.gz;\
    );\
    time (\
        python -W ignore ideeps.py\
            --predict=True\
            --data_file=/out/test.fa.gz\
            --model_dir=/out/model\
            --out_file=/out/out.txt;\
        mv -f /out/structure.gz /out/test.structure.gz;\
    );\
    "\
2> >(tee stderr.txt >&2) | tee stdout.txt
