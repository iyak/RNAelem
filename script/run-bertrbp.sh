#!/bin/bash
# usage: run-bertrbp.sh <positive.fa> <negative.fa> <test.fa> <outdir>


set -eu
function realpath { echo $(cd $(dirname $1); pwd)/$(basename $1); }


P_fa=`realpath $1`
N_fa=`realpath $2`
T_fa=`realpath $3`
OUTDIR=`realpath $4`
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


mkdir -p $OUTDIR/{training_sample_finetune,test_sample_finetune,finetuned_model}
cd $OUTDIR


TIMEFORMAT='ELEM_TIME %3R';
time (\
    python $CWD/_run-bertrbp-h0.py --positive $P_fa --negative $N_fa > $OUTDIR/training_sample_finetune/train.tsv ;\
    python $CWD/_run-bertrbp-h0.py --no-info  $T_fa                  > $OUTDIR/test_sample_finetune/dev.tsv ;\
    docker run --gpus all --volume $OUTDIR:/root/github/kkyamada/bert-rbp/sample_dataset/DATA bert-rbp:dev ;\
) \
2> >(tee stderr.txt >&2) | tee stdout.txt


python $CWD/_run-bertrbp-h1.py --test $T_fa --pred $OUTDIR/finetuned_model/pred_results.npy > $OUTDIR/res0 ;\
