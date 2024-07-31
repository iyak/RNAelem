#!/bin/bash
# usage: run-graphprot.sh <positive.fa> <negative.fa> <test.fa> <outdir>


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
    --gpus all \
    --volume $P_fa:/p.fa \
    --volume $N_fa:/n.fa \
    --volume $T_fa:/t.fa \
    --volume $OUTDIR:/out \
    --env TIMEFORMAT='ELEM_TIME %3R' \
    --entrypoint /bin/bash \
    --rm \
    brianyee/graphprot:1.1.7 -c "\
    cd /out; \
    set -x; \
    time ( \
        perl /opt/GraphProt-1.1.7/GraphProt.pl \
            --action train --fasta /p.fa --negfasta /n.fa ;\
    ); \
    time ( \
        perl /opt/GraphProt-1.1.7/GraphProt.pl \
           --action predict --fasta /t.fa --model GraphProt.model ;\
        perl /opt/GraphProt-1.1.7/GraphProt.pl \
           --action predict_profile --fasta /t.fa --model GraphProt.model ;\
        perl /opt/GraphProt-1.1.7/GraphProt.pl \
           --action predict_has --fasta /t.fa --model GraphProt.model ;\
    )" \
2> >(tee stderr.txt >&2) | tee stdout.txt


python $DIR/_run-graphprot-h0.py $OUTDIR/GraphProt.profile $T_fa >$OUTDIR/res0 2>>$OUTDIR/stderr.txt
