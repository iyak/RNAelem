RNAelem is a tool for learning sequence-structure motifs from a set of RNA sequences.

# Installation

```sh
git clone https://github.com/iyak/RNAelem.git
cd RNAelem
./waf configure --prefix=$HOME/local
./waf build test
./waf install
```
RNAelem is developed using Python3 and C++14, utilizing the `numpy` and `scipy` Python libraries. `waf` is a Python-based framework for configuring, compiling, and installing applications.

# Toy Example

First, organize the input data, which includes:

- A fasta file containing the positive sequence set
- A list of search patterns

Search patterns are candidates for local secondary structures of motifs, expressed in dot-bracket notation with `*` indicating insertion region. For example, the search pattern `((.*.))` denotes a local secondary structure with an arbitrary internal structure flanked by a stem of at least two bases and an opening loop.

To demonstrate, consider finding a motif from 76 seed sequences in Rfam's tRNA family (RF00005) with 'CAU' as their anticodon. These sequences exhibit the tRNA clover structure. Specify a single hairpin loop as the search pattern:
```sh
cat pattern
(.....)
```
Execute RNAelem to train a motif:
```sh
elem pipeline -p positive.fa -m pattern
```

|sequence logo|graph diagram|
----|----
|![tRNA profile](https://github.com/iyak/RNAelem/blob/master/material/prf-0.png)|![tRNA secondary structure](https://github.com/iyak/RNAelem/blob/master/material/rss-0.png)|

The results include sequence logos and graph diagrams displaying conserved sequential profiles and stable local secondary structures simultaneously. The sequence logo indicates base and base pair compositions, with the y-axis measuring information in bits. The trained results default to a directory named `elem_out/model-*` in the current working directory:

```sh
train.model   # contains model parameters
prf.[png/svg] # sequence logo
rss.[png/eps] # graph diagram
train.raw     # final alignment of the model to input data
train.interim # interim parameters
log           # log file
```

To align the best-trained model to new data:
```bash
elem scan -s sequences.fa -m elem_out/model-1/train.model
```
The alignment results are stored in `scan_out/scan.raw` and include detailed probabilistic and structural alignments:
```
id: @1                                # sequence id
start: [-38.5222,-20.6088, ...        # log probability that motif starts at the position
end: [-inf,-inf,-inf,-inf, ...        # log probability that motif ends at the position
inner: [-38.5222,-20.6088, ...        # log probability that the position is inside motif
psihat: [0,0,0,1,...                  # viterbi alignment of Profile CFG
motif region: 11 - 22                 # estimation of motif region
exist prob: 0.699748                  # probability that motif exists in the sequence
seq: AUAAUAUUUAGGUGCAACUCCUAAAUCCGCUA # sequence
rss: OOOOOLLLLLLLHHHHHHHRRRRRRROOOOOO # viterbi alignment of secondary structure
mot:           ((.......))            # viterbi alignment of motif
```

# Usage

To conduct de-novo motif discovery using RNAelem, multi-node parallelization on a cluster machine is required. This section describes the procedure for motif discovery on a single dataset, such as eCLIP data for an RBP. The training for each search pattern is distributed across nodes.

## Manual Parallelization

This section explains how to perform manual parallel computing without using Grid Engine. Start by preparing a FASTA file containing the positive sequence set and a text file with the search patterns. On the master node, initialize the training with the following command:

```sh
elem init --positive a.fa --pattern-list p.txt
```

Here, `a.fa` represents the positive sequence set, and `p.txt` contains the search patterns. By default, a directory named `elem_out` is created to store the training results. Subsequently, conduct the training for each search pattern on parallel nodes:

```sh
elem train --elem-out elem_out --pattern-index 1 # 1...N
```

Here, `N` is the total number of search patterns in the list. Each node trains using a specific index from `1` to `N`, and the outcomes are saved in `elem_out`. If using a distributed filesystem, consolidate the `elem_out` contents at a single location to proceed with model selection:

```sh
elem select --elem-out elem_out --num-motifs 3
```

An enrichment score is calculated for each motif. These scores are compared, and the motifs with the highest scores are selected. Further refine the chosen model using all the training data. If multiple motifs are selected, distribute their training across the nodes:

```sh
elem refine --elem-out elem_out --pattern-index 1 # 1...M
```

Here, `M` is the number of motifs previously selected. Each node conducts parallel training using an index from `1` to `M`. Consolidate the training outputs in `elem_out` if the filesystem is distributed. The optimal motif is stored under `elem_out/model-1`, and the suboptimal motifs as `elem_out/model-{2,...,M}`.

## Parallelization using Grid Engine

If the cluster machine supports job scheduling via Grid Engine, the `train` and `refine` steps can be replaced with the following commands:

```sh
elem train --elem-out elem_out --array
elem refine --elem-out elem_out --array
```

Alternatively, replace the entire procedure from `init` to `refine` with the following single command:

```sh
elem pipeline --positive a.fa --pattern-list p.txt --array
```

It is essential to configure the job execution commands based on the variant and version of the Grid Engine. At the beginning of `script/elem`, edit a Python dictionary to set these configurations:

```python
ge={
    "cmd":"qsub",                               # command to submit a job
    "array":"-t 1-:N:",                         # :N: represents number of array-jobs
    "cwd":"-cwd",                               # option to set cwd as working directory
    "envvar":"-V",                              # option to pass environment variables
    "stdout":"-o :file:",                       # :file: is the file to write stdout
    "stderr":"-e :file:",                       # :file: is the file to write stderr
    "shell":"-S $SHELL",                        # option to set current shell as interpreter
    "memory":"-l mem_req=:mem:G,s_vmem=:mem:G", # :mem: is memory to use
    "cpu":"-pe def_slot :cpu:",                 # :cpu: is num of cpu to use
    "task_id":"SGE_TASK_ID",                    # env var to store task ID
    "sync":"-sync y",                           # option to wait until job terminates
    "other":"",                                 # any other options you want to specify
    }
```