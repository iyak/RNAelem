![RNAelem logo for Github](https://github.com/iyak/RNAelem/blob/master/material/github_logo.png)

RNAelem is a tool for learning secondary structural motifs from a set of RNA sequences.
A secondary structure motif is a pair of local secondary structures and partial sequences that are conserved in RNA sequences.
Rfold model [1] assesses the fee energy stability of local secondary structures, and Profile CFG evaluates the conservation of partial sequences.
The motifs are extracted and reported in the sequence logo and graph diagram.

## Installation

```bash
$ git clone htpps://github.com/iyak/RNAelem.git
$ cd RNAelem
$ ./waf configure --prefix=$HOME/local # for example
$ ./waf build test
$ ./waf install
```
RNAelem is written using Python3 and C++ 14.
Python library `numpy` is used.

## Toy example

First, you need to organize the input data.
The following files are required for learning.

- fasta file that contains positive set
- List of search patterns

A list of search patterns is candidates for local secondary structures of the motif.
It is written in the following manner:
```
pattern 1
pattern 2
...
```
Search pattern is in dot bracket notation of the secondary structure, added by `*` which indicates any structure.
For example, search pattern `((.*.))` represents a local secondary structure consisting of an arbitrary structure which is sandwiched between a stem of length 2 or greater and an opening loop of length 1 or greater.
Specifically, internal loop and hairpin loop are applicable.

You can optionally pass negative set sequences in fasta format.

##
Let's find a motif from 76 seed sequences in the Rfam's tRNA family (RF00005) that have 'CAU' as their anticodon.
These sequences contain the tRNA clover structure.
Specify a single hairpin loop as the search pattern.
```bash
$ cat pattern
(.....)
```
Run RNAelem to train a motif.
```bash
$ elem -p positive.fa -m patterns
```
As a result, you will obtain the following motifs.

|sequence logo|graph diagram|
----|---- 
|![tRNA profile](https://github.com/iyak/RNAelem/blob/master/material/prf-0.png)|![tRNA secondary structure](https://github.com/iyak/RNAelem/blob/master/material/rss-0.png)|

The sequence logo shows both the conserved sequential profile and the stable local secondary structure at the same time.
In addition to the normal base composition, the base pair composition is also displayed.
For base pairs, the same compositions are displayed, and the partner's column is grayed-out.
The y axis indicates the information in bit, which ranges from 0 to 2 for base composition, and 0 to log 6 for base pair.
The column index shows the base pairs and the average number of repeats for each secondary structural context.
Interpreting the results based on the above, the anticodon arm has been extracted that contains an anticodon 'CAU' and the conserved 'U' and 'A' on both sides.

The trained results are stored in a folder named `elem_out/summary` by default in the current working directory.
```bash
$ ls elem_out/summary
model-*.txt           # contains model parameters
prf-*.[png/svg]       # sequencial logo
rss-*.[png/eps]       # graph diagram
raw-*.txt             # final alignment of the model to input data
model-*-interim.txt   # interim parameters
log-*.txt             # log file
```
`*` is replaced by a number between 0 and 2 by default, and a lower number corresponds to a search pattern with a higher fitness.

You can also align the trained model to the new data.
For example, if you want to align the best model to a sequence set `sequences.fa`, you can do the following:
```bash
$ elem scan -s sequences.fa -m elem_out/summary/model-0.txt
```

The aligned results are stored in `scan_out/scan.raw` by default in the current working directory.
This file is a sequence of the following blocks:
```
id: @1                                 # sequence id
start: [-38.5222,-20.6088, ...         # log probability that motif starts at the position
end: [-inf,-inf,-inf,-inf, ...         # log probability that motif ends at the position
inner: [-38.5222,-20.6088, ...         # log probability that the position is inside motif 
psihat: [0,0,0,1,...                   # viterbi alignment of Profile CFG
motif region: 11 - 22                  # estimation of motif region
exist prob: 0.699748                   # probability that motif exists in the sequence
seq: AUAAUAUUUAGGUGCAACUCCUAAAUCCGCUA  # sequence
rss: OOOOOLLLLLLLHHHHHHHRRRRRRROOOOOO  # viterbi alignment of secondary structure
mot:           ((.......))             # viterbi alignment of motif
```

## Settings for grid engine

In order to fully operate RNAelem and explore de novo motifs, parallelization is basically required in the cluster machine.
How you submit a job for a cluster machine depends on your environment.
You must configure it by directly modifying the `elem`, the executable file, so that you can submit jobs from RNAelem.
Line 16 of elem is as follows:
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
As you can see, the python dictionary is used to manage the configuration.
For options that require arguments, use variables such as `:N:` to specify the location of the arguments.
Set them and start the tool with `-a` option to make sure the job is appropriately queued.

## Reference
1. Kiryu, H., Kin, T., & Asai, K. (2007). Rfold: an exact algorithm for computing local base pairing probabilities. Bioinformatics, 24(3), 367-373.
