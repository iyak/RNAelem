#! /usr/bin/env python
# usage: python _run-graphprot2-h0.py <pred_out> > res0


import sys
from pathlib import Path


out      = Path(sys.argv[1]) / "profiles.winext20.out"
bed      = Path(sys.argv[1]) / "profiles.winext20.peak_pos.ext.bed"
seq_wise = Path(sys.argv[1]) / "whole_site_scores.out"


with open(out) as out_fi, open(bed) as bed_fi, open(seq_wise) as sw_fi:
    # Seq-wise score.
    seq_wise_scores = {}
    lines = sw_fi.readlines()
    for line in lines:
        header, score = line.strip().split("\t")
        seq_wise_scores[header] = float(score)

    # Peak region.
    peaks = {}
    lines = bed_fi.readlines()
    for line in lines:
        header, beg, end, name, score, strand = line.strip().split("\t")
        if (header not in peaks) or (peaks[header][2] < score):
            peaks[header] = [int(beg), int(end), float(score)]

    # Position-wise score.
    lines = out_fi.readlines()
    for line in lines:
        line = line.strip()
        if ">" == line[0]:
            header = line[1:]
            continue
        pos, base, pos_score = line.split()
        pos, pos_score = int(pos), float(pos_score)
        seq_wise_score = seq_wise_scores.get(header, -99999)

        # Print position-wise scores.
        sys.stdout.write("\t".join(map(str,[
            header,
            pos - 1, # GraphProt2 outputs 1-origin.
            1 if header in peaks and peaks[header][0] <= pos-1 < peaks[header][1] else 0,
            seq_wise_score,
            pos_score,
            ]))+"\n")
