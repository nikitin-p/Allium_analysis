#!/usr/bin/python3
from Bio import AlignIO
import csv
import sys
in_file_name=sys.argv[1]
align = AlignIO.read(in_file_name, "fasta")
out_fh = open("window_count.csv",'w')

def count_window(aln):
    out_dict={nuc:0 for nuc in "ACGT-"}
    for rec in aln:
        uc_seq=rec.seq.upper()
        for nuc in "ACGT-":
            out_dict[nuc]+=uc_seq.count(nuc) 
    return(out_dict)

d_w=csv.DictWriter(out_fh,delimiter=';',fieldnames=['Pos']+[s for s in "ACGT-"])
d_w.writeheader()

step=100
for i in range(step,align.get_alignment_length(),step):
    out_d=count_window(align[:,i-step:i])
    out_d['Pos']=i
    d_w.writerow(out_d)
