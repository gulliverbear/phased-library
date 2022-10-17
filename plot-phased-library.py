"""
Plot base distributions of sequencing library with phased primers
Agilent library ends are: CCACCTGGT ... TCGATACTC
then add on rBC: CCGNNNNNNNNCANNNNNNNNAGGACGACTCTATCAGTCGG
R1 will have an A and then the phased bases
R2 with just go directly to the phased bases
"""

import collections
from dataclasses import dataclass
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import random

def reverse_complement(s):
    rc_dict = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    return "".join([rc_dict[i.upper()] for i in s[::-1]])

def replace_n_with_random_base(s):
    new_s = ""
    for base in s:
        if base.upper() == "N":
            new_s += random.choice("ACGT")
        else:
            new_s += base
    return new_s

def add_reads(illumina_reads, library, r1_phase, r2_phase, rc_flag):
	"""
	add each read from the library with given flanking phased adapter sequences
	"""
    for s in library:
        r = "A" + r1_phase + s + rbc + r2_phase
        r = replace_n_with_random_base(r)
        if rc_flag:
            r = reverse_complement(r)
        illumina_reads.append(r)

def plot_reads(reads, read_length, rc_flag):
    pass
    # To Do
        
@dataclass(eq=True, frozen=True) # set eq and frozen to True so it will be hashable
class Primer():
    read_number: str
    phased_length: str
    orientation: str
    frozen=True

rbc = "CCGNNNNNNNNCANNNNNNNNAGGACGACTCTATCAGTCGG"

# read in library
library_file = pathlib.Path(".")

with open(library_file) as f:
    library = f.read().splitlines()

s = """s2.PCR2.R1.0bp.REV	CTTTCCCTACACGACGCTCTTCCGATCTA CCGACTGATAGAGTCGTC
s2.PCR2.R1.1bp.REV	CTTTCCCTACACGACGCTCTTCCGATCTA A CCGACTGATAGAGTCGTC
s2.PCR2.R1.2bp.REV	CTTTCCCTACACGACGCTCTTCCGATCTA TA CCGACTGATAGAGTCGTC
s2.PCR2.R1.3bp.REV	CTTTCCCTACACGACGCTCTTCCGATCTA ATT CCGACTGATAGAGTCGTC
s2.PCR2.R2.0bp.FOR	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT CCACCTGGTCGAGATATC
s2.PCR2.R2.1bp.FOR	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT A CCACCTGGTCGAGATATC
s2.PCR2.R2.2bp.FOR	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT TA CCACCTGGTCGAGATATC
s2.PCR2.R2.3bp.FOR	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT ATT CCACCTGGTCGAGATATC
s2.PCR2.R2.0bp.REV	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT CCGACTGATAGAGTCGTC
s2.PCR2.R2.1bp.REV	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT A CCGACTGATAGAGTCGTC
s2.PCR2.R2.2bp.REV	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT TA CCGACTGATAGAGTCGTC
s2.PCR2.R2.3bp.REV	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT ATT CCGACTGATAGAGTCGTC
s2.PCR2.R1.0bp.FOR	CTTTCCCTACACGACGCTCTTCCGATCTA CCACCTGGTCGAGATATC
s2.PCR2.R1.1bp.FOR	CTTTCCCTACACGACGCTCTTCCGATCTA A CCACCTGGTCGAGATATC
s2.PCR2.R1.2bp.FOR	CTTTCCCTACACGACGCTCTTCCGATCTA TA CCACCTGGTCGAGATATC
s2.PCR2.R1.3bp.FOR	CTTTCCCTACACGACGCTCTTCCGATCTA ATT CCACCTGGTCGAGATATC"""

p = """R1.0bp.FOR	R2.3bp.REV
R1.1bp.FOR	R2.2bp.REV
R1.2bp.FOR	R2.1bp.REV
R1.3bp.FOR	R2.0bp.REV
R2.0bp.FOR	R1.3bp.REV
R2.1bp.FOR	R1.2bp.REV
R2.2bp.FOR	R1.1bp.REV
R2.3bp.FOR	R1.0bp.REV"""

# generate dictionary to get phased adapter sequence from a primer dataclass key
d = {}
for line in s.split("\n"):
    name, seq = line.split("\t")
    
    read_number, phased_length, orientation = name.split(".")[-3:]
    d[Primer(read_number, phased_length, orientation,)] = seq

illumina_reads = []
for l in p.split("\n"):
    p1, p2 = l.split("\t")
    p1 = Primer(*p1.split("."))
    p2 = Primer(*p2.split("."))
    
    s1 = d[p1]
    s2 = d[p2]
    
    if len(s1.split()) == 2:
        phase1 = ""
    else:
        phase1 = s1.split()[1]
    if len(s2.split()) == 2:
        phase2 = ""
    else:
        phase2 = s2.split()[1]
    
    rc_flag = (p1.read_number, p1.orientation) in (("R1","FOR"), ("R2","REV"))
    
    if p1.read_number == "R1":
        add_reads(illumina_reads, library, phase1, phase2, rc_flag)
    else:
        add_reads(illumina_reads, library, phase2, phase1, rc_flag)
        
plot_reads(illumina_reads, 250, rc_flag=False)
plot_reads(illumina_reads, 250, rc_flag=True)