# NW.py
# Project1, Bioinformatics introduction, Fall 2020
# Eman Diyab

import sys
# Import pairwise2 module
from Bio import pairwise2
from Bio import SeqIO
from Bio import Align
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
import itertools


X = ""
Y = ""
i = j = 1
match = mismatch = gap_open = gap_ext = 0
seq1 = seq2 = ""


def file_read(input_file):
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequence = sequence.replace("\n", "")
        sequence = sequence.replace("\r", "")
        return sequence


def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


def read_seq(inputfile):
    with open(inputfile, "r") as f:
        seq = f.read()
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq


def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))


def alignment_print(sequence, alignment):
    seq = file_read(sequence)
    Ahmed = alignment.split("\n")
    chunk_count = 0

    seq_lines = (i.strip() for i in seq.splitlines())
    align_lines0 = (al1.strip() for al1 in Ahmed[0].splitlines())
    align_lines1 = (al2.strip() for al2 in Ahmed[1].splitlines())
    align_lines2 = (al3.strip() for al3 in Ahmed[2].splitlines())

    for (line, line0, line1, line2) in itertools.zip_longest(seq_lines, align_lines0, align_lines1, align_lines2):
        for (chunk, chunk0, chunk1, chunk2) in itertools.zip_longest((chunkstring(line, 90)), (chunkstring(line0, 30)), (chunkstring(line1, 30)), (chunkstring(line2, 30))):
            chunk_count += 1
            print(chunk)
            print(chunk0)
            print(chunk1)
            print(chunk2)
            print("")


# local alignments match score =1, no penalties for gaps or mismatch
def local_alignment(x, y):
    alignments = pairwise2.align.localxx(x, y)
    print("All possible local alignments:")
    j = 0
    for a in alignments:
        print("Alignment :" + str(j) + ":\n")
        print("Score: is " + str(a.score))
        print("~~~~~~~~~~~~~~~~~~~")
        alignment_print(seq1, format_alignment(*a))
        j += 1
        print("~~~~~~~~~~~~~~~~~~~")
    anything_else()


# global alignments match score =1, no penalties for gaps or mismatch
def global_alignment(x, y):
    alignments = pairwise2.align.globalxx(x, y)
    print("All possible global alignments:")
    i = 0
    for a in alignments:
        print("Alignment " + str(i) + ":\n")
        print("Score: is " + str(a.score))
        print("~~~~~~~~~~~~~~~~~~~")
        alignment_print(seq1, format_alignment(*a))
        i += 1
        print("~~~~~~~~~~~~~~~~~~~")
    anything_else()


# this function is for global alignment with specific scores
# inputs are (seq1,seq2,matching award,mismatching penalty,opening a gap,gap extension)
def global_with_different_scores(x, y, match, mismatch, gap_open, gap_ext):
    alignments = pairwise2.align.globalms(x, y, match, mismatch, gap_open, gap_ext)
    # Use format_alignment method to format and print the alignments in the list
    print("All possible global alignments:")
    i = 0
    for a in alignments:
        print("Alignment " + str(i) + ":\n")
        print("Score: is " + str(a.score))
        print("~~~~~~~~~~~~~~~~~~~")
        alignment_print(seq1, format_alignment(*a))
        i += 1
        print("~~~~~~~~~~~~~~~~~~~")


# this function is for local alignment with specific scores
# inputs are (seq1,seq2,matching award,mismatching penalty,opening a gap,gap extension)
def local_with_different_scores(x, y, match1, mismatch1, gap_open1, gap_ext1):
    alignments = pairwise2.align.localms(x, y, match1, mismatch1, gap_open1, gap_ext1)

    print("All possible global alignments:")
    i = 0
    for a in alignments:
        print("Alignment " + str(i) + ":\n")
        print("Score: is " + str(a.score))
        print("~~~~~~~~~~~~~~~~~~~")
        alignment_print(seq1, format_alignment(*a))
        i += 1
        print("~~~~~~~~~~~~~~~~~~~")


def local_handler_independent():
    match1 = float(input("Please enter the match value: "))
    print("Great\n")
    mismatch1 = float(input("Please enter the mismatch value: "))
    print("Great\n")
    gap_open1 = float(input("Please enter the gap open value: "))
    print("Great\n")
    local_with_different_scores(X, Y, match1, mismatch1, gap_open1, gap_open1)
    anything_else()


def local_handler_affine():
    match1 = float(input("Please enter the match value: "))
    print("Great\n")
    mismatch1 = float(input("Please enter the mismatch value: "))
    print("Great\n")
    gap_open1 = float(input("Please enter the gap open value: "))
    print("Great\n")
    gap_ext1 = float(input("Please enter the gap extension value: "))
    print("Great\n")
    local_with_different_scores(X, Y, match1, mismatch1, gap_open1, gap_ext1)
    anything_else()


def global_handler_independent():
    match1 = float(input("Please enter the match value: "))
    mismatch1 = float(input("Please enter the mismatch value: "))
    gap_open1 = float(input("Please enter the gap penalty value: "))
    global_with_different_scores(X, Y, match1, mismatch1, gap_open1, gap_open1)
    anything_else()


def global_handler_affine():
    match1 = float(input("Please enter the match value: "))
    mismatch1 = float(input("Please enter the mismatch value: "))
    gap_open1 = float(input("Please enter the gap open value: "))
    gap_ext1 = float(input("Please enter the gap extension value: "))
    global_with_different_scores(X, Y, match1, mismatch1, gap_open1, gap_ext1)
    anything_else()


def anything_else():
    yes = {'yes', 'y', 'ye', ''}
    no = {'no', 'n'}

    choice = input("\n\nDo you want another alignment with the same files??\t Y or N\n").lower()
    if choice in yes:
        print_menu()
    elif choice in no:
        call_function(6)
    else:
        sys.stdout.write("Please respond with 'yes' or 'no'")
        anything_else()


def print_menu():
    choice = ""
    print("Please choose from the following menu:")
    print("......................................\n")
    print("For local alignments with the default settings (match=1, no gap or mismatch penalties), press 0")
    print("")
    print("For local alignments with independent gaps, press 1")
    print("")
    print("For local alignments with affine gaps, press 2")
    print("")
    print("For global alignments with the default settings (match=1, no gap or mismatch penalties), press 3")
    print("")
    print("For global alignments with independent gaps, press 4")
    print("")
    print("For global alignments with affine gaps, press 5")
    print("\nTo exit press 6\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    choice = int(input("  \n"))
    call_function(choice)


def call_function(i1):
    if i1 == 0:
        local_alignment(X, Y)

    elif i1 == 1:
        local_handler_independent()

    elif i1 == 2:
        local_handler_affine()

    elif i1 == 3:
        global_alignment(X, Y)

    elif i1 == 4:
        global_handler_independent()

    elif i1 == 5:
        global_handler_affine()

    else:
        print("~~~~~~~~~~~~~~~~~~~~~~~~~ Thanks, BYE ~~~~~~~~~~~~~~~~~~~~~~~~~")
        exit()


if __name__=="__main__":

    arg_count = len(sys.argv)
    if arg_count != 3:
        print("Please follow this pattern ( python align_Ediyab.py $DNA_sequence$ $Protein_sequence$ fasta files )")

    else:
        seq1 = sys.argv[1]
        seq2 = sys.argv[2]
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~ Welcome ~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("your DNA sequence is: " + str(seq1) + " & Protein sequence is: " + str(seq2))
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        dna = file_read(seq1)
        print("Input DNA sequence is :\n")
        print(dna)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        X = translate(dna)
        print("Translated DNA sequence is :\n")
        print(X)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        Y = file_read(seq2)
        print("Input protein sequence is :\n")
        print(Y)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        # X = translate(seq1)
        # Y = seq2.upper()
        print_menu()




# TCTTTTGGGCGCTAACCCATTGATCATCTTGACTCTTCGACGTTCAACAGCGAAGTAACG
# SFGR_PIDHLDSSTFNSEVT