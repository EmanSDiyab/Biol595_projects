# NW.py
# Bioinformatics introduction, Fall 2020
# Eman Diyab

import sys

# You may define any helper functions for Needleman-Wunsch algorithm here

# Scoring values
gap_penalty = -1
match_award = 0
mismatch_penalty = -2

seq1 = ""
seq2 = ""

# Create an new matrix with zeros in every cell
def new_matrix(rows, cols):
    # Define an empty list
    matrix = []
    # Set up the rows of the matrix
    for x in range(rows):
        # For each row, add an empty list
        matrix.append([])
        # Set up the columns in each row
        for y in range(cols):
            # Add a zero to each column in each row
            matrix[-1].append(0)
    # Return the new matrix
    return matrix


# A function for determining the score between any two bases in alignment
def match_score(first, second):
    if first == second:
        return match_award
    elif first == '-' or second == '-':
        return gap_penalty
    else:
        return mismatch_penalty


# Do not change this function signature
def needleman_wunsch(seq1, seq2):
    # Store length of two sequences
    n = len(seq1)
    m = len(seq2)
    final_score = 0

    # Generate matrix of zeros to store scores
    score = new_matrix(m + 1, n + 1)

    # Calculate score table

    # Fill out first column
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i

    # Fill out first row
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j

    # Fill out all other values in the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate the score by checking the top, left, and diagonal cells
            match = score[i - 1][j - 1] + match_score(seq1[j - 1], seq2[i - 1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            # Record the maximum score from the three possible scores calculated above
            score[i][j] = max(match, delete, insert)
            # print(str(score[i][j]))

    # Traceback and compute the alignment

    # Create variables to store alignment
    alignment1 = ""
    alignment2 = ""

    # Start from the bottom right cell in matrix
    i = m
    j = n

    # We'll use i and j to keep track of where we are in the matrix, just like above
    while i > 0 and j > 0:  # end touching the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        # Check to figure out which cell the current score was calculated from,
        # then update i and j to correspond to that cell.
        if score_current == score_diagonal + match_score(seq1[j - 1], seq2[i - 1]):
            alignment1 += seq1[j - 1]
            alignment2 += seq2[i - 1]
            i -= 1
            j -= 1
            final_score += match_score(seq1[j - 1], seq2[i - 1])
        elif score_current == score_up + gap_penalty:
            alignment1 += seq1[j - 1]
            alignment2 += '-'
            j -= 1
            final_score += gap_penalty

        elif score_current == score_left + gap_penalty:
            alignment1 += '-'
            alignment2 += seq2[i - 1]
            i -= 1
            final_score += gap_penalty

    # Finish tracing up to the top left cell
    while j > 0:
        alignment1 += seq1[j - 1]
        alignment2 += '-'
        j -= 1
    while i > 0:
        alignment1 += '-'
        alignment2 += seq2[i - 1]
        i -= 1

    # Since we traversed the score matrix from the bottom right, our two sequences will be reversed.
    # These two lines reverse the order of the characters in each sequence.
    alignment1 = alignment1[::-1]
    alignment2 = alignment2[::-1]

    return (final_score,alignment1, alignment2)


if __name__=="__main__":
	#the sequence that need to be aligned
   	# seq1 = "AGCT-A"
   	# seq2 = "AGGTCA"

    arg_count = len(sys.argv)
    if arg_count == 1:
        print("Please write the file name after NW.py ")

    elif arg_count == 3:
        seq1 = sys.argv[1]
        seq2 = sys.argv[2]
        print("\n~~~~~~~~~~~~~~  Welcome :)  ~~~~~~~~~~~~~~")

        print("your input sequences are :" + str(seq1) + " and " + str(seq2))

    else:
        print("Please write the two sequences ONLY after NW.py ")


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# print(needleman_wunsch(seq1, seq2))
score, output1, output2 = needleman_wunsch(seq1, seq2)

#    '<alignment score>\n<alignment in seq1>\n<alignment in seq2>'
print(str(score) + "\n" + output1 + "\n" + output2)
print("~~~~~~~~~~~~~ THANKS BYE :) ~~~~~~~~~~~~~~")
