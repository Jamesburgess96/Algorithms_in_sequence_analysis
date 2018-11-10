#!/usr/bin/python3

"""Template for dynamic programming assignment.

The code in the template is compatible with both Python 2 and Python 3
When you finish this code, it should be at least compatible with Python 3.
"""

# Packages for commandline options:
import argparse
import sys
import pickle as pk


# Built-in exchange matrices.
with open('substitution_matrices/identity.pkl', 'rb') as f:
    identity = pk.load(f)

with open('substitution_matrices/pam250.pkl', 'rb') as f:
    pam250 = pk.load(f)

with open('substitution_matrices/blosum62.pkl', 'rb') as f:
    blosum62 = pk.load(f)


def get_args():
    """Collect the inputs."""
    parser = argparse.ArgumentParser(
        prog='PROG',
        usage='%(prog)s [options]',
        description='Aligning two sequences',
        epilog='The code was co-opted from Anton Feenstra\'s and'
        'modified by Cico Zhang'
    )
    parser.add_argument('-f', '--fasta', dest='fasta', metavar='FILE',
                        required=True, help='input alignment file (fasta)')
    parser.add_argument('-e,', '--exchange_matrix', dest='exchange_matrix',
                        metavar='SUBSTITUTION MATRIX NAME', help='Substitution '
                        'matrix: pam250, blosum62 or identity',
                        default='pam250')
    parser.add_argument('-l', '--local', dest='align_local',
                        action='store_true', help='Local alignment',
                        default=False)
    parser.add_argument('-g', '--global', dest='align_global',
                        action='store_true', help='Global alignment',
                        default=False)
    parser.add_argument('-s', '--semi_global', dest='align_semiglobal',
                        action='store_true', help='Semi-global alignment',
                        default=False)
    parser.add_argument('-p', '--penalty', dest='gap_penalty', type=int,
                        help='Gap penalty', default=2)
    parser.add_argument('-o', '--output', dest='alignment', metavar='FILE',
                        default='output.align', help='The file to store the alignment')
    parser.add_argument('-m', '--score_matrix', dest='score_matrix',
                        metavar='FILE', default='output.align',
                        help='The file to store the score matrix')
    parser.add_argument('-v', dest='print_on_screen', action='store_true',
                        help='Print the output (alignment(s) and score '
                        'matrix) on the screen', default=False)

    args = parser.parse_args()

    if args.fasta is None:
        sys.exit('Error: no input file (fasta)')

    if not (args.align_local or args.align_global or args.align_semiglobal):
        sys.exit('Error: No alignment strategy is given: global, local or '
                 'semi-global')
    if args.align_local + args.align_global + args.align_semiglobal > 1:
        sys.exit('Error: More than one alignment strategy is given.')

    if args.exchange_matrix not in ['pam250', 'blosum62', 'identity']:
        sys.exit('Unknown exchange matrix ' + args.exchange_matrix)

    return args


class Sequence:
    """Stores a sequence object."""

    def __init__(self, Label="", Sequence=""):
        """Initialize a new Sequence object.

        Label -- identifier of sequence (text)
        Sequence -- sequence string in single-letter alphabet
        """
        self.Label = Label
        self.Sequence = Sequence

    # this makes that you can do 'print sequence' and get nice output:
    def __str__(self):
        """Return string representation of a Sequence object."""
        # newline-delimited values of all the attributes
        return ">%s\n%s" % (self.Label, self.Sequence)


def readSequences(lines):
    """Return Sequences object.

    lines -- list of lines or any object that behaves like it

    This routine parses a fasta file and returns a list of Sequence objects
    containing the sequences with label and sequence data set
    """
    seqs = []
    label = None
    seq_lines = []
    for line in lines:
        line = line.strip()      # strip off white space
        if not line:             # skip empty lines
            continue
        if line.startswith(';'):  # ignore comment lines
            continue
        # check for start of next sequence:
        if line.startswith('>'):  # label line
            # first, store the previous sequence if we had one:
            if seq_lines:
                seqs.append(Sequence(label, ''.join(seq_lines)))
                seq_lines = []
            # get the label (name) for the next sequence
            label = line[1:].strip()
        else:
            # collect all lines with sequence information for this sequence:
            seq_lines.append(line)
    # take care of the last sequence in the file
    seqs.append(Sequence(label, ''.join(seq_lines)))
    return seqs

def do_global_alignment(sequences, matrix, penalty):
    """Do pairwise global alignment using DP."""
    #########################
    # INSERT YOUR CODE HERE #
    #########################
    seq1 = list(sequences[0].Sequence)
    seq2 = list(sequences[1].Sequence)
    diagonal = float('-inf')
    horizontal = float('-inf')
    vertical = float('-inf')
    t_matrix = [[0 for x in range(0, len(seq2) + 1)] for y in range(0, len(seq1) + 1)]  # creates a 0 filled matrix which can be filled in via scoring.
    #print_matrix_on_screen(t_matrix)

    s_matrix = [[0 for x in range(0, len(seq2) + 1)] for y in range(0, len(seq1) + 1)] #creates a 0 filled matrix which can be filled in via scoring.
    for y in range(0, len(seq1)+1):  #iterates through row numbers
        for x in range(0, len(seq2)+1): #iterates through the elements of each row.
             s_matrix[0][0] = 0 #sets the first cell of matrix as 0.
             s_matrix[0][x] = 0 - (x * penalty) #initialises the horizontal row with gap penalties.
             s_matrix[y][0] = 0 - (y * penalty) #initialies the vertical with gap penalties.
             if x >= 1 and y >= 1: #Ensures that the 0 cell at the start of the matrix is not iterated over.
                 diagonal = s_matrix[y - 1][x - 1] + matrix[ord(seq1[y - 1]) - ord('A')][ord(seq2[x - 1]) - ord('A')] #
             if x >= 1:
                vertical = s_matrix[y-1][x] - penalty #iterates through
             if y >= 1:
                 horizontal = s_matrix[y][x-1] - penalty
             maximum_score = max(diagonal, horizontal, vertical)
             t_matrix[0][x] = ('<')
             t_matrix[y][0] = ('^')
             t_matrix[0][0] = 0
             if x >= 1 and y >= 1:
                if horizontal == maximum_score:
                    t_matrix[y][x]= '<' #[y] then [x] is row then column.
                if diagonal == maximum_score:
                    t_matrix[y][x] = '\\'
                if vertical == maximum_score:
                    t_matrix[y][x] = '^'
             s_matrix[y][x] = maximum_score
    print_matrix_on_screen(t_matrix)
    seq1.insert(0, '-')
    seq2.insert(0, '-')
    seq1.insert(0, '')#creates the complete sequence with - signs and whitespace.
    s_matrix.insert(0, seq2)
    new_list=[]
    for i in range(0, len(s_matrix)):#iterates through each row number
        new_list.append(seq1[i])#creates a new list with each element of sequence 1 in.
        s_matrix[i].insert(0, new_list[i])


    #####Traceback method####
    seq1=seq1[2:]
    seq2=seq2[2:]
    print(seq1, seq2)
    print(len(seq1),len(seq2),len(t_matrix),len(t_matrix[0]))

    x=len(seq1)
    y=len(seq2)
    string1=''
    string2=''

    while(t_matrix[x][y] is not 0): #ensures the zero cell is not iterated over.
        if(t_matrix[x][y]=='^'): #determines vertical
            string1=string1+seq1[x-1] #adds the element from the vertical and assigns gap to horizontal position.
            string2=string2+'-'
            x=x-1
        if(t_matrix[x][y]=='<'): #determines the direction is horizontal.
            string1=string1='-'
            string2=seq2[y-1] #adds the element from the vertical and assigns gap to horizontal position.
            y=y-1
        if(t_matrix[x][y]=='\\'):
            string1=string1+seq1[x-1] #adds the characters from both sequences to the strings as it is a match.
            string2=string2+seq2[y-1]
            x=x-1
            y=y-1 #decreases the value of x and y so that the while loop does not

    alignment1 = (string1[::-1]) #reverses the string as the traceback was calculated from bottom to top and right to left.
    alignment2 = (string2[::-1])
    print(alignment1)
    print(alignment2)
    final_list=[]
    string3 = '-' * len(alignment1)#alignment 1 and 2 same length, so creates empty string of equal length to alignments.
    for i,j,z in zip(alignment1, alignment2, string3):
        if i == j:
            print(i, j, z)
            final_list.append("|")
        elif i != j:
            final_list.append('')
    print("Hello",final_list) #correct positions for matching sequences. Need to get this into a string at same position and align in middle of string1 and 2.
    final_string = ' '.join(final_list)
    print(alignment1)
    print(final_string)
    print(alignment2)

            #string3.replace(string3[z],"|",)

            #z = "|"

    #print(string3)







    """
    print(seq1, seq2)
    x = len(seq1) #n of rows.
    y = len(seq2) #n of columns.
    print(x, y)
    print(t_matrix[x][y])
    string1 = ''
    string2 = '' #creates the empty strings for alignment.
    while t_matrix[x][y] != 0:
        print(x,y,string1, string2)
        sign= t_matrix[x][y]
        print(sign)
        if t_matrix[x][y] == "^":
            print("hallo", seq1[y+1])
            string1= string1+seq1[y+1]
            string2= string2+ "-"
            x= x-1 #each iteration decreases x by 1 so that is not stuck in infite while loop.
        if t_matrix[x][y] == "<":
            string1= string1+ "-"
            string2= string2+ seq2[y+1]
            y= y-1
        if t_matrix[x][y] == "\\":
            string1= string1+ seq1[y+1]
            string2= string2+ seq2[x+1]
            x= x-1
            y= y-1
        print(string1)
        print(string2)
    #string1= string1[::-1]
    #string2= string2[::-1]
    #print(string1, string2)

    """


    return "x", s_matrix


#3 if statements check if >1 >> code on thing below in equation.

    #########################
    #   END YOUR CODE HERE  #
    #########################


def do_local_alignment(sequences, matrix, penalty):
    """Do pairwise local alignment using DP."""
    #########################
    # INSERT YOUR CODE HERE #
    #########################

    #########################
    #   END YOUR CODE HERE  #
    #########################


def do_semiglobal_alignment(sequences, matrix, penalty):
    """Do pairwise semi-global alignment using DP."""
    #########################
    # INSERT YOUR CODE HERE #
    #########################

    #########################
    #   END YOUR CODE HERE  #
    #########################


def print_matrix_to_file(matrix, fileName):
    """Write a matrix into file.

    matrix: a list of list in Python, storing a score
    matrix.
    fileName: str, a file name (with a path) to store the matrix.
    It is not recommended to tinker with this function.
    """
    with open(fileName, 'w') as f:
        for row in matrix:
            print('\t'.join(map(str, row)), file=f)


def print_alignment_to_file(alig, fileName):
    """Write a matrix into file.

    alig: a list of list in Python, storing an alignment.
    fileName: str, a file name (with a path) to store the alignment.
    It is not recommended to tinker with this function.
    """
    with open(fileName, 'w') as f:
        for row in alig:
            print(''.join(map(str, row)), file=f)


def print_matrix_on_screen(matrix, width=5):
    """Print a matrix on the screen.

    matrix: a list of list in Python, storing an alignment or a score
    matrix.
    width: that of the space one cell occupies.
    This will facilitate your testing.
    """
    for row in matrix:
        print(''.join(['{0:>{w}}'.format(item, w=width) for item in row]))


def main():
    """Main function.

    Please change it accordingly to make the program work.
    """
    # Get command line options
    args = get_args()

    # Set substitution matrix:
    if args.exchange_matrix == "pam250":
        exchangeMatrix = pam250
    elif args.exchange_matrix == "blosum62":
        exchangeMatrix = blosum62
    else:
        exchangeMatrix = identity

    # Read sequences from fasta file, and catch error reading file
    try:
        sequences = readSequences(open(args.fasta))
    except OSError as e:
        print("ERROR: cannot open or read fasta input file:", e.filename)

    for seq in sequences:
        print(seq)

    # Call alignment routine(s):
    if args.align_global:
        alignment, score_matrix = do_global_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    elif args.align_local:
        alignment, score_matrix = do_local_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    elif args.align_semiglobal:
        alignment, score_matrix = do_semiglobal_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    else:
        sys.exit("BUG! this should not happen.")


    #Print the result to files
    if args.alignment:
        print_alignment_to_file(alignment, args.alignment)
    if args.score_matrix:
        print_matrix_to_file(score_matrix, args.score_matrix)

    # Print the result on screen
    if args.print_on_screen:
        print_matrix_on_screen(alignment)
        print_matrix_on_screen(score_matrix)


if __name__ == "__main__":
    main()

# last line
