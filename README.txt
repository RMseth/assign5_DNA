#################
#     PreReqs
#################

Python 2.7
re - regex module
numpy - numpy module

###################
#  sequence files
###################

just copy and paste the code into a text file with all the numbers and spaces
my regex will take care of removing the numbers and spaces and then convert
it to 1 long string.

#################
#    DP O(n^2)
#################

If you uncomment these lines it will run the DP code

Matrix = dynamic(seq1, seq2)
score = Matrix[len(seq1)-1][len(seq2)-1]
traceBackStr = traceBack(Matrix, seq1, seq2)


#################
#  DivConq O(n)
#################

Each sequence has it's own call, that will run the algo and output the info to it's appropriate text file

h = divConq(sequence1, sequence1)
printing('seq1', 'seq2', h)

####################
#      report
####################

For the DP O(n^2) for my algorithm, I can only get about a 1/3 of the sequence before it starts to get too crazy
to keep using it.

For one of my DivConq algorithms it takes about 20min to get through and print a result.

human mitochondrial sequences:
 **** to see table go to "human" directory ****

