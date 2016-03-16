import sys
import re
import numpy as np

##############################################
#           dynamic program O(n^2)
##############################################
def traceBack(Matrix, dna1, dna2):
    string1 = ''
    string2 = ''

    i = len(Matrix) - 1
    j = len(Matrix[0]) - 1

    finished = False
    while(finished == False):
        if i <= 0 or j <= 0:
            finished = True

        # if(i > 0 and j > 0):
        ifMatch = Matrix[i-1][j-1]
        ifDelete = Matrix[i][j-1]
        ifInsert = Matrix[i-1][j]

        if(ifMatch >= max(ifDelete, ifInsert)):
            # no change
            string1 = string1 + dna1[i-1]
            string2 = string2 + dna2[j-1]
            i = i-1
            j = j-1

        elif(ifInsert >= max(ifMatch, ifDelete)):
            # needed to insert
            string1 = string1 + dna1[i-1]
            string2 = string2 + '_'
            i = i-1

        elif(ifDelete >= max(ifMatch, ifInsert)):
            # needed to delete
            string1 = string1 + '_'
            string2 = string2 + dna2[j-1]
            j = j-1

        else:
            print('somethine when wrong')
            finished = True

    return string1, string2

################################################
#             NWScore
################################################
def NWScore(dna1,dna2):
    LastLine = []
    Matrix = dynamic(dna1, dna2)

    for j in range(len(dna1)):
        LastLine.append(Matrix[len(dna2)-1][j])
    return LastLine

################################################
#               PartitionY
################################################
def PartitionY(scoreL, scoreR):
    return np.argmax(scoreL + scoreR[::-1])

################################################
#             NeedlemanWunsch
################################################
def NeedlemanWunsch(dna1, dna2):

    Matrix = dynamic(dna1, dna2)
    return traceBack(Matrix, dna2, dna1)

################################################
#           Divide and Conquer O(n)
################################################
def divConq(X, Y):
    Z = ''
    W = ''

    if len(X) == 0:
        for i in range(len(Y)):
            Z = Z + '_'
            W = W + Y[i]
    elif len(Y) == 0:
        for i in range(len(X)):
            Z = Z + X[i]
            W = W + '_'
    elif len(X) == 1 or len(Y) == 1:
        Z, W = NeedlemanWunsch(X, Y)
    else:
        xmid     = len(X)/2
        xmid2    = X[xmid:]
        revXmid2 = xmid2[::-1]

        ScoreL = NWScore(X[:xmid], Y)
        ScoreR = NWScore(revXmid2, Y[::-1])
        ymid   = PartitionY(ScoreL, ScoreR)

        Z1, W1 = divConq(X[:xmid], Y[:ymid])
        Z2, W2 = divConq(X[xmid:], Y[ymid:])

        Z = Z1 + Z2
        W = W1 + W2

    return Z, W

##################################################
#           creates matrix  for algorithms
#               to use
##################################################
def dynamic(dna1, dna2):
    dna1 = list(dna1)
    dna2 = list(dna2)

    Matrix = [[0 for n in range(0, len(dna1)+1)] for m in range(0, len(dna2)+1)]

    for x in range(len(dna2)): # adds the 2nd dna sequence into the matrix (vertical x)
        Matrix[x][1] = dictionary(dna2[x-1], '_') + Matrix[x-1][0]

    for y in range(len(dna1)): # adds the 1st dan sequence into the matrix (horizontal y)
        Matrix[1][y] = dictionary('_', dna1[y-1]) + Matrix[0][y-1]

    for x in range(1, len(dna2)+1): # fills in the matrix with the comparison of dna1 to dna2
        for y in range(1,len(dna1)+1):
            # if((x-1) >= 0 and (y-1) >= 0):
            ifMatch = Matrix[x-1][y-1] + dictionary(dna2[x-1], dna1[y-1])
            ifDelete = Matrix[x][y-1] + dictionary('_', dna1[y-1])
            ifInsert = Matrix[x-1][y] + dictionary(dna2[x-1], '_')

            Matrix[x][y] = max(ifDelete, ifInsert, ifMatch)

    return Matrix

###############################################
#           Dictionary
###############################################
def dictionary(letter1, letter2):
    if(letter1 == letter2):
        return 5
    elif(letter1 == '_' and letter2 == '_'):
        return -6
    elif(letter1 == 'A'):
        if(letter2 == 'C'):
            return -1
        elif(letter2 == 'G'):
            return -2
        elif(letter2 == 'T'):
            return -1
        elif(letter2 == '_'):
            return -3
    elif(letter1 == 'C'):
        if(letter2 == 'A'):
            return -1
        elif(letter2 == 'G'):
            return -3
        elif(letter2 == 'T'):
            return -2
        elif(letter2 == '_'):
            return -4
    elif(letter1 == 'G'):
        if(letter2 == 'A'):
            return -2
        elif(letter2 == 'C'):
            return -3
        elif(letter2 == 'T'):
            return -2
        elif(letter2 == '_'):
            return -2
    elif(letter1 == 'T'):
        if(letter2 == 'A'):
            return -1
        elif(letter2 == 'C'):
            return -2
        elif(letter2 == 'G'):
            return -2
        elif(letter2 == '_'):
            return -1
    elif(letter1 == '_'):
        if(letter2 == 'A'):
            return -3
        elif(letter2 == 'C'):
            return -4
        elif(letter2 == 'G'):
            return -2
        elif(letter2 == 'T'):
            return -1

########################################################
#           read in sequences and store as lists
########################################################
# sequenceStrGorrilla = ""
# with open('gorillaMitochondrial.txt', 'r') as ins:
#     for line in ins:
#         m = re.search('(( ?\D{1,10}){1,6})', line)
#         sequenceStrGorrilla += m.group(0).rstrip('\n').replace(' ','')
#     # print("\n\nGorrilla\n" + sequenceStrGorrilla)
#
# sequenceStrHomoSapiens = ""
# with open('homoSapiensMitochondrion.txt', 'r') as ins:
#     for line in ins:
#         m = re.search('(( ?\D{1,10}){1,6})', line)
#         sequenceStrHomoSapiens += m.group(0).rstrip('\n').replace(' ','')
#     # print("\n\nhomoSapiense\n" + sequenceStrHomoSapiens)

# sequenceStrHomoSapiensEnglish = ""
# with open('homoSapiensMitochondrionEnglish.txt', 'r') as ins:
#     for line in ins:
#         m = re.search('(( ?\D{1,10}){1,6})', line)
#         sequenceStrHomoSapiensEnglish += m.group(0).rstrip('\n').replace(' ','')
#     # print("\n\nsequenceStrHomoSapiensEnglish\n" + sequenceStrHomoSapiensEnglish)
#
# sequenceStrNeanderthal = ""
# with open('homoSapiensNeanderthalensisMitochondrion.txt', 'r') as ins:
#     for line in ins:
#         m = re.search('(( ?\D{1,10}){1,6})', line)
#         sequenceStrNeanderthal += m.group(0).rstrip('\n').replace(' ','')
#     # print("\n\nsequenceStrNeanderthal\n" + sequenceStrNeanderthal)

seqTest1 = ""
with open('test1.txt', 'r') as ins:
    for line in ins:
        m = re.search('(( ?\D{1,10}){1,6})', line)
        seqTest1 += m.group(0).rstrip('\n').replace(' ','').upper()
    # print("\n\seqTest1\n" + seqTest1)

seqTest2 = ""
with open('test2.txt', 'r') as ins:
    for line in ins:
        m = re.search('(( ?\D{1,10}){1,6})', line)
        seqTest2 += m.group(0).rstrip('\n').replace(' ','').upper()
    # print("\n\seqTest2\n" + seqTest2)

####################################################
#       function calls
####################################################

Matrix = dynamic(seqTest1, seqTest2)
print('Dynamic: ', Matrix[len(seqTest1)-1][len(seqTest2)-1])
print('trackBack: ', traceBack(Matrix, seqTest1, seqTest2))

h = divConq(seqTest1, seqTest2)
print 'trackBack: ', h
score = 0
for i in range(0, len(h[0])):
    score = score + dictionary(h[0][i], h[1][i])

print 'DivConq: ', score

# print'changes', matrixFill(seqTest1, seqTest2)
# print'changes', matrixFillmatrixFill(sequenceStrGorrilla, sequenceStrHomoSapiens)
# print'changes', matrixFill(sequenceStrGorrilla, sequenceStrHomoSapiensEnglish)
# print'changes', matrixFill(sequenceStrGorrilla, sequenceStrNeanderthal)
# print'changes', matrixFill(sequenceStrHomoSapiens, sequenceStrHomoSapiensEnglish)
# print'changes', matrixFill(sequenceStrHomoSapiens, sequenceStrNeanderthal)
# print'changes', matrixFill(sequenceStrHomoSapiensEnglish, sequenceStrNeanderthal)