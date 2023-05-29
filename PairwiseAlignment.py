#!/usr/bin/env python3
import SubstitutionMatrix as SM 
import FastaReader as FR #Used for inputing and reading fastafiles for the program
import math
import sys


def matrixCreator(seqA,seqB):
    # Creates None-matrix with size of the two sequences seqA and seqB
    rows, cols = (len(seqA)+1, len(seqB)+1)
    matrix=[[None for _ in range(cols)] for _ in range(rows)]
    # Fills out the first column and first row with zeros inorder to initialise matrix
    for y in range(len(seqA)+1):
        matrix[y][0]=0
    for x in range(len(seqB)+1):
        matrix[0][x]=0
    return matrix

def matrixFiller(emptyMatrix):
    # Takes empty matrix and fills it with values according to the substitution dictionary
    # We start from pos 1,1 since the first row and column have already been filled
    # Boundary is defined for optimization
    boundary = 15 + int((1/2)*math.sqrt(math.sqrt((len(seqA)+1)*(len(seqB)+1))))
    gapCount = 0
    # The flag is true if the two sequences are "too" asymetric for optimization
    try:
        flag = (len(seqA)/len(seqB) > 1.05 or len(seqA)/len(seqB) < 1/1.05)
        for y in range(1,len(seqA)+1):
            for x in range(1,len(seqB)+1):
                # If the sequences are within the boundary or to assymetric for opptimization
                if y > x - boundary and y < x + boundary or flag: #Optimization step
                    gapPenalty = gapScore(gapCount)
                    # Calculates gap scores and substitution scores for each position in matrix
                    if emptyMatrix[y][x-1] is not None:
                        gapScore1 = emptyMatrix[y][x-1] - gapPenalty     
                    if emptyMatrix[y-1][x] is not None:
                        gapScore2 = emptyMatrix[y-1][x] - gapPenalty
                    subScore = emptyMatrix[y-1][x-1] + int(subDict[seqA[y-1],seqB[x-1]])
                    if gapScore1 > subScore or gapScore2 > subScore:
                        gapCount += 1
                    else:
                        gapCount = 0
                    # Finds the maximum of the three scores and assigns this as the score at the current position
                    newScore = max(0,gapScore1,subScore,gapScore2) 
                    emptyMatrix[y][x] = newScore
    except KeyError:
        print("The substitution matrix used is not valid for this alignment.")
        sys.exit(1)
    
    filledMatrix = emptyMatrix
    return filledMatrix

def traceBack(matrix):
    # Function taking filled out matrix and traces the alignment path back from the largest score in matrix
    traceList = []
    #Defines different variables
    largestNumber = 0
    alignmentScore = 0
    #Finds the largest number and assigns the position of this number to the first element in traceList
    for y in range(1,len(seqA)+1):
        for x in range(1,len(seqB)+1):
            if matrix[y][x] is not None and matrix[y][x] > largestNumber:
                largestNumber = matrix[y][x]
                traceList = [[y, x]]
                
    y = traceList[0][0]
    x = traceList[0][1]
    while matrix[y][x] > 0:
        # Finds the coordinate to the left, diagonal and up direction of the current coordinate (x,y)
        left = [matrix[y-1][x],[y-1,x]]
        up = [matrix[y][x-1],[y,x-1]]
        diagonal = [matrix[y-1][x-1],[y-1,x-1]]

        left0 = left[0]
        up0 = up[0]
        diagonal0 = diagonal[0]
        # Checks if diagonal is zero and breaks while loop if true
        if diagonal0 == 0:
            traceList.insert(0,[y-1,x-1]) 
            break
        
        # Series of elif statements that inserts coordinates of either the diagonal,
        # left or right to the first element in traceList and pushes the other elements.
        elif diagonal0 >= up0 and diagonal0 >= left0:
            traceList.insert(0,diagonal[1])
            alignmentScore += diagonal0
        elif left0 >= up0 and left0 >= diagonal0:
            traceList.insert(0,left[1])
            alignmentScore += left0
        elif up0 >= left0 and up0 >= diagonal0:
            traceList.insert(0,up[1])
            alignmentScore += up0
        
        y = traceList[0][0]
        x = traceList[0][1]

    return traceList, alignmentScore

def traceTranslator(traceList,seqA,seqB):
    # Translates the traceBack coordinates from the list to an alignment
    # Returns the two alignments, an identity score, identity, alignmentLength and totalGapCount
    
    #Initializing parameters for the Trace Translation
    alignmentA = ""
    alignmentB = ""
    identityScore = 0
    identity = ""
    alignmentLength = 0 
    totalGapCount = 0

    # The trace back considers two different base cases:
    ### Case 1: If the traceback at a given position was diagonal (substitution)
    ### Case 2: or if the two elements was not a diagonal (gaps)
    for i in range(1,len(traceList)):
        elem = traceList[i]
        
        #Case 1: Substitution
        if elem[0] - 1 == traceList[i-1][0] and elem[1] - 1 == traceList[i-1][1]:
            alignmentA += seqA[traceList[i][0]-1]
            alignmentB += seqB[traceList[i][1]-1]
            
            if seqA[traceList[i][0]-1] == seqB[traceList[i][1]-1]: # Identical 
                identityScore += 1
                identity += "|"
            else: # Substitution, the two bases/AA are not identical
                identity +=" " 
        
        #Case 2: Gaps either tracing up or to the left
        elif elem[0] - 1 == traceList[i-1][0]: 
            alignmentA += seqA[traceList[i][0]-1]
            alignmentB += "-" #Inserts gap
            identity +=" " 
            totalGapCount += 1
        elif elem[1] - 1 == traceList[i-1][1]:
            alignmentA += "-" #Inserts gap
            alignmentB += seqB[traceList[i][1]-1] #Inserts base / AA
            identity +=" "
            totalGapCount += 1 

        alignmentLength += 1
    
    return alignmentA, alignmentB, identityScore, identity, alignmentLength, totalGapCount

def subPenalty(seqA,seqB):
    # Function that calls a dict containing the substitution matrix
    subPen = subDict[(seqA,seqB)]
    return subPen

def gapScore(gaps):
    # Calculates an affine gap penalty with the following formula:
    gapPenalty = (gaps) * gapExtend + gapOpenPenalty
    return gapPenalty

### Global variables ###
inputs = FR.inputMenu() # Takes input from console
# Opens fasta file, outputs the two sequences and 0 or 1 telling if input was dna or protein
try:
    sequences = FR.fastaReader(inputs[0])  
except TypeError as error:
    print(str(error))
    sys.exit(1)
#Saves the two sequences to two variables seqA and seqB
seqA = sequences[0]
seqB = sequences[1]
dnaorprotein = sequences[2] # sequences[2] outputs 0 or 1 (0 is dna, 1 is protein)

#Initializes the standard open and extend penalties used below
stdOpenPenalty = 10
stdExtendPenalty = 0.5

# Assigns the correct gap penalties based on the user input and assigns 
# the standard if no input from user was given.
if len(inputs) > 2:
    try:
        # Checks that inputs[i] is not none, if it was None, it means that the user 
        # did not enter a proper input for the penalty and assigns std penalty.
        if inputs[1] is not None:
            gapOpenPenalty = float(inputs[1])
        else:
            gapOpenPenalty = stdOpenPenalty

        if inputs[2] is not None:
            gapExtend = float(inputs[2])
        else:
            gapExtend = stdExtendPenalty
        
        #Makes sure no negative or 0 gap penalties are accepted
        if gapOpenPenalty <= 0 or gapExtend <= 0: 
            raise ValueError
    except ValueError as error:
        print("You did not input a positive value, but something else...")
        print("Please make sure that you input only positive integers for gapOpenPenalty and gapExtend.")
        sys.exit(1)
    
else: # Standard penaltys if the 
    gapOpenPenalty = 10 #Applies both to DNA and protein sequence gaps
    gapExtend = 0.5

#Loading substitution matrix based on if the input sequences was DNA or protein
if len(inputs) > 2 and inputs[3] == 0: #DNA (overwritten by user)
    subDict = SM.matrixSubs("DNAFULL.txt")
    dnaorprotein = 0
elif len(inputs) > 2 and inputs[3] == 1: #Protein (overwritten by user)
    subDict = SM.matrixSubs("Blosum62.txt")
    dnaorprotein = 1
elif dnaorprotein == 0: #DNA as detected by the type algorithm
    subDict = SM.matrixSubs("DNAFULL.txt")
elif dnaorprotein == 1: #Protein as detected by the type algorithm
    subDict = SM.matrixSubs("Blosum62.txt")


### Main function of the program
def main():
    # Runs the two sequences through the 4 algorithms
    emptyMatrix = matrixCreator(seqA,seqB) # Creates an empty matrix
    scoringMatrix = matrixFiller(emptyMatrix) # Scores the matrix
    # Traces the alignment back and outputs the
    #  coordinates of the alignment and alignment score
    trace = traceBack(scoringMatrix) 
    traceCoord = trace[0] # Saves the coordinates to new variable
    result = (traceTranslator(traceCoord,seqA,seqB)) # Saves result of the alignment

    # Different alignment metrics derived from "result"
    identityScore = result[2] 
    alignmentLength = result[4]
    gapCount = result[5]

    # Calculates different scores for the alignment
    identityPercent = (identityScore / alignmentLength) * 100
    gapPercent = (gapCount / alignmentLength) * 100
    

    i = 0 # Initializing iterator starting at 0 for the while loop
    lineLength = 50 # Variable telling how many characters of the alignment should be printed on each line
    #While loop that prints the alignment to the console and to an output file
    while result[0][i*lineLength:(i+1)*lineLength] != "" or result[1][i*lineLength:(i+1)*lineLength] != "":
        # Gets the correct length of seq1 and 2 coords from the tracecoord
        seq1Coord = (traceCoord[i*lineLength][0]+1)
        seq2Coord = (traceCoord[i*lineLength][1]+1)
        
        #Tries to extract the end coordinates from trace
        try:
            seq1CoordEnd = (traceCoord[i*lineLength+lineLength][0])
            seq2CoordEnd = (traceCoord[i*lineLength+lineLength][1])
        except IndexError:  
        #If index exceeds trace list, the last line must have been reached 
        # thus print the last coordinates from trace
            seq1CoordEnd = traceCoord[-1][0]
            seq2CoordEnd = traceCoord[-1][1]

        #Prints to console
        print("Seq1:\t",seq1Coord," ",result[0][i*lineLength:(i+1)*lineLength]," ",seq1CoordEnd,sep="")
        print("Iden:\t",(len(str(seq1Coord))+1)*" ",result[3][i*lineLength:(i+1)*lineLength],sep="")
        print("Seq2:\t",seq2Coord, " ",result[1][i*lineLength:(i+1)*lineLength]," ",seq2CoordEnd,sep="")
        
        i += 1

    # Alignment information
    
    print("\nSeq1 ID: %s. Seq2 ID: %s." % (sequences[3],sequences[4])) #Seq ID's
    print("\nAlignment length: ",alignmentLength) #Alignment length
    print("Alignment Type:", "DNA alignment" if dnaorprotein == 0 else "AA alignment") # Alignment type
    print("Identity:    ", identityScore, "/", alignmentLength, " (",round(identityPercent,1),"%)", sep="") # Identity score
    print("Gap count:   ", gapCount, "/", alignmentLength, " (",round(gapPercent,1),"%)", sep="") # Gap count
    print("Gap penalties used: Open = ",round(gapOpenPenalty,1),", Extend = ",(round(gapExtend,1)),".",sep="") # Gap penalties used
    
main()