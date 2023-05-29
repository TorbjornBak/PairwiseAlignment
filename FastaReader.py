#!/usr/bin/env python3

import sys
import re
def inputMenu():
    #
    inputtypes = ["file for processing","Gap Open penalty (std: 10)",
    "Gap Extend Penalty (std: 0.5)","DNA or protein (enter for auto detect)"]
    
    inputcount = len(inputtypes)
    userinput = [] # Creates the empty userinput list
    if len(sys.argv) > inputcount + 1: #Checks if there was two many inputs
        print("You should enter only %d input%s!" % (inputcount,"s" if inputcount > 1 else ""))
        sys.exit(1) 
    
    #If only a file name has been given, the program assumes 
    # the user wishes the standard open penalty and extension penalty
    elif (len(sys.argv) == 2): 
        userinput.append(sys.argv[1])
    
    # Checks if the number of inputs was 4 
    elif (len(sys.argv) == inputcount): 
        for i in range(inputcount-1):
            userinput.append(sys.argv[i+1]) #Appends the command inputs one by one to the userinput list
        userinput.append(None) #Appends None to last place in list, makes sure we don't run in to index errors
     # Checks if the number of inputs was 5
    elif (len(sys.argv) == inputcount + 1): 
        for i in range(inputcount-1):
            userinput.append(sys.argv[i+1]) #Appends the command inputs one by one to the userinput list
        
        # Overwriting the automatic detection of protein or dna
        if (str(sys.argv[inputcount]).upper() == "PROTEIN"):
            userinput.append(1)
        elif str(sys.argv[inputcount]).upper() == "DNA":
            userinput.append(0)
        else:
            userinput.append(None)

    # Else the user is prompted to input filename, as well as the two penalty scores and the protein/dna overwrite
    else: 
        print("Please enter %d input%s..." % (inputcount,"s" if inputcount > 1 else ""))
        for i in range(inputcount):
            userInput = input("Please input %2s: " % (inputtypes[i]))
            if len(userInput) == 0:
                userInput = None
            userinput.append(userInput) #Appends the command inputs one by one to the userinput list
        
    return userinput # Returns list of userinputs
        

def fastaReader(file):
    # Reads a file containing two fasta files
    # and returns two sequences A and B along with seqStatus,
    # which is telling if the sequences were DNA or protein.
    # The function also returns the names of the two fasta entries.
    try:
        infile = open(file,'r')
    except FileNotFoundError:
        print("File not found!")
        sys.exit(1)
    except NameError:
        print("File not found!")
        sys.exit(1)
    except TypeError:
        print("File not inputted!")
        sys.exit(1)

    #Check lines for info using ">" as criteria.
    seq = ""
    seqCount = 0
    for line in infile:
        if line[0] == ">":
            seqCount += 1

            if seqCount > 2:
                print("ERROR in file, more than 2 sequences submitted!")
                sys.exit(1)

            if seq != "":
                A = seq
                seq = ""
            
            # Finds the name / acc. id of the fasta entry
            name = line.split(" ")[0].replace(">","").strip()
            if seqCount == 1:
                seqAName = name
            elif seqCount == 2:
                re.search
                seqBName = name

        #If line doesn't contain info other than DNA
        else:           
            #Following gets the dna seq
            seq += line.strip()
    infile.close()

    # for the last entry:        
    if seq != "":
        B = seq

    if seqCount < 2 :
        print("ERROR in file, less than 2 sequences submitted!")
        sys.exit(1)

    
    
    # Making two sets to check if strings are DNA or Protein and that the strings are compatible
    codons = set(["A","R","N","D","C","Q","E","G","H","I","L",
    "K","M","F","P","S","T","W","Y","V","B","Z","X"])
    dna = set(["A","T","C","G","N"])
    
    seq1 = set(A)
    seq2 = set(B)

    # Checks if dna or codons are superset of seq1 and 2
    # and stores boolean
    superDNA1 = dna.issuperset(seq1) 
    superDNA2 = dna.issuperset(seq2)
    superprot1 = codons.issuperset(seq1)
    superprot2 = codons.issuperset(seq2)

    if superDNA1 and superDNA2:
        #If the DNA set is a superset of both strings, the input must be DNA
        seqStatus = 0 
    elif superprot1 and superprot2:
        # If codons is a superset of both strings and the DNA set is not a 
        # superset of any of the two sequences, the input must be protein
        if superDNA1 or superDNA2: 
            # If DNA is a superset of either of the strings,
            # the two sequences cannot be aligned.
            print("ERROR in file, sequences should either be both protein or both DNA!")
            sys.exit(1)
        else: #The sequence is protein
            seqStatus = 1 
    else:
        # If the input strings are neither DNA or protein, we print the following...
        print("ERROR in file, sequences should either be both protein or both DNA!")
        sys.exit(1)
        
    
    return A,B,seqStatus,seqAName,seqBName







