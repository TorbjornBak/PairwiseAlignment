#!/usr/bin/env python3


def matrixSubs(subsFile):
    #Function creating substitution dictionary for substitution matrix
    subsDict = {}
    infile = open(subsFile,'r')
    
    #Loop through substitution matrix file
    lineY = 0 
    for line in infile:
        #Save contents of first line (will be used as key)
        if lineY == 0:
            line = line.split()
            key = "".join(line)
        #Following lines must contain the values for the dict
        else: 
            line = line.split()
            lineX = 0
            #Loop throug all elements in line
            for i in line:
                #The first element is not a value and will be ignored
                if lineX != 0:
                    #The dict is updated with the pair keys and their corresponding value
                    subsDict[key[lineY-1],key[lineX-1]] = int(i)
                lineX += 1
        lineY += 1
    infile.close()
    return subsDict

