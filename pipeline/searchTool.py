#! /usr/bin/python
from compareTool import compare
from sys import stdout

def searchFasta(fName, resultsDir, motif, motifnum):
    print "~",motifnum,"#"
    sInstances = {0:[], -1:[]}
    threshold = 0
    sequenceName = ""
    mLen = len(motif)

    with open(fName, "r") as f:
        fileText = f.read().split("\n")
    fLen = len(fileText)
    lineNum = -1

    while lineNum < fLen - 1:
        # retrieve line
        lineNum += 1
        line = fileText[lineNum]

        if len(line) == 0:
            # ignore blank lines
            continue

        if line[0] == '>':
            # beginning of new sequence, save sequence name
            sequenceName = line[1:]
            continue
        
        # retrieve sequence
        sequence = ""
        while lineNum < fLen and len(fileText[lineNum]) != 0:
            sequence += fileText[lineNum]
            lineNum += 1
        
        sLen = len(sequence)
        
        # search sequence
        sPos = 0
        while sPos < sLen - mLen:
            match = sequence[sPos:sPos+mLen]
            score = compare(motif, match, vals=[-1, 0, 2])
            if score > min(sInstances.keys()):
                if not score in sInstances:
                    # add score to found instances and remove previous min
                    sInstances[score] = []
                    del sInstances[min(sInstances.keys())]
                sInstances[score].append([match, score, sequenceName, str(sLen-sPos), sLen])
                sPos += mLen
            #end if
            sPos += 1
        #end while
    #end while

    print motif + "\n Instances:"
    instances = []
    for key in sInstances:
        for instance in sInstances[key]:
            instances.append(instance)
    
    # write to file for visualisation
    with open(resultsDir + "motif" + str(motifnum) + ".instances", "w") as fout:
        for instance in instances:
            fout.write(">"+instance[2]+"\t"+str(instance[3])+"|"+str(instance[4])\
                    +"\n"+instance[0]+"\n\n")
    makePWM(mLen, [x[0] for x in instances], stdout)
    
def makePWM(mLen, instances, fout):
    print "~" + str(instances) + "!"
    pwm = [[0.0 for x in xrange(mLen)] for y in xrange(4)]
    letters = {'A':0, 'C':1, 'G':2, 'T':3}
    inc = 1.0/len(instances)
    
    for instance in instances:
        #fout.write(str(instance)+"\n")
        for pos, char in enumerate(instance):
            pwm[letters[char]][pos] += inc
    
    for let in "ACGT":
        row = let + ": "
        for posVal in pwm[letters[let]]:
            row += "{:.2f} | ".format(posVal)
        fout.write(row+"\n")
#        print pwm



def main():
    fastaFile = "data/positive_genes.fasta"
    searchFasta(fastaFile, "./", "ACTCTATG")

if __name__ == "__main__":
    main()
