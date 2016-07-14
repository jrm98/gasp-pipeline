#! /usr/bin/python

import threading
import subprocess
import time
import argparse
import fasta2bp
from compareTool import compare
from searchTool import searchFasta
from visTool import visualizeOutput

# define constants and paths, could be read from command line
NUM_OF_TOOLS = 5
MOTIF_LEN = 12
NUM_MOTIFS = 1 #number of motifs expected per sequence
DECOD_ITER = 5 #number of iterations to perform
SEQ_CARD = 5 #****number of total motif occurences (?)****

DATA_DIR = "data/"
POS_SEQ = "mad2.txt"
NEG_SEQ = "negative_mad2.txt"
RES_DIR = "results/"

times = []
foundMotifs = {}
mlist = []

# functions for running each tool
########################

def runCMF():
    times[0]["Name"] = "CMF"
    times[0]["Start"] = time.time()
    with open("results/cmf.out", "w") as fout:
        cmdLine = "./cmf/cmf -i1 " + POS_SEQ + " -i2 "\
              + NEG_SEQ + " -l " + str(MOTIF_LEN - 2)\
              + " -u " + str(MOTIF_LEN + 2) + " -w " + str(MOTIF_LEN - 2)\
              + " -f results -o results/cmfSeeds"
        subprocess.call([cmdLine], stdout=fout, shell=True)
    times[0]["End"] = time.time()
    print times[0]

def runWeeder():
    times[1]["Name"] = "Weeder"
    times[1]["Start"] = time.time()
    with open(RES_DIR+"weeder.out", "w") as fout:
        cmdLine = "./weeder/weeder2 -f " + POS_SEQ + " -O SP -maxm 10"
        subprocess.call([cmdLine], stdout=fout, stderr=fout, shell=True)
    times[1]["End"] = time.time()
    print times[1]

def runDECOD():
    times[2]["Name"] = "DECOD"
    times[2]["Start"] = time.time()
    with open("results/decod.out", "w") as fout:
        cmdLine = "java -jar ./DECOD/DECOD-20111024.jar -nogui -pos " \
            + POS_SEQ + " -neg " + NEG_SEQ + " -w " + str(MOTIF_LEN)\
            + " -nmotif " + str(NUM_MOTIFS) + " -niter " + str(DECOD_ITER) + " -c "\
            + str(SEQ_CARD) + " -o " + RES_DIR + "decod_found_motifs.txt"
        subprocess.call([cmdLine], stdout=fout, shell=True)
    times[2]["End"] = time.time()
    print times[2]

def runMEME():
    times[3]["Name"] = "MEME"
    times[3]["Start"] = time.time()
    with open("results/meme.out", "w") as fout:
        cmdline = "./meme/bin/meme -dna -oc " + RES_DIR + "meme " + POS_SEQ\
             + " -nmotifs " + str(NUM_MOTIFS) + " -w " + str(MOTIF_LEN)
        subprocess.call([cmdline], stdout=fout, stderr=fout, shell=True)
    times[3]["End"] = time.time()
    print times[3]

def runBioProspector():
    times[4]["Name"] = "BioProspector"
    times[4]["Start"] = time.time()
    with open("results/bioprospector.out", "w") as fout:
        outfile = fasta2bp.convert(POS_SEQ, "bp_"+get_filename(POS_SEQ))
        cmdline = "./BioProspector/BioProspector -i "+outfile \
            + " -o " + RES_DIR + "bp_output.txt"
        subprocess.call([cmdline], stdout=fout, stderr=fout, shell=True)
    times[4]["End"] = time.time()
    print times[4]

# runs all tools in parallel
########################

def runTools():
    global times
    times = [{} for x in range(NUM_OF_TOOLS)]
    t = [threading.Thread() for x in range(NUM_OF_TOOLS)]
    t[0].run = runCMF
    t[1].run = runWeeder
    t[2].run = runDECOD
    t[3].run = runMEME
    t[4].run = runBioProspector
    #t[5].run = runDREME
    # DREME hasn't been producing any results, and we aren't parsing it
    for th in t:
        th.start()
    print "Started pipeline"
    for th in t:
        th.join()
    print "Finished pipeline"



# parses the output from each tool and combines it into consistantly formatted file
########################

def parseCMF(fout):
    with open(RES_DIR + "output.txt", "r") as CMFResults:
        foundMotifs["CMF"] = []
        alreadyFound = {}
        for line in CMFResults:
            if line[0:7] == "MOTIF:\t" and not line[7:] in alreadyFound:
                    fout.write("CMF\t"+line[7:-1]+"\n")
                    foundMotifs["CMF"] += [line[7:-1]]
                    alreadyFound[line[7:]] = 0

def parseWeeder(fout):
    print "--|"+RES_DIR+POS_SEQ+".w2"
    with open(RES_DIR + POS_SEQ + ".w2", "r") as _weederResults:
        foundMotifs["Weeder"] = []
        weederResults = _weederResults.readlines()[6:]
        for result in weederResults:
            result = result.split()
            if len(result) == 0:
                break
            fout.write("weeder\t"+result[1]+"\n")
            foundMotifs["Weeder"] += [result[1]]

def parseDECOD(fout):
    with open(RES_DIR + "decod_found_motifs.txt", "r") as _decodResults:
        foundMotifs["DECOD"] = []
        alreadyFound = {}
        decodResults = _decodResults.readlines()
        chars = ["A", "C", "G", "T"]
        for motif_num in range(NUM_MOTIFS):
            PWM = {}
            for line in decodResults:
                if len(line) > 1 and line[0] in chars:
                    PWM[line[0]] = line.strip("ACGT []\n").split()
                    if line[0] == "T":
                        break
            motif = ""
            for i in range(MOTIF_LEN):
                col = [PWM[x][i] for x in chars]
                motif += chars[col.index(max(col))]
            if not motif in alreadyFound:
                fout.write("DECOD\t"+motif+"\n")
                foundMotifs["DECOD"] += [motif]
                alreadyFound[motif] = 0
#        currentMotifNum = "0"
#        for line in decodResults:
#            if line[0:6] == ">Motif":
#                currentMotifNum = line[6]
#            # Find only motif occurrences in positive sequences
#            # Negative sequences are preceded by gi|
#            # >Motif indicates the PWM
#            elif line[0] == ">" and line[1:4] != "gi|":
#                line = line.split()
#                fout.write("DECOD\t"+line[6]+"\n")


def parseMEME(fout):
    with open(RES_DIR + "/meme/meme.txt", "r") as memeResults:
        foundMotifs["MEME"] = []
        skip = -1
        for line in memeResults:
            # if skip > 0:
            #     skip -= 1
            #     continue
            # elif skip == 0:
            #     skip -= 1
            #     fout.write("MEME\t"+line+"\n")
            # print '"'+line[7:12]+'"'
            # if "Motif" in line and "regular expression" in line:
            #     skip = 1
            if "Multilevel" in line:
                motif = line.strip().replace("Multilevel","").replace(" ","")
                fout.write("MEME\t"+motif+"\n")
                foundMotifs["MEME"] += [motif]

def parseBioProspector(fout):
    with open(RES_DIR + "bp_output.txt", "r") as bpResults:
        foundMotifs["BioProspector"] = []
        for line in bpResults:
            if "Motif #" in line:
                motif = line[line.index("(")+1:line.index("/")]
                fout.write("BioProspector\t"+motif+"\n")
                foundMotifs["BioProspector"] += [motif]


def parseResults():
    with open("results/results.out", "w") as fout:
        parseCMF(fout)
        parseWeeder(fout)
        parseDECOD(fout)
        parseMEME(fout)
        parseBioProspector(fout)
        parseDREME(fout)

def combineMotifs():
    print "Running"
    results = []
    for tool in foundMotifs:
        others = list(set(foundMotifs) - set([tool]))
        for motif in foundMotifs[tool]:
            mscore = 0
            for oTool in others:
                for oMotif in foundMotifs[oTool]:
                    mscore += compare(motif, oMotif)
            results += [[mscore, motif]]
    results.sort()
    for m in results:
        print m[1], m[0]
    # search for top scoring motif
    mlist.append(results[-1][1])
    searchFasta(POS_SEQ, RES_DIR, results[-1][1], 1)

def get_filename(s):
    index = s.rfind('/')
    if index < 0:
        return s
    else:
        return s[index+1:]

def main():
    opt = argparse.ArgumentParser()
    opt.add_argument('positive_seq', help='the fasta file for the positive sequences')
    opt.add_argument('negative_seq', help='the fasta file for the negative sequences')
    opt.add_argument('-w', default=12, type=int, help='motif width')
    #opt.add_argument('-o', default='', type=str, help='output directory')
    opt.add_argument('--nmotifs', default=1, type=int, help='max number of distinct motifs to look for')
    opt.add_argument('--iter', default=5, type=int, help='max number of iterations (for DECOD tool)')
    args=opt.parse_args()



    global MOTIF_LEN
    global NUM_MOTIFS
    global DECOD_ITER 
    global SEQ_CARD
    global POS_SEQ
    global NEG_SEQ

    MOTIF_LEN = opt.w
    NUM_MOTIFS = opt.nmotifs #number of motifs expected per sequence
    DECOD_ITER = opt.iter #number of iterations to perform
    SEQ_CARD = 5 #****number of total motif occurences (?)****
    POS_SEQ = opt.positive_seq # "mad2.txt"
    NEG_SEQ = opt.negative_seq # "negative_mad2.txt"


    runTools()
    parseResults()
    print [tool["Name"] + " " + str(int(tool["End"]) - int(tool["Start"])) for tool in times]
    print foundMotifs
    combineMotifs()
    visualizeOutput(RES_DIR, mlist)

if __name__ == "__main__":
    main()
