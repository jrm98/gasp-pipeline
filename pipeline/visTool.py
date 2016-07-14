#! /usr/bin/python
import subprocess
from searchTool import makePWM
from stats import *

def createSeqLogo(RES_DIR, mnum):
    print "Starting visualization"
    cmdline = "./weblogo/seqlogo -F PNG -c -f " +\
        RES_DIR + "motif" + str(mnum) + ".instances" + " > " + RES_DIR +\
        "motif_result_" + str(mnum) + ".png"
    subprocess.call([cmdline], shell=True)

def parseInstances(RES_DIR, mnum):
    instances = {}
    maxSeqLen = 0
    with open(RES_DIR + "motif"+str(mnum+1)+".instances", "r") as instanceFile:
        for line in instanceFile:
            if len(line) < 2:
                pass
            elif line[0] == '>':
                sequenceName = line[1:line.rfind("\t")]
                seqPos, seqLen = line.split("\t")[-1].split("|")
                maxSeqLen = max(maxSeqLen, int(seqLen))
                if not sequenceName in instances:
                    instances[sequenceName] = []
            else:
                instances[sequenceName].append(line.strip()+"\t"+seqPos+'|'+seqLen)
    for i in instances:
        print i, instances[i]
    return (maxSeqLen, instances)

def createHTML(RES_DIR, mlist, pv):
    with open(RES_DIR+"results.html", "w") as rfile:
        rfile.write('<html><head><link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"></head>\n<body>\n')
        idx = 0
        for mnum, motif in enumerate(mlist):
            rfile.write('<div class="container"><h1>Motif: ' + motif + "</h1>\n")
            rfile.write("length: " + str(len(motif))+"<br>")
            rfile.write("p-value: " + str(pv[idx])+"<br>")
            idx += 1
            rfile.write('<img src="motif_result_'+str(mnum+1)+'.png"></div>\n<br>\n<br>')
            rfile.write('<div class="container"><h2>Position Weight Matrix</h2>\n')
            rfile.write('<pre style="margin-top: 0px; margin-bottom: 0px;">')
            maxLen, instances = parseInstances(RES_DIR, mnum)
            ins = []
            for i in instances:
                ins += [x.split("\t")[0] for x in instances[i]]
            makePWM(len(ins[0]), ins, rfile)
            rfile.write("</pre></div>\n")
            rfile.write('<div class="container"><h2>' + str(len(ins)) + " instances found</h2>\n")
            for sequence in instances:
                rfile.write("<h4>"+sequence + "</h4>\n")
                for instance in instances[sequence]:
                    insPos, seqLen = instance.split("\t")[-1].split('|')
                    seqLen = int(seqLen)
                    insPos = seqLen - int(insPos)
                    vis = (seqLen)/10 * "_"
                    vis = vis[:insPos/10] + '|' + vis[insPos/10+1:]
                    vis = (maxLen - seqLen)/10 * "&nbsp" + vis
                    rfile.write(instance[:instance.rfind("\t")] + "\n<br>")
                    rfile.write('<pre style="margin-top: 0px; margin-bottom: 0px;">')
                    rfile.write(vis+" + "+str(insPos)+" - "+str(seqLen - insPos)+" | "+str(seqLen))
                    rfile.write("</pre>\n<br>\n")
                rfile.write("\n<br>")
            rfile.write("</div>")
        rfile.write("</body>\n</html>")

def recordResult(filename, motif, pvalue, output):
    record = open(filename, 'a')
    record.write(output+' '+motif+' '+str(pvalue)+'\n')
    record.close()

def visualizeOutput(RES_DIR, mlist, filename, output):
    pv = [0 for x in mlist]
    idx = 0
    for mnum, motif in enumerate(mlist):
        pv[idx] = pval(motif, filename)
        createSeqLogo(RES_DIR, mnum+1)
        #if pv[idx] < 1e-7:
        recordResult('res.txt', motif, pv[idx], output)
        idx += 1
    createHTML(RES_DIR, mlist, pv)

def main():
    createHTML("results/", ["ACGT", "ACTG", "AGCT", "AGTC"])

if __name__ == "__main__":
    main()
