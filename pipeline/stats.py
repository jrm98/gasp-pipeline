# bio stats

from math import *
from compareTool import *
#from scipy import stats

def pval(motif,sequences):
	sInstances = {0:[], -1:[]}
	threshold = 0
	sequenceName = ""
	mLen = len(motif)

	with open(sequences, "r") as f:
		fileText = f.read().split("\n")
	fLen = len(fileText)
	lineNum = -1

	num_seq = 0
	num_spots = 0
	seq_matched = {}

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
		
		num_seq += 1
		sLen = len(sequence)
		
		# search sequence
		sPos = 0
		while sPos < sLen - mLen:
			num_spots += 1
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

	for key in sInstances:
		seq_matched[key] = []
		for i in sInstances[key]:
			if i[2] not in seq_matched:
				seq_matched[key] += [i[2]]


	matches = seq_matched[max(seq_matched.keys())] # gets best score
	p = float(len(matches) / float(num_spots))
	k = len(set([x[2] for x in matches]))
	n = num_seq


	print('k='+str(k)+'\tn='+str(n)+'\tP='+str(p))
	#return binomial_pdf(k,n,p)
	#return stats.binom_test(k,n,p)
	return 'stats.binom_test('+str(k)+','+str(n)+','+str(p)+')' #required libs not on rlogin, so just give parameters

def comb(n,r):
    f = factorial
    return f(n) / f(r) / f(n-r)

def normal_pdf(x, m, v):
	return 1.0/sqrt(2*pi*v) * exp(-(x-m)**2/(2*v))

def binomial_pdf(k, n, p):
	if n < 100:
		return comb(n, k) * p**k * p**(n-k)  # Fall back to your current method
	return normal_pdf(k, n*p, n*p*(1.0-p))


