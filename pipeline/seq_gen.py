#! /usr/bin/python
# randomized motif-embedding sequence generator
# by: jake martinez (jrm98)


import random
import argparse
import time


def main():
	opt = argparse.ArgumentParser()
	opt.add_argument('-w', default=8, type=int, help='motif width')
	opt.add_argument('-N', default=20, type=int, help='number of sequences')
	opt.add_argument('-n', default=100, type=int, help='sequence length')
	opt.add_argument('-e', default=0, type=int, help='max # of SNPs in motif')
	opt.add_argument('--eprob', default=.65, type=float, help='probability of a SNP existing in the motif')
	opt.add_argument('--maxloc', default=-1, type=float, help='max ending location for motifs (-1 represents the end of the sequence)')
	opt.add_argument('--minloc', default=0, type=float, help='min starting location for motifs')
	opt.add_argument('-P', default=1, type=float, help='probability motif is present in sequence')
	opt.add_argument('-o', default='seq_'+str(int(time.time()))+'.fasta', type=str, help='output filename')
	opt.add_argument('--dir', default='', type=str, help='output directory')
	opt.add_argument('--dyad', help='spaced dyad mode (overrides nmotifs parameter, -w option becomes the total width of the combined dyad)', action='store_true')
	opt.add_argument('--nmotifs', default=1, type=int, help='max number of distinct motifs in a single sequence')
	opt.add_argument('--nsites', default=1, type=int, help='max number of appearances of the motif in a single sequence (not implemented)')
	opt.add_argument('--negative', help='used to generate a sequence for a negative discriminative set', action='store_true')
	args=opt.parse_args()

	N = args.N
	n = args.n
	e = args.e
	ep = args.eprob
	P = args.P
	w = args.w
	minloc = args.minloc
	maxloc = args.maxloc
	output = args.o
	directory = args.dir
	negative = args.negative
	nmotifs = args.nmotifs
	dyad = args.dyad

	generate(w,N,n,e,ep,P,minloc,maxloc,output,directory,negative,nmotifs,dyad)


def generate(w,N,n,e=0,ep=.65,P=1,minloc=0,maxloc=-1,
	output='seq_'+str(int(time.time()))+'.fasta',directory='',negative=False,
	nmotifs=1,dyad=False):
	if len(directory) > 0 and directory[-1:] != '/':
		directory += '/'

	if negative and not dyad:
		print('seq_gen: NEGATIVE=true')
	elif dyad:
		print('seq_gen: DYAD=true')
		print('seq_gen: w='+str(w)+'\tN='+str(N)+'\tn='+str(n)+'\te='+str(e)+'\tP='+str(P))
	else:
		print('seq_gen: w='+str(w)+'\tN='+str(N)+'\tn='+str(n)+'\te='+str(e)+'\tP='+str(P))
	print('\toutput => '+directory+output)

	sequences = []
	motif = []


	if dyad:
		w1 = w - random.randint(2,w-2) #requires total motif width of at least 4
		w2 = w - w1
		motif.append(gen_seq(w1))
		motif.append(gen_seq(w2))
		for i in range(N):
			sequence = gen_seq(n)
			result = embed_motif(sequence, motif, maxerror=e, errorprob=ep, P=P, minstart=minloc, maxstart=maxloc, nmotifs=2)
			sequences.append((result,motif,i))
			# print result, motif, i
		write_to_fasta(directory+output, sequences, negative)
		return

	for i in range(nmotifs):
		motif.append(gen_seq(w))
	for i in range(N):
		sequence = gen_seq(n)
		if negative:
			result = embed_motif(sequence, motif, maxerror=0, errorprob=0, P=0, minstart=minloc, maxstart=maxloc)
		else:
			result = embed_motif(sequence, motif, maxerror=e, errorprob=ep, P=P, minstart=minloc, maxstart=maxloc, nmotifs=nmotifs)
		sequences.append((result,motif,i))
		# print result, motif, i
	write_to_fasta(directory+output, sequences, negative)
	return (sequences,motif)

def gen_seq(length):
	return ''.join(random.SystemRandom().choice('ACGT') for _ in range(length))

def embed_motif(seq, motif, maxerror=0, errorprob=.20, minstart=0, maxstart=-1, P=1,nmotifs=1):
	if random.random() < P and (maxstart > minstart or maxstart == -1):
		if maxstart < 0:
			maxstart = len(seq) - len(motif[0]) - 1

		# pick starting location for motif
		loc = [random.randint(minstart,maxstart)]

		# introduce artificial error in motif
		errors = [0]
		errorind = [[]]
		r = random.random()
		while r < errorprob and errors[0] < maxerror:
			letters = 'ACGT'.replace(motif[0][int(r*len(motif[0]))],'')
			snp = random.randint(0,len(motif[0])-1)
			while snp in errorind[0]:
				snp = random.randint(0,len(motif[0])-1)
			motif[0] = motif[0][:snp] \
				+ random.SystemRandom().choice(letters) \
				+ motif[0][snp+1:]
			errors[0] += 1
			errorind[0].append(snp)
			r = random.random()

		# embed modified motif in sequence
		newseq = seq[:loc[0]] + motif[0] + seq[loc[0] + len(motif[0]):]
		
		errorind[0].sort() # sort snp locations in motif before returning

		# try to add more motifs
		if nmotifs > 1 and random.random() < P:
			res = embed_motif(newseq, motif[1:],maxerror,errorprob,loc[0]+len(motif[0]),maxstart,P,nmotifs-1)
			newseq = res[0]
			loc += res[1]
			errors += res[3]
			errorind += res[4]
		return (newseq, loc, motif, errors, errorind)
	else:
		return (seq, [-1], None, [None], [None])

def write_to_fasta(filename, sequences, negative):
	f = open(filename, 'w')
	for sequence in sequences:
		if negative:
			f.write('>artificial'+str(sequence[2])+'-'+str(int(time.time()))+' NEGATIVE sequence'+'\n')
		else:
			f.write('>artificial'+str(sequence[2])+'-'+str(int(time.time()))+' seq w/ embedded motif:'+str(sequence[1])+' start_loc:'+str(sequence[0][1])+' errors:'+str(sequence[0][3])+' errorind:'+str(sequence[0][4])+'\n')
		f.write(sequence[0][0]+'\n\n')
	f.close()
	pass

if __name__ == '__main__':
	main()