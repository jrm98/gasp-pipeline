#! /usr/bin/python
# batch tester for real data
# by: jake martinez (jrm98)

import subprocess
import os

def test(w, posseq, negseq='data/n_genes.fasta'):
	cmdLine = "python orange_pipeline.py -w "+str(w)+" "+posseq+" "+negseq
	subprocess.call([cmdLine], shell=True)

	f = open('results/results.html','r')
	html = f.read()
	f.close()
	f = open('results/motif_result_1.png','r')
	img1 = f.read()
	f.close()
	f = open('results/motif_result_2.png','r')
	img2 = f.read()
	f.close()
	f = open('results/motif_result_3.png','r')
	img3 = f.read()
	f.close()
	f = open('results/motif1.instances','r')
	inst1 = f.read()
	f.close()
	f = open('results/motif2.instances','r')
	inst2 = f.read()
	f.close()
	f = open('results/motif3.instances','r')
	inst3 = f.read()
	f.close()

	directory = 'results/res/'+get_filename(posseq)+'_'+str(w)+'/'
	try:
		os.mkdir(directory)
	except OSError:
		pass
	copy = open(directory+'results.html','w')
	copy.write(html)
	copy.close()
	copy = open(directory+'motif_result_1.png','w')
	copy.write(img1)
	copy.close()
	copy = open(directory+'motif_result_2.png','w')
	copy.write(img2)
	copy.close()
	copy = open(directory+'motif_result_3.png','w')
	copy.write(img3)
	copy.close()
	copy = open(directory+'motif1.instances','w')
	copy.write(inst1)
	copy.close()
	copy = open(directory+'motif2.instances','w')
	copy.write(inst2)
	copy.close()
	copy = open(directory+'motif3.instances','w')
	copy.write(inst3)
	copy.close()

	#see if we found a motif that was generated
	return parse_results()



def parse_results():
	try:
		with open('results/results.html','r') as f:
			lines = f.readlines()
			results = []
			for line in lines:
				if 'Motif:' in line:
					m = line.replace('Motif:','').replace(' ','').replace('\n','').replace('<divclass="container"><h1>','').replace('</h1>','')
			 		results += [m]
			return results
	except Exception:
		print 'file not found...'
	return results

def get_filename(s):
	index = s.rfind('/')
	end = s.rfind('.')
	if index < 0:
		return s
	else:
		return s[index+1:end]

def main():
	print('beginning test...')

	fout = open('output_real.txt','w')

	for w in range(6,16):
		print('  w='+str(w))
		#test results
		res1 = []
		res2 = []

		species = ['dm','hs','pneumocystis','scerevisiae','scryophilus',
			'sjaponicus','soctosporus','sp']

		#species = ['sp']
		genes = ['mad1','mad2','mad3','bub1','bub3','mph1']

		for s in species:
			res1 += [test(w, 'data/'+s+'.fasta')]
		print('by species: '+str(res1))

		for g in genes:
			res2 += [test(w, 'data/'+g+'.fasta')]
		print('by gene: '+str(res2))

		fout.write('w='+str(w)+'\n')
		fout.write(str(res1)+'\n')
		fout.write(str(res2)+'\n')

	fout.close()
	
	pass

if __name__ == '__main__':
	main()