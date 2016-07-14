# reformat fasta for bioprospector
# author: jake martinez
# python fasta2bp.py <input file> [<output file>]

import sys

def main():
	outfilename = 'bp_'+sys.argv[1]
	if len(sys.argv) >= 3:
		outfilename = sys.argv[2]
	convert(sys.argv[1],outfilename=outfilename)

def convert(infilename, outfilename=None):
	if outfilename is None:
		outfilename = 'bp_'+infilename
	fin = open(infilename,'r')
	fout = open(outfilename,'w')
	mappings = open(outfilename+'.mappings','w')
	lines = fin.read().split('\n')
	seqcount = 0
	for line in lines:
		if '>' in line:
			if seqcount > 0:
				fout.write('\n')
			fout.write('>seq'+str(seqcount)+'\n')
			mappings.write('seq'+str(seqcount)+'>'+line+'\n')
			seqcount += 1
		else:
			fout.write(line.lower())

	fin.close()
	fout.close()
	mappings.close()
	return outfilename

if __name__ == '__main__':
	main()