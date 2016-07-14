from scipy import stats

def main():
	fin = open('res.txt','r')
	fout1 = open('critical.txt','w')
	fout2 = open('all.txt','w')

	lines = fin.readlines()
	fin.close()

	results = []
	for line in lines:
		a = line.split(' ')

		if len(a) < 3:
			continue

		name = a[0]
		motif = a[1]
		args = a[2].split('(')[1].replace(')','').split(',')

		pvalue = stats.binom_test(int(args[0]), int(args[1]), float(args[2]))

		print name,motif,pvalue
		results.append( (name,motif,pvalue) )

	results = sorted(results, key=lambda x: float(x[2]))
	for x in results:
		name = x[0]
		motif = x[1]
		pvalue = x[2]

		#if below certain threshold, output to critical.txt
		if pvalue < 1e-7:
			fout1.write(str(name)+' '+str(motif)+' '+str(pvalue)+'\n')

		fout2.write(str(name)+' '+str(motif)+' '+str(pvalue)+'\n')

	fout1.close()
	fout2.close()
	
	pass

if __name__ == '__main__':
	main()