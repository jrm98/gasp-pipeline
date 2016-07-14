def parse_results():
	# try:
	f = open('results/results.html','r')
	print('read file...')
	lines = f.readlines()
	results = []
	for line in lines:
		if 'Motif:' in line:
			m = line.replace('Motif:','').replace(' ','').replace('\n','')
			results += [m]
	# except Exception:
	# 	print 'file not found...'
	return results

def main():
	print parse_results()
	pass

if __name__ == '__main__':
	main()