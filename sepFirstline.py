#!/usr/bin/env python
import sys,getopt

def main(argv):
	print('seprate the first line with \\t, saved as names.txt. the first element is not includedv \n')
	inputfile=""
	outputfile=""
	try:
		opts, args = getopt.getopt(argv, "hi:o:",["infile=", "outfile="])
	except getopt.GetoptError:
		print('Error: test_arg.py -i <inputfile> -o <outputfile>')
		print('or: test_arg.py --infile=<inputfile> --outfile=<outputfile>')
		sys.exit(2)

	for opt, arg in opts:
		if opt == "-h":
			print('test_arg.py -i <inputfile> -o <outputfile>')
			print('or: test_arg.py --infile=<inputfile> --outfile=<outputfile>')
		 	sys.exit()
		elif opt in ("-i", "--infile"):
			inputfile = arg
		elif opt in ("-o", "--outfile"):
			outputfile = arg
	print(inputfile)
	print(outputfile)

	#print('seprate the first line with \\t, saved as names.txt. the first element is not included')
	#'/Users/jiajiepeng/workspace/p1/genes.xls'
	f = open(inputfile)
	firstLine = f.readlines()[0]
	ids = firstLine.split('\t')
	length_ids = len(ids)
	with open(outputfile,'w') as outf:
		for x in range(2,length_ids):
			outf.write("%s\n" % ids[x])

if __name__ == "__main__":
	main(sys.argv[1:])
