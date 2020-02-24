'''
For given motif sequence and PSSM, the script
calcualte probablity of that motif sequence
'''
import argparse 
import pandas as pd
import math
def iupac():
	# Nucleotide ambiguity code
	#code = {}
	code={
		"A":['A'], 
		"U":['U'], 
		"G":['G'], 
		"C":['C'], 
		"Y":['C', 'U'], 
		"R":['A', 'G'],
		"W":['A', 'U'],
		"S":['G', 'C'],
		"K":['U', 'G'],
		"M":['C', 'A'],
		"D":['A', 'G', 'U'],
		"V":['A', 'C', 'G'],
		"H":['A', 'C', 'U'],
		"B":['C', 'G', 'U'],
		"X":['A', 'C', 'G', 'U'],
		"N":['A', 'C', 'G', 'U']
	}
	return code
	#print len(code['N'])

def PSSM(pssm_file, motif):
	code = iupac() # define IUPAC 
	motif = motif.upper()
	pssm = pd.read_table(pssm_file)
	pssm = pssm.iloc[:,1:]

	# define consensus seq based the the base with highest prob in PSSM
	bases = pssm.idxmax(axis=1).tolist()	
	cons = ''.join(bases)

	# determine seq prob
	p_motif = p_cons = 1 
	for pos in range(len(motif)):
		# motif seq
		motif_base = motif[pos] # nt code
		motif_bases = code[motif_base] # base(s) for this code 
		p = pssm[motif_bases].iloc[pos].sum() # sum the prob if 2 or more bases for the same code
		p_motif *= p
		# cons seq
		cons_base = cons[pos]
		p = pssm[cons_base].iloc[pos]
		p_cons *= p 	
	print motif + '\t' +  str(p_motif) + '\t' + cons + '\t' + str(p_cons)#+  str((math.log10(p_motif))*-1)	

def main():
	# defing input/outout parametes
	parser = argparse.ArgumentParser()
	parser.add_argument("pssm", help="pssm file")
	parser.add_argument("motif", help="motif sequence")
	args = parser.parse_args()
	PSSM(args.pssm, args.motif)
	
if __name__ == '__main__':
	main()
