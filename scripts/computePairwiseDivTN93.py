import argparse, sys, os
from Bio import AlignIO, SeqIO
from tn93 import tn93

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument("input_alignment")
	parser.add_argument("--mode", nargs="?", const=1, default=2)
	parser.add_argument("--list_addInfo", nargs="?", const=1, default=False)
	args = parser.parse_args()


	if args.input_alignment in ["-","stdin"]:
		infile = sys.stdin
	else:
		infile = args.input_alignment

	alignment = AlignIO.read(infile, "fasta")

	if args.list_addInfo:
		l_addinfo = eval(args.list_addInfo)
	
	alignment1 = str(alignment[0].seq.upper())
	alignment2 = str(alignment[1].seq.upper())

	l_out = []
	propGC_seq1 = round((alignment1.count("G") + alignment1.count("C")) / (len(alignment1) - alignment1.count("-")), 3)
	propGC_seq2 = round((alignment2.count("G") + alignment2.count("C")) / (len(alignment2) - alignment2.count("-")), 3)
	propGC = round((propGC_seq1 + propGC_seq2)/2, 3)
	distance_tn93 = tn93(alignment1, alignment2, len(alignment1), args.mode, 100)
	l_out.extend(l_addinfo)
	l_out.extend([propGC_seq1, propGC_seq2, propGC, distance_tn93])
	print("\t".join([str(x) for x in l_out]))

	


