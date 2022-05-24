import argparse, sys
from Bio import AlignIO, SeqIO

def hamming_distance(seq1, seq2):
	agreeBP = sum([base1 == base2 for base1, base2 in zip(seq1, seq2) if base1 != "-" and base2 != "-"])
	disagreeBP = sum([base1 != base2 for base1, base2 in zip(seq1, seq2) if base1 != "-" and base2 != "-"])

	return (agreeBP, disagreeBP)

def window(seq, width=500, step=100):
	seqlen = len(seq)
	for i in range(0, seqlen, step):
		if i + width > seqlen:
			j = seqlen
		else:
			j = i + width
		yield (i, j)
		if j == seqlen:
			break

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("input_alignment")
	parser.add_argument("--window", nargs="?", const=1, default=1000)
	parser.add_argument("--step", nargs="?", const=1, default=100)
	parser.add_argument("--list_addInfo", nargs="?", const=1, default=False)
	args = parser.parse_args()

	# Load the sequence alignment.

	if args.input_alignment in ["-","stdin"]:
		infile = sys.stdin
	else:
		infile = args.input_alignment

	alignment = AlignIO.read(infile, "fasta")

	if args.list_addInfo:
		l_addinfo = eval(args.list_addInfo)

	for win in window(range(alignment.get_alignment_length()), width=int(args.window), step=int(args.step)):
		start, end = win
		if set(alignment[0][start : end]) != {'-'} and set(alignment[1][start : end]) != {'-'}:
			agree_bp, disagree_bp = hamming_distance(alignment[0][start : end], alignment[1][start : end])
			total_bp = agree_bp + disagree_bp
			frac_identity = (int(total_bp) - int(disagree_bp)) / float(total_bp)
			list_tmp = [str(x) for x in [start, end, agree_bp, disagree_bp, frac_identity, "w%s_s%s" % (args.window, args.step)]]
			list_tmp.extend(l_addinfo)
			print('\t'.join(list_tmp))
		
	

