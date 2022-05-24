import os, re, glob, itertools
import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)

if config == {}:
	configfile: "%s/config.yaml" % SNAKEMAKE_DIR

if not os.path.exists("log"):
	os.makedirs("log")

dir_inFASTA = config["dir_inFASTA"]
list_inFASTA = glob.glob1(dir_inFASTA, "*.fasta")

list1_ids = config["list1_ids"]
list2_ids = config["list2_ids"]
		

rule dummy: 
	input: expand("identityCount/{id1}-{id2}.identityCount", id1 = list1_ids, id2=list2_ids),
			expand("summary_divTN93/{id1}.TN93.txt", id1=list1_ids)


rule cat2:
	input: in1=expand("{{id1}}-{{id2}}/fasta/{FASTA}.fa", FASTA=list_inFASTA),
			in2=expand( "{{id1}}-{{id2}}/fasta/{FASTA}.aln.fa", FASTA=list_inFASTA),
			in3=expand("{{id1}}-{{id2}}/identityCount/output_{FASTA}.identityCount", FASTA=list_inFASTA)
	output: "identityCount/{id1}-{id2}.identityCount"
	params: sge_opts="-l mfree=4G -l h_rt=48:00:00"
	priority: 6
	shell:
			" cat {input.in3} > {output}"

rule cat1:
	input: expand("summary_divTN93/{{id1}}_perFASTA/output_{FASTA}.divTN93", FASTA=list_inFASTA)
	output: "summary_divTN93/{id1}.TN93.txt"
	params: sge_opts="-l mfree=4G -l h_rt=48:00:00"
	priority: 7
	shell:
			" cat {input} > {output}"



rule summary_divTN93:
	input: expand("{{id1}}-{id2}/divTN93/output_{{FASTA}}.divTN93", id2=list2_ids)
	output: "summary_divTN93/{id1}_perFASTA/output_{FASTA}.divTN93"
	params: sge_opts="-l mfree=4G -l h_rt=48:00:00"
	priority: 8
	shell:
			""" cat {input} | bedtools merge -i - -c 4,5,6,7,8 -o distinct,mean,mean,mean,mean > {output} """

#			""" awk '{ total += $6; count++ } END { print $1,$2,$3,$4,$5,total/count }' > {output} """


rule divTN93:
	input:  "{id1}-{id2}/fasta/{FASTA}.fa",
			"{id1}-{id2}/fasta/{FASTA}.aln.fa",
			"{id1}-{id2}/identityCount/output_{FASTA}.identityCount"
	output: "{id1}-{id2}/divTN93/output_{FASTA}.divTN93"
	params: sge_opts="-l mfree=4G -l h_rt=48:00:00"
	priority: 9
	shell:
			""" numSeq=`grep \">\" {input[1]} | wc -l` || true ;  if [[ $numSeq -eq 2 ]] ; then basename {wildcards.FASTA} .fasta | awk -F'_' '{{print $1,$2,$3}}' | while read c1 c2 c3 ; do python scripts/computePairwiseDivTN93.py {input[1]} --list_addInfo \"['${{c1}}', '${{c2}}', '${{c3}}', '{wildcards.id1}-{wildcards.id2}']\" > {output} ; done ; else touch {output} ; fi """


rule identityCount:
	input: "%s/{FASTA}" % dir_inFASTA
	output: "{id1}-{id2}/fasta/{FASTA}.fa",
			"{id1}-{id2}/fasta/{FASTA}.aln.fa",
			"{id1}-{id2}/identityCount/output_{FASTA}.identityCount"
	params: sge_opts="-l mfree=4G -l h_rt=48:00:00"
	priority: 10
	shell:
			""" file=`basename {input} .fasta`; seqID=`grep \">\" {input} | grep -P "{wildcards.id1}|{wildcards.id2}"` ; numSeq=`echo ${{seqID}} | tr -cd '>' | wc -c` ; if [[ $numSeq -eq 2 ]]; then perl scripts/extractFASTAseqs.pl -n {dir_inFASTA}/{wildcards.FASTA} "${{seqID}}" > {output[0]} ; arr+=($(echo ${{seqID}} | sed 's/>//g' | sed 's/[:|-]/ /g')) ; mafft --maxiterate 100 --globalpair --adjustdirectionaccurately  {output[0]} > {output[1]} ; python scripts/calc_windowIdentity.v2.py --window 10000000000 --list_addInfo \"['${{arr[0]}}','${{arr[1]}}','${{arr[2]}}','${{arr[3]}}','${{arr[4]}}','${{arr[5]}}', '{wildcards.id1}-{wildcards.id2}']\" {output[1]} > {output[2]} ; else touch {output[0]} && touch {output[1]} && touch {output[2]} ;fi """
			


