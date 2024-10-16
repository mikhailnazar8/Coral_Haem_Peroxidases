#!/usr/bin/env python3

import argparse
import sys
import math
import itertools
import re
import logging
import subprocess
from operator import itemgetter


log = logging.getLogger()
logging.basicConfig(stream=sys.stderr,level=logging.DEBUG)

parser = argparse.ArgumentParser(description='Extract Interpro domain sequences and export to fasta')

parser.add_argument('iprin',metavar='<TSV FILE>',nargs=1,type=argparse.FileType('r'),help='Interproscan tsv file')

parser.add_argument('fasta',metavar='<FASTA FILE>',nargs=1,type=str,help='Protein fasta file')

parser.add_argument('domain',type=str,nargs='?',help='Domain Accession',default="PF03098")

args  = parser.parse_args()

# Fields in the interproscan output
#id	md5	length	analysis	sig_acc	sig_desc	start	end	score	status	 date	ipr_acc	ipr_desc	goterm

domain_sequence_dict={}

for line in args.iprin[0]:

	line_values = line.strip().split("\t")

	prot_id=line_values[0]
	sig_acc=line_values[4].strip()

	if sig_acc==args.domain:


		start=int(line_values[6])-1
		end=int(line_values[7])-1


		samtools_args=["samtools","faidx",args.fasta[0],prot_id]

		samtools_result = subprocess.Popen(samtools_args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

		fasta_records={}

		
		id_line=samtools_result.stdout.readline().decode('ascii').strip()


		line_parts = [line.decode('ascii').strip() for line in iter(samtools_result.stdout.readline, b"")]
		sequence = ''.join(line_parts)
		domain_sequence = sequence[start:end]


		current_seq_list = domain_sequence_dict.get(prot_id,[])
		current_seq_list.append(domain_sequence)

		domain_sequence_dict[prot_id] = current_seq_list


# This line invokes the debugger.  If you want to set a breakpoint for debugging uncomment this line and put it at the point in code where you want to debug.
#import pdb;pdb.set_trace()

for key,values in domain_sequence_dict.items():
	vc=0
	for v in values:
		if vc==0:
			v_id = key
		else:
			v_id = key + "_" + str(vc)
		vc+=1
		sys.stdout.write(">"+v_id+"\n"+v+"\n")

