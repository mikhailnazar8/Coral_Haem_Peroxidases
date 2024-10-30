import argparse
import logging
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import warnings

# Suppress specific Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)

# Set up logging
log = logging.getLogger()
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

# Argument parser setup
parser = argparse.ArgumentParser(description='Translate DNA or RNA sequences to protein sequences and extract ORFs.')
parser.add_argument('fasta', metavar='<FASTA FILE>', type=str, help='DNA/RNA fasta file')
parser.add_argument('-f', '--frame', type=int, choices=[0, 1, 2, 3, -1, -2, -3], help='Translation Frame.', default=1)
parser.add_argument('-s', '--split', action="store_true", help='Split into separate sequences (ORFs) based on stop codons')
parser.add_argument('-lo', '--longest_orf', action="store_true", help="Find the longest Open Reading Frame for each protein sequence")

# Parse arguments
args = parser.parse_args()

# Function to generate ORFs for the entire translated sequence
def get_orfs(translated_seq, parent_accession):
    start = 0
    N = 0
    for orf in translated_seq.split("*"):
        if orf:
            os = start
            oe = start + len(orf) * 3
            start = oe + 3

            N += 1
            yield {
                'parent': parent_accession,
                'start': os + 1,
                'end': oe,
                'orfnum': N,
                'length': len(orf),
                'seq': orf
            }

try:
    with open(args.fasta, "r") as file:
        for sequence in SeqIO.parse(file, "fasta"):
            longest_orf = None
            frames = [args.frame] if args.frame != 0 else [1, 2, 3, -1, -2, -3]

            for frame in frames:
                seq = sequence.seq
                if frame < 0:
                    seq = seq.reverse_complement()
                seq = seq[abs(frame) - 1:]
                protein_seq = seq.translate()

                if args.longest_orf:
                    for orf in get_orfs(str(protein_seq), sequence.id):
                        if longest_orf is None or orf['length'] > longest_orf['length']:
                            longest_orf = orf
                elif args.split:
                    for orf in get_orfs(str(protein_seq), sequence.id):
                        orf_id = f"{orf['parent']}_ORF{orf['orfnum']}_Frame{frame}"
                        print(f">{orf_id} [start:end = {int((orf['start'] + 2) / 3)}:{int((orf['end'] + 3) / 3)} | length = {orf['length']} | frame = {frame}]\n{orf['seq']}\n")
                else:
                    print(f">{sequence.id}_Frame{frame}\n{protein_seq}\n")

            if args.longest_orf and longest_orf:
                orf_id = f"{longest_orf['parent']}_ORF{longest_orf['orfnum']}"
                print(f">{orf_id} [start:end = {int((longest_orf['start'] + 2) / 3)}:{int((longest_orf['end'] + 3) / 3)} | length = {longest_orf['length']} | frame = {frame}]\n{longest_orf['seq']}\n")

except Exception as e:
    logging.error(f"An error occurred: {e}")
