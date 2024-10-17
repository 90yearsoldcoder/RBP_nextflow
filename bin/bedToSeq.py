# python bedToSeq.py -i ../example/test.bed -p 5 -r /restricted/projectnb/benwol/jmh/script/ref/human/hg38_111/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
from pybedtools import BedTool
import argparse

import pybedtools

def get_sequence(chromosome, start, end, strand, reference):
    """
    Retrieve a DNA sequence from a specified region of a reference genome using pybedtools.
    
    Args:
    chromosome (str): The chromosome identifier (e.g., 'chr1', 'chr2', etc.).
    start (int): The start position in the chromosome (0-based index).
    end (int): The end position in the chromosome (exclusive).
    strand (str): Strand information ('+' for forward, '-' for reverse).
    reference (str): Path to the reference genome in FASTA format.
    
    Returns:
    str: The DNA sequence from the specified region, adjusted for strand.
    """
    # Create a BedTool object from a string defining the region
    region = f"{chromosome}\t{start}\t{end}\t.\t0\t{strand}"
    bedtool = pybedtools.BedTool(region, from_string=True)
    
    # Get the sequence associated with the region
    sequence = bedtool.sequence(fi=reference, s=True)
    
    # Read the sequence from the output file
    with open(sequence.seqfn) as f:
        next(f)  # skip header
        sequence = f.read().strip()
    
    return sequence


def read_bed_and_write_fasta(bed_file, fasta_file, reference, extend):
    with open(bed_file, 'r') as bed, open(fasta_file, 'w') as fasta:
        for line in bed:
            fields = line.strip().split()
            if len(fields) < 6:
                continue  # Ensuring there's enough fields for chr, start, end, name, score, strand
            
            chromosome, start, end, name, score, strand = fields[:6]
            print(chromosome, start, end, name, score, strand)
            sequence = get_sequence(chromosome, int(start) - int(extend), int(end) + int(extend), strand, reference)
            fasta.write(f"> {name}\n{sequence}\n")


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Convert BED to FASTA using specified genomic reference.")
    parser.add_argument("-i", "--input_bed", required=True, help="The input BED file.")
    parser.add_argument("-o", "--output_fasta", required=False, default="../processing/processing.fasta", help="The output FASTA file.")
    parser.add_argument("-r", "--reference", required=True, help="The reference genome file or identifier.")
    parser.add_argument("-p", "--append", required=False, default=0, help="The length extending seq")
    
    args = parser.parse_args()
    
    read_bed_and_write_fasta(args.input_bed, args.output_fasta, args.reference, args.append)