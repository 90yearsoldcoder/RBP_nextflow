#python Motif.py -mr /restricted/projectnb/casa/bu_brain_rnaseq/hippo_rawdata/circrna/homer/data/knownTFs/known.rna.motifs
from rk import motifScore
import argparse
import json

def local_align(input_seq, RBPname, RBPhelper7, RBPhelper8,  cutoff_score, detail:bool=False):
    RBPhelper7.initSeq()
    RBPhelper8.initSeq()
    score_array = []
    motif_count = 0
    position_array = []
    ind = 0

    for base in input_seq:
        ind += 1
        score7 = RBPhelper7.askScore(base, 7, RBPname)
        score8 = RBPhelper8.askScore(base, 8, RBPname)
        if score7 and score7 >= cutoff_score:
            score_array.append(score7)
            position_array.append((ind - 6, ind))
            motif_count += 1
        if score8 and score8 >= cutoff_score:
            score_array.append(score8)
            position_array.append((ind - 7, ind))
            motif_count += 1

    if detail == True:
        return motif_count, score_array, position_array
    return motif_count, score_array

def read_fasta(filepath):
    """
    Reads a FASTA file and returns a dictionary mapping gene names to sequences.

    Args:
    filepath (str): The path to the FASTA file.

    Returns:
    dict: A dictionary where keys are gene names and values are sequences.
    """
    with open(filepath, 'r') as file:
        fasta_dict = {}
        gene_name = ''
        sequence = ''

        for line in file:
            line = line.strip()
            if line.startswith('>'):  # New gene/sequence
                if gene_name:  # Save the previous gene and sequence
                    fasta_dict[gene_name] = sequence
                gene_name = line[1:].strip()  # Remove the '>' and take the rest as the gene name
                sequence = ''  # Reset sequence for the new gene
            else:
                sequence += line  # Continue reading the sequence

        # Add the last gene and sequence to the dictionary
        if gene_name:
            fasta_dict[gene_name] = sequence

    return fasta_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get the possible RBP motifs on each sequence")
    parser.add_argument("-i", "--input_fasta", required=False, default="../processing/processing.fasta", help="The input fasta file")
    parser.add_argument("-o", "--output_json", required=False, default="../result/result.json", help="The output json file.")
    parser.add_argument("-mr", "--motif_reference", required=True, help="The reference file for motifs")
    parser.add_argument("-cf", "--cutoff_score", required=False, default=5, help="The cutoff score for selecting RBP motifs")
    args = parser.parse_args()

    #matrixPath = "/restricted/projectnb/casa/bu_brain_rnaseq/hippo_rawdata/circrna/homer/data/knownTFs/known.rna.motifs"
    JSONpath = "./RBPdic.json"
    RBPhelper7 = motifScore(7, args.motif_reference, JSONpath)
    RBPhelper8 = motifScore(8, args.motif_reference, JSONpath)
    RBPlist = RBPhelper7.getRBPlist()
    
    RNA_dic = read_fasta(args.input_fasta)

    record = {}
    for RNA in RNA_dic:
        record[RNA] = {}
        record[RNA]['type'] = 'circular' if 'circ' in RNA else 'linear'
        record[RNA]['seq'] = RNA_dic[RNA]
        record[RNA]['motifs'] = []
        for RBP_ind in range(len(RBPlist)):
            RBP = RBPlist[RBP_ind]
            print("Working on circRNA:{RNA} and RBP:{RBP}".format(RBP=RBP, RNA=RNA))
            motif_count, score_array, position_array = local_align(RNA_dic[RNA], RBP, RBPhelper7, RBPhelper8, args.cutoff_score, detail=True)
            if motif_count == 0:
                continue
            record[RNA]['motifs'].append(
                {
                    'motif_name': RBP,
                    'RBP_key': RBP_ind,
                    'count': motif_count,
                    'scores': score_array,
                    'positions': position_array,
                }
            )
        record[RNA]['motifs'] = sorted(record[RNA]['motifs'], key=lambda x: x['count'], reverse=True)


    with open(args.output_json, 'w') as f:
        json.dump(record, f, indent=4)