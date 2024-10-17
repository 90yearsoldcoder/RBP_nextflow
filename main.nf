#!/usr/bin/env nextflow

params.bed_files = "data/input_bed_files/*.bed"
params.extra_length = 5
params.reference_fasta = "data/reference/reference.fasta"
params.motif_library = "data/motif_library/known.rna.motifs"

process ConvertBedToSeq {
    input:
    file bed_file from file(params.bed_files)
    file ref_fasta from file(params.reference_fasta)

    output:
    file "sequences_${bed_file.simpleName}.fasta" into fasta_sequences

    """
    bedToSeq.py -i ${bed_file} -p ${params.extra_length} -r ${ref_fasta} -o sequences_${bed_file.simpleName}.fasta
    """
}

process FindRBPMotifs {
    input:
    file fasta_file from fasta_sequences
    file motif_lib from file(params.motif_library)

    output:
    file "motif_results_${fasta_file.simpleName}.txt" into motif_results

    """
    Motif.py -i ${fasta_file} -mr ${motif_lib} -o motif_results_${fasta_file.simpleName}.txt
    """
}

workflow {
    ConvertBedToSeq()
    FindRBPMotifs()
}
