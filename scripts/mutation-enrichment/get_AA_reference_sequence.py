import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from gtfparse import read_gtf


def get_dict_gene_coordinate_from_gtf(reference_gtf_path):
    # Load GTF
    reference_gtf = read_gtf(reference_gtf_path).to_pandas()
    # Get list of gene names
    list_gene_names = list(reference_gtf['gene_id'].unique())

    # Get dictionary of gene coordinates
    dict_gene_coordinates = {}
    for gene in list_gene_names:
        start_start_codon = int(reference_gtf[(reference_gtf['gene_id'] == gene) & (reference_gtf['feature'] == "start_codon")]['start'].iloc[0])
        end_start_codon = int(reference_gtf[(reference_gtf['gene_id'] == gene) & (reference_gtf['feature'] == "start_codon")]['end'].iloc[0])
        start_end_codon = int(reference_gtf[(reference_gtf['gene_id'] == gene) & (reference_gtf['feature'] == "stop_codon")]['start'].iloc[0])
        end_end_codon = int(reference_gtf[(reference_gtf['gene_id'] == gene) & (reference_gtf['feature'] == "stop_codon")]['end'].iloc[0])
        dict_gene_coordinates[gene] = {'start_codon': [start_start_codon, end_start_codon], 'stop_codon': [start_end_codon, end_end_codon]}

    return dict_gene_coordinates


def subset_sequence(seq, start_position, stop_position):
    # start_position and stop_position are >= 1
    nuc_string = str(seq.seq)[(start_position - 1):stop_position]
    return Seq(nuc_string)


def dict_AA_by_gene(reference_gtf_path, reference_alignment_path):
    # Read alignment fasta
    reference_seq = SeqIO.read(reference_alignment_path, "fasta")

    # Get gene coordinates
    dict_gene_coordinates = get_dict_gene_coordinate_from_gtf(reference_gtf_path)

    # Get dictionary with gene alignments
    dict_gene_fasta = {}
    for gene, gene_dict in dict_gene_coordinates.items():
        nuc_alignment = subset_sequence(reference_seq, gene_dict['start_codon'][0], gene_dict['stop_codon'][1])
        nuc_translate = nuc_alignment.translate()
        dict_gene_fasta[gene] = {'alignment': nuc_alignment, 'translate': nuc_translate}

    return dict_gene_fasta


def convert_dict_AA_by_gene_to_df(dict_AA):
    rows = []
    for gene, sequence in dict_AA.items():
        gene_AA_seq = str(sequence['translate'])
        for index, aa in enumerate(gene_AA_seq, start=1):
            rows.append({
                "gene": gene,
                "AA_pos": index,
                "AA": aa
            })
    return pd.DataFrame(rows)


if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Extract amino acid sequence from reference genome for each gene",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--reference_gtf', type=str, required=True, help='Reference GTF file with gene coordinates')
    parser.add_argument('--reference_fasta', type=str, required=True, help='Reference FASTA file')
    parser.add_argument('--output_df_aa_reference_sequence', type=str, required=True, help='Output TSV file with AA sequence of reference genome per gene')

    args = parser.parse_args()

    # Get AA sequence by gene from reference
    dict_AA_ref = dict_AA_by_gene(args.reference_gtf, args.reference_fasta)
    df_aa_reference_sequence = convert_dict_AA_by_gene_to_df(dict_AA_ref)

    # Write results to output file
    df_aa_reference_sequence.to_csv(args.output_df_aa_reference_sequence, sep="\t", index=False)
