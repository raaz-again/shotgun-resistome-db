#!/usr/bin/env python3
"""
build_database.py

Build a k-mer index for the custom antibiotic resistance gene (ARG)
database (database/arg_database.fasta + database/arg_metadata.tsv).

The index maps each k-mer to the set of ARG gene IDs it belongs to,
which scripts/screen_sequences.py later uses to screen shotgun
metagenomic contigs/reads against the database without needing an
external aligner such as BLAST or DIAMOND.

Usage:
    python scripts/build_database.py --fasta database/arg_database.fasta \
        --k 21 --output database/kmer_index.json
"""
import argparse
import json
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Build a k-mer index for the custom ARG database")
    parser.add_argument("--fasta", required=True, help="Path to the ARG reference FASTA file")
    parser.add_argument("--k", type=int, default=21, help="K-mer size (default: 21)")
    parser.add_argument("--output", required=True, help="Path to write the k-mer index (JSON)")
    return parser.parse_args()


def read_fasta(path):
    sequences = {}
    header = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if header is not None:
            sequences[header] = "".join(chunks)
    return sequences


def build_index(sequences, k):
    kmer_to_genes = defaultdict(set)
    gene_lengths = {}
    for gene_id, seq in sequences.items():
        gene_lengths[gene_id] = max(len(seq) - k + 1, 0)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            kmer_to_genes[kmer].add(gene_id)
    return kmer_to_genes, gene_lengths


def main():
    args = parse_args()
    sequences = read_fasta(args.fasta)
    kmer_to_genes, gene_lengths = build_index(sequences, args.k)

    index = {
        "k": args.k,
        "gene_lengths": gene_lengths,
        "kmer_to_genes": {kmer: sorted(genes) for kmer, genes in kmer_to_genes.items()},
    }

    with open(args.output, "w") as f:
        json.dump(index, f)

    print(f"Indexed {len(sequences)} ARG sequences into {len(kmer_to_genes)} unique {args.k}-mers -> {args.output}")


if __name__ == "__main__":
    main()
