#!/usr/bin/env python3
"""
screen_sequences.py

Screen shotgun metagenomic contigs/reads (FASTA) against the custom ARG
k-mer index built by build_database.py, and write a raw hit table of
per-gene k-mer matches and coverage.

Usage:
    python scripts/screen_sequences.py --query examples/example_contigs.fasta \
        --index database/kmer_index.json --output hits.tsv
"""
import argparse
import json
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Screen sequences against the custom ARG k-mer index")
    parser.add_argument("--query", required=True, help="FASTA of shotgun contigs/reads to screen")
    parser.add_argument("--index", required=True, help="Path to the k-mer index JSON built by build_database.py")
    parser.add_argument("--output", required=True, help="Path to write the raw hit table (TSV)")
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


def screen(sequences, index):
    k = index["k"]
    kmer_to_genes = index["kmer_to_genes"]
    gene_lengths = index["gene_lengths"]

    hits = defaultdict(lambda: defaultdict(int))
    for query_id, seq in sequences.items():
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            for gene_id in kmer_to_genes.get(kmer, []):
                hits[query_id][gene_id] += 1

    rows = []
    for query_id, gene_hits in hits.items():
        for gene_id, matched_kmers in gene_hits.items():
            total_kmers = gene_lengths.get(gene_id, 0)
            coverage = matched_kmers / total_kmers if total_kmers else 0.0
            rows.append((query_id, gene_id, matched_kmers, total_kmers, round(min(coverage, 1.0), 4)))
    return rows


def main():
    args = parse_args()
    sequences = read_fasta(args.query)

    with open(args.index) as f:
        index = json.load(f)

    rows = screen(sequences, index)

    with open(args.output, "w") as f:
        f.write("query_id\tgene_id\tmatched_kmers\tgene_total_kmers\tcoverage_fraction\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"Wrote {len(rows)} raw hits to {args.output}")


if __name__ == "__main__":
    main()
