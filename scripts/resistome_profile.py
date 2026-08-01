#!/usr/bin/env python3
"""
resistome_profile.py

Turn the raw hit table produced by screen_sequences.py into a resistome
abundance profile, aggregating hits by gene, gene family, and drug
class using the ARG metadata table.

Usage:
    python scripts/resistome_profile.py --hits hits.tsv \
        --metadata database/arg_metadata.tsv --output resistome_profile.tsv \
        --min-coverage 0.5
"""
import argparse
import csv
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate raw ARG hits into a resistome profile")
    parser.add_argument("--hits", required=True, help="Raw hit table TSV from screen_sequences.py")
    parser.add_argument("--metadata", required=True, help="ARG metadata TSV (database/arg_metadata.tsv)")
    parser.add_argument("--output", required=True, help="Path to write the resistome profile TSV")
    parser.add_argument("--min-coverage", type=float, default=0.5, help="Minimum coverage_fraction to count a hit (default: 0.5)")
    return parser.parse_args()


def load_metadata(path):
    metadata = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            metadata[row["gene_id"]] = row
    return metadata


def aggregate(hits_path, metadata, min_coverage):
    gene_counts = defaultdict(int)
    with open(hits_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if float(row["coverage_fraction"]) >= min_coverage:
                gene_counts[row["gene_id"]] += 1

    family_counts = defaultdict(int)
    drug_class_counts = defaultdict(int)
    for gene_id, count in gene_counts.items():
        meta = metadata.get(gene_id, {})
        family_counts[meta.get("gene_family", "Unknown")] += count
        drug_class_counts[meta.get("drug_class", "Unknown")] += count

    return gene_counts, family_counts, drug_class_counts


def main():
    args = parse_args()
    metadata = load_metadata(args.metadata)
    gene_counts, family_counts, drug_class_counts = aggregate(args.hits, metadata, args.min_coverage)

    with open(args.output, "w") as f:
        f.write("level\tname\tcount\n")
        for gene_id, count in sorted(gene_counts.items(), key=lambda x: -x[1]):
            gene_name = metadata.get(gene_id, {}).get("gene_name", gene_id)
            f.write(f"gene\t{gene_name}\t{count}\n")
        for family, count in sorted(family_counts.items(), key=lambda x: -x[1]):
            f.write(f"gene_family\t{family}\t{count}\n")
        for drug_class, count in sorted(drug_class_counts.items(), key=lambda x: -x[1]):
            f.write(f"drug_class\t{drug_class}\t{count}\n")

    print(f"Wrote resistome profile ({len(gene_counts)} genes) to {args.output}")


if __name__ == "__main__":
    main()
