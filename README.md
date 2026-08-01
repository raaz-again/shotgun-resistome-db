# shotgun-resistome-db

Custom curated antibiotic resistance gene (ARG) database with tools to screen shotgun metagenomic sequences and build resistome abundance profiles.

## Overview

This repository provides a small, self-contained resistome analysis pipeline built around a custom ARG reference database. It screens assembled shotgun metagenomic contigs (or long reads) for antibiotic resistance genes using a pure-Python k-mer matching approach, then aggregates the resulting hits into a resistome profile broken down by gene, gene family, and drug class. No external aligner (BLAST, DIAMOND) is required.

## Database

The custom ARG database lives in the `database/` folder. `arg_database.fasta` contains reference sequences for each resistance gene, and `arg_metadata.tsv` describes each gene's family, drug class, and resistance mechanism. The bundled sequences are illustrative placeholders meant to demonstrate the pipeline end-to-end; swap them out for curated real reference sequences (e.g. from CARD, ResFinder, or MEGARes) before using this for real analysis. The metadata schema is: gene_id, gene_name, gene_family, drug_class, resistance_mechanism, notes.

## Installation

```
pip install -r requirements.txt
```

No external Python packages are required; every script uses only the standard library.

## Usage

First build the k-mer index from the ARG database.

```
python scripts/build_database.py --fasta database/arg_database.fasta --k 21 --output database/kmer_index.json
```

Then screen your shotgun metagenomic contigs against the index to get a raw hit table.

```
python scripts/screen_sequences.py --query examples/example_contigs.fasta --index database/kmer_index.json --output hits.tsv
```

Finally aggregate the raw hits into a resistome profile by gene, gene family, and drug class.

```
python scripts/resistome_profile.py --hits hits.tsv --metadata database/arg_metadata.tsv --output resistome_profile.tsv --min-coverage 0.5
```

## Input format

Query sequences should be assembled shotgun metagenomic contigs (or long reads) in FASTA format. See `examples/example_contigs.fasta`, which contains one exact match, one partial match, and one non-matching contig to demonstrate the different outcomes.

## Output format

`screen_sequences.py` writes a raw hit table with columns query_id, gene_id, matched_kmers, gene_total_kmers, and coverage_fraction. `resistome_profile.py` writes a long-format table with columns level (gene, gene_family, or drug_class), name, and count.

## License

Released under the MIT License. See [LICENSE](LICENSE).
