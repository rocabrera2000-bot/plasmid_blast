# plasmid_blast

Utility scripts for pairwise comparison of plasmid sequences using MUMmer.

## all_vs_all_nucmer.py

Run nucmer for every combination of FASTA files in a directory and collate the
results from show-coords into a single table.

### Example

```bash
python all_vs_all_nucmer.py plasmids/ -o results/ --nucmer-args "--maxmatch"
```

The command above creates `results/` with `.delta` and `.coords` files for each
pair and an `all_coords.tsv` summary of every alignment.
