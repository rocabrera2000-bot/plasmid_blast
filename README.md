# plasmid_blast

Utility scripts for pairwise comparison of plasmid sequences using MUMmer.

## all_vs_all_nucmer.py

Run nucmer for every combination of FASTA files in a directory and collate the
results from show-coords into a single table.

### Example

```bash
python all_vs_all_nucmer.py plasmids/ -o results/ --nucmer-args "--maxmatch" -t 4
```

The command above runs pairwise comparisons using four threads, creating
`results/` with `.delta` and `.coords` files for each pair and an
`all_coords.tsv` summary of every alignment.

## plasmid_cliques.py

Identify groups of plasmids that share alignment blocks of at least a specified
length and identity from the `all_coords.tsv` table produced by
`all_vs_all_nucmer.py`.  Each group represents a clique in which all plasmids
share the block with every other member.

### Example

```bash
python plasmid_cliques.py results/all_coords.tsv plasmids/ --min-length 5000 --plot
```

This command finds cliques connected by blocks of at least 5 kbp at ≥90% identity,
creates `group_05kb_*` directories containing member FASTA files, writes reports
summarising the groups and, with `--plot`, generates an HTML network graph.
