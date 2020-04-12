# gyeet

"gee-yeet": the sequence to [variation graph](https://pangenome.github.io/) [yeeter](https://www.urbandictionary.com/define.php?term=Yeet)

## overview

`gyeet` is a sequence to graph mapper that exploits several properties common to most genome variation graphs to achieve runtime similar to best-of-class methods for read alignment to reference genomes.

1. Most genome variation graphs have a "manifold" linear property. That is, they may have large scale structural variation globally, but they are locally usually very linear.
2. We can sort these graphs so that their colinear regions are represented contiguously within a sort order.
3. By projecting the graph into a string, but remembering its topolgy, we can use fast alignment algorithms for strings to achieve graphical alignments. We just have to fix up the alignment after the fact.

## process

We first built a k-mer index over all k-mers in the graph (with pruning of complex regions to remove redundant k-mers) and a transformation of the graph designed to support efficient processing during read mapping.
The k-mers are hashed (this simplifies index construction and allows them to be of any length) and written into a [minimal perfect hash function](https://github.com/rizkg/BBHash).
We record their positions in the graph.

To align a sequence to the graph, we apply a two-stage clustering method.
The first stage is similar to that used in [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778), and chains with respect to the target sequence (our linearized graph).
Because our graph can contain structural rearrangements, we then add a second pass that combines the target-relative chains into "superchains" in a similar banded process that proceeds over the query sequence.

Finally, we align each sequence to the graph by progressing through the chains in each superchain and locally aligning them.
We apply a high-performance pairwise alignment algorithm (currently [edlib](https://github.com/Martinsos/edlib)) to map subsequences of our query to subsequences of the graph string, and then corrects the alignment to describe a path through the actual graph.
Aligning to a string greatly simplifies our process relative to methods like those in `vg` (`vg map`, `vg mpmap`, `vg giraffe`).
However, it also introduces issues, with invalid graph path being possible.

Current development is focused on integrating [dozeu](https://github.com/ocxtal/dozeu) to align through complex regions of the graph, and only using edlib for alignment to long contigs (or nodes), which will result in a hybrid local alignment approach.
Once integrated, the model will fully respect local graph structures except for small cyclic and inverting components.

## operation

`gyeet` reads [odgi](https://github.com/vgteam/odgi) format graphs, indexes them, and maps reads from FASTA or FASTQ into a subset of the [GAF](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) graph alignment format.

```
odgi build -g g.gfa -s -o g.odgi
gyeet index -i g.odgi -k 15 -e 3 -t 4
gyeet map -i g.odgi -f reads.fq >aln.gaf
```

Downstream processing of GAF records in enabled by [gaffy](https://github.com/ekg/gaffy), which can project these into various graph matrix formats.

## considerations

This is currently bleeding-edge research software.

## author

Erik Garrison

## license

MIT
