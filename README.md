# dozyg

"dozey-g": a sequence to [variation graph](https://pangenome.github.io/), [dozeu](https://github.com/ocxtal/dozeu)-based [yeeter](https://www.urbandictionary.com/define.php?term=Yeet)

## overview

`dozyg` is a sequence to graph mapper that exploits several properties common to most genome variation graphs to achieve runtime similar to best-of-class methods for read alignment to reference genomes.

1. Most pangenome variation graphs have a "manifold" linear property. That is, they may have large scale structural variation globally, but they are locally usually linear or partially orderable.
2. We can sort these graphs so that their colinear regions are represented contiguously within a sort order (e.g. using `odgi sort`).
3. This lets us use efficient collinear chaining methods to find target mappings, and apply POA locally (with `dozeu`) to obtain a base-level alignment.
4. A hierarchical chaining model lets us split the alignment across different regions of the graph, while maximizing our coverage of our mapped sequence.

## process

We first built a k-mer index over all k-mers in the graph (with pruning of complex regions to remove redundant k-mers) and a transformation of the graph designed to support efficient processing during read mapping.
The k-mers are hashed (this simplifies index construction and allows them to be of any length) and written into a [minimal perfect hash function](https://github.com/rizkg/BBHash).
We record their positions in the graph, and store linearized versions of the forward and reverse complement of the graph.

To align a sequence to the graph, we apply a two-stage clustering method.
The first stage is similar to that used in [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778), and chains with respect to the target sequence (our linearized graph).
Because our graph can contain structural rearrangements, we then add a second pass that combines the target-relative chains into "superchains" in a similar banded process that proceeds over the query sequence.

Finally, we align each sequence to the graph by progressing through the chains in each superchain and locally aligning them.
The final alignment is derived by applying [dozeu](https://github.com/ocxtal/dozeu) partial order alignment.
This process does not respect local complex structures like small inversions and cycles.
But, because of its two-stage chaining process, `dozyg` is able to align to graphs with complex large structural variation of all types.

## operation

`dozyg` reads [odgi](https://github.com/vgteam/odgi) graphs, indexes their kmers, and then maps reads from FASTA or FASTQ into a subset of the [GAF](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) graph alignment format.

```
odgi build -g g.gfa -o g.odgi
dozyg index -i g.odgi -k 15 -e 3 -t 4
dozyg map -i g.odgi -f reads.fq >aln.gaf
```

Downstream processing of GAF records in enabled by numerous algorithms in [vg](http://github.com/vgteam/vg), including `vg pack` and `vg call`.
Other methods working on GAF alignments can be applied, such as [gaffy](https://github.com/ekg/gaffy), which can project these into various graph matrix formats.

Future ergonomic improvements will allow the direct indexing of any pangenome graph in [GFAv1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format, and indexing only kmers that occur in actual genomes.

## scope

`dozyg` is designed to map sequences of any length, both small and large.
Different indexing patterns benift shorter versus longer reads.

## considerations

This is currently bleeding-edge research software.

`dozyg` doesn't map.
It yeets your sequence against the graph and hopes that it sticks.
It's designed to go fast, and not ask hard questions.
It will not get stuck in weird universal graph motifs.
Its hierarchical chaining model means that it is not afraid of complex pangenome graphs.

## author

Erik Garrison

## license

MIT
