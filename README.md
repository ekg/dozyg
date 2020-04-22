# edyeet

the [edlib](https://martinsos.github.io/edlib/)-based sequence to sequence [yeeter](https://www.urbandictionary.com/define.php?term=Yeet)

## overview

`edyeet` is a sequence to sequence mapper based on the edlib alignment algorithm, minimal perfect k-mer hashing, and two-stage alignment chaining

## process

We first built a k-mer index over all k-mers in the input sequence set.
The k-mers are hashed (this simplifies index construction and allows them to be of any length) and written into a [minimal perfect hash function](https://github.com/rizkg/BBHash).
We record their positions in the input sequences.

To align, we apply a two-stage clustering method.
The first stage is similar to that used in [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778), and chains with respect to the target sequence (our linearized graph).
To account for structural rearrangements, we then add a second pass that combines the target-relative chains into "superchains" in a similar banded process that proceeds over the query sequence.

Finally, we align each sequence to the target by progressing through the chains in each superchain and locally aligning them with edlib.
Each chain alignment is emitted in PAF format.

## operation

```
edyeet -i target.fa -f query.fq >aln.paf
```

## considerations

For an established, full-featured tool, use minimap2.
If you need to align things really, really fast, use edyeet.

## author

Erik Garrison

## license

MIT
