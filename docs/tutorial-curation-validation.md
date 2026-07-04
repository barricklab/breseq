[Back to the Main Curation Tutorial Page](tutorial-curation.md)

In the curation cycle, you first run _breseq_. Then, you edit a GenomeDiff file that will be your list of mutations for further analysis by inspecting _breseq_ output. Finally, you need to check whether your edits are correct. This page covers the last step.

## Simulate the evolved/engineered genome

The `gdtools APPLY` command takes as input your initial reference genome file(s) and a GenomeDiff file describing mutations. It outputs a genome that has been mutated according to your GenomeDiff. If it has all the mutations, this should be the genome of your sequenced sample.

```
gdtools APPLY -f GENBANK -o hypothetical_genome.gbk -r input.gbk input.gd
```

You can pick several other output file formats, but GenBank or GFF3 are best since they preserve gene/feature annotations.

!!! note
    The evidence lines in the GenomeDiff file are NOT USED by `gdtools APPLY`. Why? Because something like a new junction or missing coverage doesn't fully describe how to exactly edit the nucleotide sequenvce of the genome. There is some ambiguity that it is up to you to fix, by manually editing the GenomeDiff file to add mutations based on this evidence.

## Re-running _breseq_ against the hypothetical genome

Now run _breseq_ again, but change both the output directory and the reference genome.

```
breseq -j 8 -l 80 -o hypothetical_genome_output -r hypothetical_genome.gbk read_file_1.fastq.gz read_file_1.fastq.gz
```

If you view the HTML output from this run, all of the mutations should now be missing!

If you still have unassigned evidence, your work is not done. See the next sections for some ways to more deeply explore the read alignments and advice regarding common cases you may encounter.

!!! note
    Be careful interpreting the output from a run against a hypothetical genome. All of the coordinates of genes and other features will be shifted if any deletions, insertions, or substitutions are being applied.

**Next:** [Exploring aligned reads](tutorial-curation-exploring-aligned-reads.md)