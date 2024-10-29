Creating Gene tree
================

First all the species domain fasta files were compiled into a single
file to be aligned using mafft

``` bash
cat $species_hp_domains.fasta >> all_proteins.fasta

mafft --maxiterate 1000 --genafpair all_proteins.fasta > haem_perox.fasta
```

a locally installed version of goalign was then used to remove loci on
the alignment that had more than 50% gaps

``` bash
cat haem_perox.fasta| singularity run goalign.sif clean sites -c 0.5 > haem_perox_cleaned.fasta
```

the cleaned alignments were then used to make a phylogenetic tree using
IQ-TREE, using UltraFastBootstrapping as a support value at the nodes

``` bash
iqtree2 -s haem_perox_cleaned.fasta --prefix gene_tree -B 1000 -nt AUTO --cptime 400000
```
