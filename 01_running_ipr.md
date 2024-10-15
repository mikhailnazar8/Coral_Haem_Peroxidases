Running interproscan
================
Mikhail Nazareth

To annotate the protein files from all the genomes (look at
Supplementary Table 1) with the domains present in each of the protein
sequences ̛

``` bash
 interproscan.sh -i proteins.fasta  --disable-precalc -f tsv -goterms -cpu 12
```
