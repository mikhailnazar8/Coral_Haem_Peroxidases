Finding Nested Domain Information
================

The domain information extracted in step 01 will be used here to find
out which domains overlap with the haem peroxidase domain and which
domains do not overlap using the subtract function in BEDtools

The first step is to extract the information of only the haem peroxidase
domains seperating it from information about the additional domains

``` bash
 grep "PF03098"  ${species}_ipr.tsv | cat >> "${species}_ipr_hp.tsv"
```

The second step is convert the tsv files into gff formnat using a custom
python script. This was done for both the entire ipr.tsv files and the
ipr_hp.tsv files

``` bash
python3 tsv2gff.py  ${species}_ipr.tsv gff_files/${species}.gff
```

The third step is to find out which domains are nested within the haem
peroxidase domain. This is done counter intuitively by finding out which
domains do not overlap with the haem peroxidase domain using BEDtools,
giving us a comprehensive list of all the non-overlapping domains.

``` bash
bedtools subtract -A -a ${species}.gff -b ${species}_hp.gff > ${species}_no_hp.gff
```

These results were then converted back into the same .tsv format as this
helped for easier downstream analysis

``` bash
python3 gff2tsv.py no_hp/${file}_no_hp.gff no_hp/${file}_ipr.tsv
```
