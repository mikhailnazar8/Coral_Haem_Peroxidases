Running interproscan
================

To annotate the protein files from all the genomes (look at
Supplementary Table 1) with the domains present in each of the protein
sequences ̛

``` bash
 interproscan.sh -i ${species}_proteins.fasta  --disable-precalc -f tsv -goterms -cpu 12
```

after running the whole genomes(protein files) through interproscan, we
first create subset of all the proteins that have the Animal Haem
Peroxidase domain signature(PF03098) and their additional domain
information

``` bash
grep "PF03098" ${species}_interproscan.tsv | cat > ${species}_hp.tsv

cut -f1 ${species}_hp.tsv | sort | uniq > ${species}_prot_id.txt

sed -i 's/[[:space:]]//g' "${species}_prot_id.txt"

while read -r seq; do
  grep "$seq"  ${species}_interproscan.tsv | cat >> "${species}_ipr.tsv"
done < "${species}_prot_id.txt"
```

we also used a python script intergrated with SAMTOOLS to extract the
haem peroxidaase domains from the protein sequences (look at
ipr2fasta.py for more details)

``` bash
python3 ipr2fasta.py ${species}_interproscan.tsv ${species}_prot.fasta >${species}_hp_domains.fasta
```

All the domain files had their sequence IDs edited to provide
information on genome of origin for each protein sequence

``` bash
sed -i 's/^>/>A_cerv_/' acrop_cerv_hp_domains.fasta
sed -i 's/^>/>A_digi_/' acrop_digi_hp_domains.fasta
sed -i 's/^>/>A_mill_/' acrop_mill_hp_domains.fasta
sed -i 's/^>/>A_tumida_/' aethi_tumi_hp_domains.fasta
sed -i 's/^>/>C_virg_/' cras_virg_hp_domains.fasta
sed -i 's/^>/>N_vect_/' nema_vect_hp_domains.fasta
sed -i 's/^>/>D_pert_/' des_pert_hp_domains.fasta
sed -i 's/^>/>O_fave_/' orbi_fave_hp_domains.fasta
sed -i 's/^>/>P_dami_/' poci_dami_hp_domains.fasta
sed -i 's/^>/>P_ever_/' pori_ever_hp_domains.fasta
sed -i 's/^>/>S_pist_/' stylo_pist_hp_domains.fasta
sed -i 's/^>/>S_cili_/' syco_cili_hp_domains.fasta
sed -i 's/^>/>Xenia_sp_/' xenia_sp_hp_domains.fasta
```
