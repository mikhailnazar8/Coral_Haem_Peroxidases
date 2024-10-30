Tree Data Frames and Metrics
================

To understand the clades in the tree a suite of functions was created to
help wrangle the data. First there was simply the function of creating a
tree data frame as was needed for the heatmap (seen in 04)

``` r
# Define the function
process_tree <- function(tree_file) {
  
  # Read the phylogenetic tree
  tree <- read.tree(tree_file)
  
  # Extract tree data for ggtree compatibility
  ggtree_data <- ggtree(tree)$data
  
  # Step to account for multiple domains in sequence IDs and add node information
  tree_tibble <- tree %>%
    as_tibble() %>%
    mutate(prot_accession = ifelse(grepl("_[1-3]$", label) & grepl("[A-Za-z]", label),
                                   sub("_[1-3]$", "", label),
                                   label)) %>%
    left_join(ggtree_data %>% select(node, isTip), by = "node") %>%
    mutate(node_number = ifelse(!isTip, node, NA))
  
  # Initialize play_tree
  play_tree <- tree_tibble
  
  # Get unique signal accessions
  signal_accessions <- unique(prot_table$sig_acc)
  
  # Add columns for each signal accession
 for (sig in signal_accessions) {
 
  temp_group <- prot_table %>% filter(sig_acc == sig)
  

  play_tree <- play_tree %>%
    mutate(!!sig := ifelse(prot_accession %in% temp_group$id, "yes", "no"))
  }
  

  play_tree <- as.data.frame(play_tree)
  
  tree_heatmap <- play_tree%>%
  mutate(label = make.unique(as.character(label))) %>%
  column_to_rownames(var = "label") %>% 
  select(!parent:node_number)
  
  return(tree_heatmap)
}

domain_df <- process_tree("data_files/run1.contree")
```

The next function needed is the function to find all the Sequence IDs in
a certain clade, for which you will need to know the node number of the
clade of interest

``` r
ftipnode <- function(tree, node) {
  descendants <- getDescendants(tree, node)
  tips <- descendants[descendants <= length(tree$tip.label)]
  return(tree$tip.label[tips])
}

tree <- read.tree("data_files/run1.contree") %>% 
  drop.tip(drop_all)
peroxinectin_clade <- ftipnode(tree, 426)
```

Using the list of Sequence IDs in a clade, you can then get various
metrics, such as domains within the clade and the Ratio of Presence (as
seen in Sec 2 of the paper)

``` r
find_domains <- function(clade, tree_df, filter_ratio) {
  # Total number of tips in the clade
  tot_num <- length(clade)
  
  # Subset the data frame to include only the rows (tips) in the clade
  clade_df <- tree_df[rownames(tree_df) %in% clade, ]
  
  # Convert "yes"/"no" to logical TRUE/FALSE
  clade_df_logical <- clade_df == "yes"
  
  # For each domain (column), count the number of tips with the domain
  domain_counts <- colSums(clade_df_logical)
  
  # Compute the ratio of tips with the domain to the total number of tips
  domain_ratios <- domain_counts / tot_num
  
  # Create a data frame to store the results
  result_df <- data.frame(
    sig_acc = names(domain_counts),
    ratio = domain_ratios
  )
  rownames(result_df) <- NULL
    result_df <- merge(result_df, prot_list, by.x = "sig_acc", by.y = "sig_acc", all.x = TRUE) %>% filter(ratio>filter_ratio) %>% 
      mutate(domain_overlap= ifelse(sig_acc %in% c(prot_list$sig_acc),F,T))
    
    tree_df_logical <- tree_df == "yes"
  total_domain_counts <- colSums(tree_df_logical) 
  result_df <- result_df %>%
    rowwise() %>%
    mutate(is_unique = domain_counts[sig_acc] == total_domain_counts[sig_acc]) %>%
    ungroup()
}

peroxinectin_clade_domain_information <- find_domains(peroxinectin_clade,domain_df,0)
```

This gives you a data frame as such

``` r
head(peroxinectin_clade_domain_information %>% arrange(ratio))
```

    ## # A tibble: 6 × 7
    ##   sig_acc            ratio sig_desc    ipr_acc ipr_desc domain_overlap is_unique
    ##   <chr>              <dbl> <chr>       <chr>   <chr>    <lgl>          <lgl>    
    ## 1 cd00037           0.0149 CLECT       -       -        FALSE          TRUE     
    ## 2 G3DSA:3.10.100.10 0.0149 -           IPR016… C-type … FALSE          TRUE     
    ## 3 PF00059           0.0149 Lectin C-t… IPR001… C-type … FALSE          TRUE     
    ## 4 PS50041           0.0149 C-type lec… IPR001… C-type … FALSE          TRUE     
    ## 5 SM00034           0.0149 CLECT_2     IPR001… C-type … FALSE          TRUE     
    ## 6 SSF56436          0.0149 C-type lec… IPR016… C-type … FALSE          TRUE
