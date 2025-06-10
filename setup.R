# Install packages
here::i_am("setup.R")
library(here)
library(ape)
library(phytools)
library(picante)
library(tidyverse)

read_mapping_info <- function(map_fp = here("00_raw_data", "Supplementary_Tables.xlsx"), sheet = "Table 1") {
    if(grepl("xlsx", map_fp)) {
      d <- readxl::read_excel(map_fp, sheet = sheet, na = "",skip = 1) %>%
        rename(SampleID = 1)      
    } else {
      d <- read_table(map_fp) %>%
        rename(SampleID = 1)
    }
    return(d)
}


read_otu_tax_table <- function() {
  d <- read_tsv(here("00_raw_data", "feature-table_sort_de_final_sub_GP2_16000.txt"), skip = 1) %>%
    rename(OTU_ID = 1)
  
  # get taxa table with renamed asvs
  taxonomy_table_all <- d %>%
    select(OTU_ID, taxonomy) %>%
    separate(taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
             sep = "; ", fill = "right") %>%
    mutate(AltOTUID = paste0("OTU_", 1:n()))

  # Simplify OTU ID, but keep track of the old ids
  OTU_IDs_all <- taxonomy_table_all %>% select(OTU_ID, AltOTUID)
  
  taxonomy_table_all <- taxonomy_table_all %>% select(AltOTUID, Domain:Species) %>%
    rename(OTU_ID = AltOTUID)
  
  # Prepare OTU table by removing taxonomy and simplifying OTU_IDs
  otu_table_all <- d %>% select(-taxonomy) %>%
    mutate(OTU_ID = paste0("OTU_", 1:n())) # for convenience overwrite OTU_ID codes with simpler names
  
  # Construct otu lists
  input_list_all <- list(otu_table = otu_table_all, taxonomy_table = taxonomy_table_all,
                         OTU_IDs = OTU_IDs_all)
  
  input_list <- list(input_list_all = input_list_all)
  return(input_list)
}


read_phy_tree <- function(otu_table, OTU_IDs, tree_fp = here("00_raw_data", "phylogenetic_tree.tre"), 
                          visualize = FALSE) {
    require(ape)
    require(phytools)
    
    tree.raw <- phytools::read.newick(tree_fp)

    # check tree merge
    if(visualize) {
        writeLines("Plotting initial tree")
        plot(tree.raw)
    }
    
    # Filter tree for ASVs in dataset (the tree was built from all the OTUs in the
    # database, but we only want to keep those that are in our otu table)
    # The tree has fewer OTUs than are in the dataset, so later we will drop otus from 
    # the table that are not in the tree.
    db_tips <- tree.raw$tip.label[!(tree.raw$tip.label %in% OTU_IDs$OTU_ID)]
    
    tree.nodbtips <- ape::drop.tip(tree.raw, tip = db_tips)
    # check tree
    if(visualize) {
        writeLines("Plotting filtered tree")
        # plot resulting tree
        plot(tree.nodbtips)
    }
    
    # Rename tips to be simpler
    tiplab <- data.frame( tiplab = tree.nodbtips$tip.label) %>%
      left_join(OTU_IDs, by = c("tiplab" = "OTU_ID"))
    
    
    tree.simpletips <- tree.nodbtips
    
    tree.simpletips$tip.label <- tiplab$AltOTUID
    
    # reroot tree if necessary (doesn't appear necessary in this case. The archaea are not the root, but
    #there doesn't seem to be much long branch attraction. Archaea seem grouped together, and mitochondria and
    # cloroplasts appear to be filtered out))
    #View(taxonomy_table_all[taxonomy_table_all$OTU_ID %in% tree.simpletips$tip.label,])
    #tree.mprt <- phytools::midpoint.root(tree.simpletips)
    
    # reroot tree with arcahaea
    
    return(list(origroottree = tree.simpletips))
}

# Currently no cleaning to be done, this is a placeholder function in case I need a cleaning function late on
clean_sample_metadata <- function(sample_metadata = sample_metadata_all) {
    sample_metadata_mod <- sample_metadata %>%
      # add a second depth column with avg. depth in meters
      separate(`Depth (m)`, into = c("lowdepth", "highdepth"), sep = "-", remove = F) %>%
      mutate(lowdepth = as.numeric(lowdepth),
             highdepth = as.numeric(highdepth)) %>%
      mutate(AvgDepth = (lowdepth + highdepth)/2)
    
    return(sample_metadata_mod)
}


# Check for matching names between datasets
check_sample_metadata <- function(otu_table, sample_metadata, 
                                  sample_id_column = "SampleID",
                                  taxonomy = NULL,
                                  tree = NULL) {
    # Checking OTU table and sample metadata
    # Find the missing samples for otu_table and sample metadata
    missing_genomic_data <- na.omit(sample_metadata$SampleID[!(sample_metadata$SampleID %in% names(otu_table[-1]))]) # Samples missing from OTU table
    writeLines(paste("Samples missing from the OTU table that are present in the",
                     "metadata:", paste(missing_genomic_data, collapse = ", ")))
    
    missing_env_data <- names(otu_table[-1])[!(names(otu_table[-1]) %in% sample_metadata$SampleID )] # Samples missing from env data
    writeLines(paste("Samples missing from the metadata that are present in the",
                     "OTU table:", paste(missing_env_data, collapse = ", ")))
    writeLines(paste0("Now filtering out missing sample(s)..."))
    
    # Filter out the samples that don't occur in both
    otu_table <- otu_table %>% select(OTU_ID, !all_of(missing_env_data))
    
    sample_metadata <- sample_metadata %>% 
        filter((!!as.name(sample_id_column) %in% names(otu_table[-1])))
    
    # After filtering samples, remove any columns which are entirely NAs (or 0s for MAGs)
    cols_to_keep <- colSums(is.na(sample_metadata))<nrow(sample_metadata)
    writeLines(paste0("Removing columns with no data for any sample: ",
                      paste(names(which(cols_to_keep == FALSE)), collapse = ", ")))
    sample_metadata <- sample_metadata[,colSums(is.na(sample_metadata))<nrow(sample_metadata)]
    
    # Reorder otu table columns to match order in sample metadata
    otu_table <- otu_table[,c("OTU_ID", sample_metadata$SampleID)]

    writeLines(paste(nrow(sample_metadata), "samples in sample_metadata"))
    writeLines(paste(ncol(otu_table[-1]), "samples in otu_table"))
    
    # Filter OTUs with 0 abundance in otu table after samples were filtered
    otus_to_remove <- otu_table$OTU_ID[which(rowSums(otu_table[,-1]) == 0)]
    writeLines(paste0("Removing ",  length(otus_to_remove), " OTUs with zero abundnace for any sample (after filtering): ",
                      paste(otus_to_remove, collapse = ", ")))
    
    otu_table <- otu_table %>%
      filter(!(OTU_ID %in% otus_to_remove))

    # Checking taxonomy, tree, and otu table (if taxonomy and tree are present)
    # and reordering to match tree order
    if(!is.null(tree) & !is.null(taxonomy)) {
      
        # After filtering 
        # reorder otu table to match order of tips in tree
        match.phylo.otu <- picante::match.phylo.data(tree, 
                                                     otu_table %>% as_tibble() %>%
                                                         column_to_rownames(var = "OTU_ID"))
        tree <- match.phylo.otu$phy
        otu_table <- match.phylo.otu$data %>%
            rownames_to_column(var = "OTU_ID") %>%
            mutate(rowname = OTU_ID) %>%
            column_to_rownames()
        
        # Check that otu table genomes are all in the tree and taxonomy file
        OTU_IDs_missing_from_taxonomy <- otu_table$OTU_ID[!otu_table$OTU_ID %in% taxonomy$OTU_ID]
        writeLines(paste("OTU_IDs missing from the taxonomy that are present in the",
                         "OTU table:", paste(OTU_IDs_missing_from_taxonomy, collapse = ", ")))
        
        OTU_IDs_missing_from_tree <- otu_table$OTU_ID[!otu_table$OTU_ID %in% tree$tip.label]
        writeLines(paste("OTU_IDs missing from the tree that are present in the",
                         "OTU table:", paste(OTU_IDs_missing_from_tree, collapse = ", ")))
                                                   
        # reorder taxonomy to match otu_table and tree tip order and filter to match tips
        taxonomy <- taxonomy %>% 
          filter(!(OTU_ID %in% otus_to_remove)) %>%
          as_tibble() %>% column_to_rownames("OTU_ID")
        taxonomy <- taxonomy[rownames(match.phylo.otu$data),] %>%
            rownames_to_column(var = "OTU_ID") %>%
            mutate(rowname = OTU_ID) %>%
            column_to_rownames()
        
        # Report the number of OTUs to the user
        writeLines(paste0(nrow(otu_table), " taxa in otu table"))
        writeLines(paste0(length(tree$tip.label), " tips in tree"))
        writeLines(paste0(nrow(taxonomy), " taxa in taxonomy"))

    }

    
    # Create list of outputs to return:
    return_list <- list(sample_metadata = sample_metadata, otu_table = otu_table)
    
    if(!is.null(taxonomy)) {
        return_list$taxonomy <- taxonomy
    }
    
    if(!is.null(tree)) {
        return_list$tree <- tree
    }
    
    return(return_list)
}


transform_perc <- function(vec) {
  # Remove 0s for CLR transformation
  # See Smithson & Verkuilen 2006 (https://doi.org/10.1037/1082-989X.11.1.54)
  (vec * (length(vec) - 1) + 0.5) / length(vec)
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

otu_tax_input <- read_otu_tax_table()
# Read in both single and average
otu_tax_input_all <- otu_tax_input$input_list_all

sample_metadata_all <- read_mapping_info(map_fp = here("00_raw_data", "Supplementary_Tables.xlsx"),
                                          sheet = "Table 1") %>%
  clean_sample_metadata(sample_metadata = .) %>%
  select(SampleID, everything())


trees_all <- read_phy_tree(otu_table = otu_tax_input_all$otu_table, 
                           OTU_IDs = otu_tax_input_all$OTU_IDs,
                           tree_fp = here("00_raw_data", "phylogenetic_tree.tre"),
                           visualize = F)


# Check sample metadata and combine into one output
input_all <- check_sample_metadata(otu_table = otu_tax_input_all$otu_table,
                                    sample_metadata = sample_metadata_all,
                                    taxonomy = otu_tax_input_all$taxonomy_table,
                                    tree = trees_all$origroottree)

print("Done with setup.R")
# For convenience parse input into different formats:
#sample_metadata <- input$sample_metadata
#otu_table <- input$otu_table # NOTE: says this should be rarefied, but the columns sums don't appear that way
#taxonomy <- input$taxonomy
#filttree <- input$tree
#otu_translate <- otu_tax_input$OTU_IDs
