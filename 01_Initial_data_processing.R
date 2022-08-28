####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 1: Performs initial data processing and formatting starting with the fasta file of plasmid
#           sequences and plasmid similarity table.
# 
#
# The following outputs are used for downstream analysis: 
# full.net.dat.Rda, layer.id.table.Rda, node.id.table.Rda, plasmid.metadat.updated.Rda
# 
# Script tested for R version 4.1.1
####################################################################################################################


####################################################################################################################
# SCRIPT SET-UP
####################################################################################################################
# Set working directory to wherever your files are located

# Load necessary packages
library(tidyverse)
library(readr)

# Starting files (in Starting_file zipped folder):
# Plasmid sequence file:
all_single_circs.fa
# Requires unzipping

# Plasmid similarity file (plasmid-to-plasmid BLAST results):
Plasmidome_to_plasmidome_similarity.txt

# Contig coverage file (read-mapping):
plasmid_coverage.csv
####################################################################################################################


####################################################################################################################
# Section 1: EXTRACT NECESSARY DATA FROM SEQUENCES
####################################################################################################################

# Determine plasmid lengths based on sequence length

library(seqinr)
# Read the fasta file containing the plasmid sequences
plas.fast <- read.fasta("Starting_files/all_single_circs.fa")

# Count the length of each plasmid sequence
plas.lengths.vector <- getLength(plas.fast)
# Note, in case of errors: phylotools package can cause problems running the 
# getLength command in seqinr package
# Make sure phylotools is NOT loaded when using seqinr functions

library(phylotools)
# Extract the names of each plasmid sequence (First line starting with > before each sequence)
plas.names <- get.fasta.name("Starting_files/all_single_circs.fa")

# Bind the plasmid name vector and the plasmid lengths vector into a single dataframe
plas.lengths.bind <- as.data.frame(cbind(plas.names, plas.lengths.vector))

# Rename the columns in the plasmid name and length dataframe and make the plasmid length numeric
plasmid_lengths <- plas.lengths.bind %>%
  dplyr::rename(plasmid.name= plas.names,
         plasmid_length = plas.lengths.vector) %>%
  mutate(plasmid_length = as.numeric(plasmid_length))
####################################################################################################################


####################################################################################################################
# Section 2: FORMAT PLASMID SIMILARITY TABLE 
####################################################################################################################

# Read the data:
plasm_sim <- read.csv("Starting_files/Plasmidome_to_plasmidome_similarity.txt",sep="\t")

# Format the data, removing unnecessary columns and renaming columns where necessary:
plasm.sim1 <- plasm_sim %<>% 
  # Creates a boolean column in order to determine if subject plasmid is duplicated
  mutate(plasmid_subject_duplicated=ifelse(str_detect(sseqid,'Duplicated_'),T,F)) %>% 
  # Moves "plasmid_subject_duplicated" column between "plasmid_subject" and "length"
  select(plasmid_query=qseqid, plasmid_subject=sseqid, plasmid_subject_duplicated, everything()) %>% 
  # Removes the "Duplicated_" prefix in subject plasmids
  mutate(plasmid_subject=str_remove_all(plasmid_subject, 'Duplicated_')) 

# Join the plasmid_lengths dataframe to plasm.sim1 dataframe.
# Match the length of both queries and subjects
prep1 <- left_join(plasm_sim, plasmid_lengths, by=c("plasmid_query"="plasmid.name")) %>%
  dplyr::rename(plasmid.length.query=plasmid_length)

prep2 <- left_join(prep1, plasmid_lengths, by=c("plasmid_subject"="plasmid.name")) %>%
  dplyr::rename(plasmid.length.subject=plasmid_length)

# Identify the cow each plasmid came from, based on the first letter of the plasmid name (corresponding to cow)
prep3 <- prep2 %>%
  mutate(cow_query = str_sub(plasmid_query,1,1))%>%
  mutate(cow_subject = str_sub(plasmid_subject,1,1))%>%
  select(plasmid_query, cow_query, plasmid_subject, cow_subject, everything())
  
# Select relevant columns and rename plasmids/cows to nodes/layers in line with network vocabulary:
prep4 <- prep3 %>%
  select(-plasmid_subject_duplicated, -qlen, -slen) %>%
  rename(layer_from = cow_query, 
         node_from = plasmid_query,
         layer_to = cow_subject, 
         node_to = plasmid_subject)
# Contains 295,334 alignments

# In some cases because of "duplication" of sequences for alignments in order to BLAST all sequences against each other
# alignments between plasmid pairs may be repeated but reversed.
# To resolve this issue, group by the unique alignment and query start and end and select the first occurrence.
net.dat1 <- prep4 %>%
  group_by(node_from, layer_from, node_to, layer_to, qstart,qend) %>%
  slice_head() %>%
  ungroup()
# Number of alignments reduced to 172,257

# Assign an ID to each unique pair of aligned plasmid (each plasmid pair may
# have multiple alignments) and a ranking number to order each alignment
# within each plasmid pair.
# This will facilitate all downstream filtering and formatting
net.dat2 <- net.dat1 %>% 
  group_by(node_from, layer_from, node_to, layer_to) %>%
  dplyr::arrange(desc(length), .by_group=TRUE) %>%
  # Assign group id to each alignment pair:
  mutate(align.id = cur_group_id()) %>%
  # Set longest alignment to the reference alignment; identify nested and separate alignments
  mutate(rank.in.group = 1:n()) %>%
  select(align.id, node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, rank.in.group, everything())
####################################################################################################################


####################################################################################################################
# Section 3: ELIMINATE NESTED ALIGNMENTS
####################################################################################################################

# Identify and remove nested alignments by ordering all alignments for each unique plasmid pair 
# by their position and then comparing each alignment to the alignment in the next row.
# This procedure is repeated several times until no nested alignments remain.

# First compare all alignments for a plasmid pair with the longest alignment:
nest.elim1 <- net.dat2 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(desc(length), .by_group = TRUE) %>% 
  # Compare each alignment to the longest alignment of the pair
  mutate(nested=ifelse(qstart > first(qstart) & qend < first(qend) |
                         qstart == first(qstart) & qend < first(qend) |
                         qstart > first(qstart) & qend == first(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < first(qstart) | qstart > first(qend), T, F)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, nested, sep.align, rank.in.group, everything()) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  ungroup()
# Reduces number of alignments to 145,629


# Change order of alignments for each plasmid pair so that they are in order of appearance in sequence.
# Note: Rank will still indicate LONGEST alignment. 
nest.elim2 <- nest.elim1 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, nested, sep.align, rank.in.group, everything()) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  ungroup()
# Reduces number of alignments to 132,542

# Repeat 
nest.elim3 <- nest.elim2 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, nested, sep.align, rank.in.group, everything()) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  ungroup()
# Reduces number of alignments to 132,232

# Repeat
nest.elim4 <- nest.elim3 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, nested, sep.align, rank.in.group, everything()) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  ungroup()
# Reduces number of alignments to 132,203

# Repeat
nest.elim5 <- nest.elim4 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, nested, sep.align, rank.in.group, everything()) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  ungroup()
# Reduces number of alignments to 132,194

# Repeat
nest.elim6 <- nest.elim5 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, nested, sep.align, rank.in.group, everything()) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  ungroup()
# Reduces number of alignments to 132,192

# A final test shows the number of alignments has plateaued at 132,192 and all nested alignments have been removed
nest.elim7 <- nest.elim6 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, nested, sep.align, rank.in.group, everything()) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  select(-nested, -sep.align) %>%
  ungroup()
# Reduces number of alignments to 132,192
####################################################################################################################


####################################################################################################################
# Section 4: IDENTIFY IDENTICAL PLASMIDS IN DIFFERENT COWS
####################################################################################################################
# Identify plasmids that were isolated from different cows but are identical to each other. 
# In such cases, the plasmids will have different names. 
# This section identifies such plasmids and assigns them the same node id

# Begin with the dataframe created above that has no nested alignments:
head(nest.elim7)

# Find plasmids with identical lengths & calculate the percent of the plasmid length represented by the alignment
id.length <- nest.elim7 %>%
  filter_at(vars(plasmid.length.query, plasmid.length.subject), 
            all_vars(plasmid.length.query == plasmid.length.subject)) %>%
  mutate(align.length.percent = length/plasmid.length.query*100) %>%
  select(node_from, layer_from, node_to, layer_to, length, align.length.percent, everything())

# Note: There are some with alignments >100% of the length, which indicates a gap in the alignment


# Select plasmids with 100% similiarity over 100% alignment
id.plas100.100 <- id.length %>%
  filter(align.length.percent == 100 & pident == 100) 


# Aggregate this dataframe so that each node_from is followed by a list of all plasmids that are identical to it. 
id.plas.list <- id.plas100.100 %>%
  select(node_from, node_to) %>%
  group_by(node_from)%>%
  summarise(id = paste(node_to, collapse = ","))

# Separate the appended list of matching plasmids into separate columns.
# Doing so requires knowing the maximum number of items because this is the number of columns needed when disaggregating. 
# Use the max number of items / columns.
# Find max number of commas in the appended lists
max(sapply(strsplit((id.plas.list$id),','),length))
# Max appended list size size = 5

# Separate the appended plasmid names (separated by commas) into separate columns. 
# Fill in empty columns with NAs. 
# (Warning is normal because most sets of identical plasmids are 2)
id.plas.list1 <- id.plas.list %>%
  separate(id, into =c("id1","id2","id3","id4","id5"), sep = ",")

# Note: Each row is repeated by the number of matches but if there's a group of 2 plasmids, A and B, there's a row A, B, and then B,A
# Get rid of duplicates in different orders by:
#  1. Transposing the dataframe
#  2. Ordering each column.
#  3. Retransposing
#  4. Removing the duplicates

# 1. Transpose the data frame:
id.plas.list2 <- t(id.plas.list1)

# 2. Reorder each column:
id.plas.list_sorted <- apply(id.plas.list2,2, sort, na.last=TRUE, decreasing=F)

# Transpose back:
id.plas.list_sorted.transp <- t(id.plas.list_sorted)

# 3. Convert lists back to data frame:
id.plas.list3 <- as.data.frame(id.plas.list_sorted.transp)

# Filter distinct rows
id.plas.list4 <- id.plas.list3 %>%
  distinct(.)

View(id.plas.list4)
# 314 plasmids --> 138 UNIQUE plasmids found in 2-5 different cows
####################################################################################################################


####################################################################################################################
# Section 5: ASSIGN SINGLE ID NUMBERS TO IDENTICAL PLASMIDS
####################################################################################################################

# Transpose the data.frame so that each set of identical plasmids is its own column
id.plas.list5 <- as.data.frame(t(id.plas.list4))

# Gather all the columns into one column called "plasmid_name"
# The key is plasmid ID and will be filled with the column names
id.plas.list6 <- id.plas.list5 %>%
  gather(.,key=plasmid_id, value=plasmid_name)

# Final plasmid ID's will a number.
# Since all the columns start with V and then a number, remove the V to obtain a number for plasmid ID.
# Then remove all NA's (which are created due to differing numbers of identical plasmids in each "set")
id.plas.list7 <- id.plas.list6 %>%
  mutate(plasmid_id=gsub("V","",plasmid_id)) %>%
  drop_na()

# There are 138 sets of identical plasmids so plasmid ID's
# go from 1 - 138.
####################################################################################################################


####################################################################################################################
# Section 6: ASSIGN ID NUMBERS TO NON-IDENTICAL PLASMIDS
####################################################################################################################
# For the rest of the (non-identical) plasmids, start numbering at 139.

# But first remove the plasmids that are identical from the plasmid names list created 
# from the fasta file in Section 1:
plas.names.df <- as.data.frame(get.fasta.name("all_single_circs.fa"))
plas.names.df1 <- plas.names.df %>% 
  rename(plasmid_name=`get.fasta.name("all_single_circs.fa")`)

plas.no.ident <- anti_join(plas.names.df1, id.plas.list7, by="plasmid_name")
# The 314 plasmids that are identical to at least one other plasmid are successfully removed

# Create a sequence for the remaining plasmids, each getting a unique number from 139 - 8565 
pl.seq <- c(139:8565)

# Renumber the plasmid_id, start at 139:
plas.no.ident1 <- plas.no.ident %>%
  select(plasmid_name)%>%
  mutate(plasmid_id = pl.seq) %>%
  select(plasmid_id, plasmid_name)

# Bind the unique plasmids with their ID's to the rest:
plasmid.idnums <- rbind(id.plas.list7, plas.no.ident1)

# Add the cow origin and plasmid length to make a metadata file
plasmid_metadat.no.cowid <- plasmid.idnums %>%
  mutate(cow_name = paste(substr(plasmid_name,1,1))) %>%
  right_join(., plasmid_lengths, by=c("plasmid_name"="plasmid.name"))%>%
  drop_na()

# Assign ID number to cows:
cow.idnums <- plasmid_metadat.no.cowid %>%
  select(cow_name) %>%
  arrange(cow_name) %>%
  distinct() %>%
  mutate(cow_id=(c(1:length(cow_name))))

# Create final metadata file containing ID's for PLASMIDS and COWS
plasmid_metadat <- plasmid_metadat.no.cowid %>%
  right_join(., cow.idnums, by="cow_name") %>%
  mutate(node_id = as.integer(plasmid_id),
         layer_id = as.integer(cow_id)) %>%
  select(node_id, plasmid_name, layer_id, everything()) %>%
  select(-plasmid_id, -cow_id)

# Node id table (each node id occurs once)
node.id.table.prelim <- plasmid_metadat %>%
  select(node_id, plasmid_length) %>%
  distinct()
####################################################################################################################


####################################################################################################################
# Section 7: REPLACE NODE AND LAYER NAMES WITH ID NUMBERS CREATE EXTENDED EDGELIST
####################################################################################################################
# Use the plasmid-plasmid similarity table without duplicates or nested alignments:
head(nest.elim7)

net.dat3 <- nest.elim7 %>%
  left_join(., plasmid.idnums, by=c("node_from"="plasmid_name")) %>%
  select(-node_from) %>%
  rename(node_from=plasmid_id) %>%
  left_join(.,plasmid.idnums, by=c("node_to"="plasmid_name")) %>%
  select(-node_to) %>%
  rename(node_to=plasmid_id) %>%
  left_join(., cow.idnums, by=c("layer_from"="cow_name")) %>%
  select(-layer_from) %>%
  rename(layer_from=cow_id) %>%
  left_join(.,cow.idnums, by=c("layer_to"="cow_name")) %>%
  select(-layer_to) %>%
  rename(layer_to=cow_id) %>%
  select(layer_from,node_from, layer_to, node_to, everything())
####################################################################################################################


####################################################################################################################
# Section 8: ELIMINATE OVERLAPPING ALIGNMENTS
####################################################################################################################
# In some cases, alignments overlap. We attempted to extend the overlaps but this did not reflect the actual 
# alignments between plasmids.
# Therefore, for each plasmid pair, keep only the longest alignment when overlap occurs.

# Continue dataframe created above in Section 7
head(net.dat3)

# Use a similar procedure to what was done for eliminating nested alignments
no.overlap <- net.dat3 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(desc(length), .by_group = TRUE) %>% 
  # Identify overlapping alignments
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(sep.align=ifelse(qend < first(qstart) | qstart > first(qend), T, F)) %>%
  mutate(overlap=ifelse(qstart < first(qstart) & qend > first(qstart) &qend < first(qend) 
                        | qstart > first(qstart) & qstart < first(qend) 
                        & qend > first(qend) , T, F)) %>%
  mutate(overlap=ifelse(rank.in.group == min(rank.in.group), "TOP", overlap)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, sep.align, overlap, rank.in.group, everything()) %>%
  # Eliminate the overlapping segments
  filter(overlap == "TOP" | overlap == FALSE)
# Reduces alignments to 106,561

# Repeat
no.overlap2 <- no.overlap %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(desc(length), .by_group = TRUE) %>% 
  # Identify overlapping alignments
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  mutate(overlap=ifelse(qstart < lag(qstart) & qend > lag(qstart) &qend < lag(qend) 
                        | qstart > lag(qstart) & qstart < lag(qend) 
                        & qend > lag(qend) , T, F)) %>%
  mutate(overlap=ifelse(rank.in.group == min(rank.in.group), "TOP", overlap)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, sep.align, overlap, rank.in.group, everything()) %>%
  # Eliminate the overlapping segments
  filter(overlap == "TOP" | overlap == FALSE)
# Reduces alignments to 104,983


# Repeat
no.overlap3 <- no.overlap2 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(desc(length), .by_group = TRUE) %>% 
  # Identify overlapping alignments
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  mutate(overlap=ifelse(qstart < lag(qstart) & qend > lag(qstart) &qend < lag(qend) 
                        | qstart > lag(qstart) & qstart < lag(qend) 
                        & qend > lag(qend) , T, F)) %>%
  mutate(overlap=ifelse(rank.in.group == min(rank.in.group), "TOP", overlap)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, sep.align, overlap, rank.in.group, everything()) %>%
  # Eliminate the overlapping segments
  filter(overlap == "TOP" | overlap == FALSE)
# Reduces alignments to 104974

# Repeat
no.overlap4 <- no.overlap3 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(desc(length), .by_group = TRUE) %>% 
  # Identify overlapping alignments
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  mutate(overlap=ifelse(qstart < lag(qstart) & qend > lag(qstart) &qend < lag(qend) 
                        | qstart > lag(qstart) & qstart < lag(qend) 
                        & qend > lag(qend) , T, F)) %>%
  mutate(overlap=ifelse(rank.in.group == min(rank.in.group), "TOP", overlap)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, sep.align, overlap, rank.in.group, everything()) %>%
  # Eliminate the overlapping segments
  filter(overlap == "TOP" | overlap == FALSE)
# Reduces alignments to 104973

# Repeat
no.overlap5 <- no.overlap4 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(node_from, layer_from, node_to, layer_to) %>% 
  # Order each alignment pair by length
  arrange(desc(length), .by_group = TRUE) %>% 
  # Identify overlapping alignments
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  mutate(overlap=ifelse(qstart < lag(qstart) & qend > lag(qstart) &qend < lag(qend) 
                        | qstart > lag(qstart) & qstart < lag(qend) 
                        & qend > lag(qend) , T, F)) %>%
  mutate(overlap=ifelse(rank.in.group == min(rank.in.group), "TOP", overlap)) %>%
  select(align.id,node_from, layer_from, node_to, layer_to, length, pident,
         qstart, qend, sep.align, overlap, rank.in.group, everything()) %>%
  # Eliminate the overlapping segments
  filter(overlap == "TOP" | overlap == FALSE) %>%
  select(-sep.align, -overlap)
# Reduces alignments to 104973. Plateau reached
####################################################################################################################


####################################################################################################################
# Section 9: ADD IN THE 'MASKED' PLASMIDS
####################################################################################################################
# Due to the way contigs are constructed, there is a bias towards INTER-layer edges
# because highly similar contigs within a single cow are likely to be collapsed into a single contig. 
# Therefore, there are less highly similar plasmids (for example 99% similar with slightly different lengths) 
# in the same cow whereas there are many such cases between cows.

# In order to get around this issue, the reads from each cow have been mapped to the contigs to detect
# additional cows in which the contig occurs.

# Read the coverage file. Each row is a different plasmid, identified by the node_id (plasmid id number). 
# Each column is a cow. 
# The values represent coverage, the percent of each contig that is covered by at least one read. 
# 100% coverage indicates that the entire contig is recovered in the plasmid.
# NOTE: In a special case, if a plasmid is nested within another longer and completely 
# identical within the nested portion, this cannot distinguish between the
# nested plasmid and the longer plasmid. Both plasmids would be counted in such a case.
plas.cov <- read.csv("plasmid_coverage.csv",header=TRUE)

# Reorganize the coverage table:
plas.cov2 <- plas.cov %>%
  # Gather all the cows into a single column.
  # Column 1 = plasmid, Column 2 = cow, Column 3 = coverage in that cow
  gather(-Sequence.Sample, key="cow", value="coverage") %>%
  # Keep only plasmid-cow combinations with 100% (full) coverage
  filter(coverage==100) %>%
  # Order the dataframe by plasmid_id (for readability)
  arrange(Sequence.Sample) %>%
  # Change Sequence.Sample column name to node_id for consistency
  mutate(node_id = Sequence.Sample) %>%
  # Eliminate Sequence.Sample column
  select(-Sequence.Sample) %>%
  # Reorder columns starting with node_id
  select(node_id, everything()) %>%
  # Add in layer id's (cow id numbers) by joining the layer.id.table
  right_join(.,layer.id.table, by=c("cow"="cow_name")) %>%
  # Eliminate cow name (letter) and coverage (because all should be 100%) 
  select(-cow, -coverage) %>%
  # Keep only layer_id and node_id
  select(layer_id, node_id) %>%
  # Change layer_id and node_id to character for compatibility with 
  # other tables.
  mutate(layer_id = as.character(layer_id),
         node_id = as.character(node_id))
# Coverage data is now prepared to be combined with the rest of the alignment data

# Manipulate the alignment data (no.overlap5) for compatability:
cover.test1 <- no.overlap5 %>%
  select(-align.id,-rank.in.group) %>%
  mutate(layer_from = as.character(layer_from),
         layer_to = as.character(layer_to),
         node_from = as.character(node_from),
         node_to = as.character(node_to)) %>%
  ungroup()


# Now join the "masked" plasmids to the alignments. This will add the plasmids has now been found via read-mapping to those layers.
# First join the masked plasmids to the node_from
masked.node.from <- cover.test1 %>%
  # Join the plasmids from the masked file to the alignment table
  # This will add another column with a new layer that each plasmid was found in
  left_join(., plas.cov2, by=c("node_from"="node_id")) %>%
  mutate(mask.layer.from = layer_id) %>%
  select(-layer_id) %>%
  drop_na(.) %>%
  select(-layer_from) %>%
  mutate(layer_from = mask.layer.from) %>%
  select(-mask.layer.from) %>%
  select(node_from, layer_from, everything())

# Repeat for node_to
masked.node.to <- cover.test1 %>%
  left_join(., plas.cov2, by=c("node_to"="node_id")) %>%
  mutate(mask.layer.to = layer_id) %>%
  select(-layer_id) %>%
  drop_na(.) %>%
  select(-layer_to) %>%
  mutate(layer_to = mask.layer.to) %>%
  select(-mask.layer.to) %>%
  select(node_from, layer_from, node_to, layer_to, everything())

# Compare the masked plasmids to each other by adding to both from and to at same time
# E.g. this means adding alignments between masked plasmids when necessary.
masked.each.other <- cover.test1 %>%
  # Join the masked plasmids/layers to the node_from
  left_join(., plas.cov2, by=c("node_from"="node_id")) %>%
  # Rename the layer_id to identify as the layer_from
  mutate(mask.layer.from = layer_id) %>%
  select(-layer_id) %>%
  # Join the masked plasmids/layers to the node_to
  left_join(., plas.cov2, by=c("node_to"="node_id")) %>%
  # Rename the layer_id to identify the layer_to
  mutate(mask.layer.to = layer_id) %>%
  select(-layer_id) %>%
  drop_na(.) %>%
  select(node_from, layer_from, mask.layer.from, node_to, layer_to, 
         mask.layer.to,everything())

# Plasmids should only link to: different plasmids in the same cow OR the same plasmids in different cows.
# Because of the joining, there will also be cases where the same plasmid is linked to itself in the same cow, 
# e.g. layer_from = layer_to & node_from = node_to.
# Identify such cases, label them, then eliminate them:
masked.each.other2 <- masked.each.other %>%
  filter(!(mask.layer.from==mask.layer.to & node_from==node_to))%>%
  distinct()

# Relabel masked to just layer_from and layer_to
masked.each.other3 <- masked.each.other2 %>%
  mutate(layer_from = mask.layer.from) %>%
  mutate(layer_to = mask.layer.to) %>%
  select(-mask.layer.from, -mask.layer.to) %>%
  select(node_from, layer_from, node_to, layer_to, everything() )

# Bind all the masked dataframes:
masked.all.bind <- bind_rows(masked.node.from, masked.node.to, masked.each.other3)  

# Bind the masked plasmids
full.align.set <- bind_rows(cover.test1, masked.all.bind) %>%
  #Filter out any node/layer combos that are duplicated in the from/to (both layer and node) because this is a self-link
  # that should not be included 
  filter(!(layer_from==layer_to & node_from==node_to))%>%
  distinct()
# ~2500 alignments added for a total of 209,197
####################################################################################################################


####################################################################################################################
# Section 10: ELIMINATE ALIGNMENTS WITH GAPS
####################################################################################################################
# There are alignments with gaps, indicated by an alignment that is longer than the difference of qstart - qend
# Eliminating these alignments will also catch cases of gaps when an entire alignment is longer than either  
# plasmid in a pair leading to plasmids >100% similar, which is not possible and only due to gaps.
# This would also count a gapped alignment as a single alignment, which 
# doesn't make sense in the context of how edges will be measured.
full.align.set1 <- full.align.set %>%
  ungroup() %>%
  # What is the actual length between qstart and qend?
  mutate(qs.qe.diff = qend - qstart + 1) %>%
  # Keep only alignments when length = difference between qend, qstart
  # Note: there are no alignments with length < difference
  filter(length == qs.qe.diff)
# Total left: 104,094

# Separate data set of the gaps for further analysis / visualization:
gap.alignments <- full.align.set %>%
  ungroup() %>%
  # What is the actual length between qstart and qend?
  mutate(qs.qe.diff = qend - qstart + 1) %>%
  # Keep alignments when length > difference between qend, qstart
  # Note: there are no alignments with length < difference
  filter(length > qs.qe.diff) %>%
  mutate(gap.length = length - qs.qe.diff)
# 103,593

# Plots of gapped alignments:
gap.hist.all.plot <- ggplot(gap.alignments, aes(x=gap.length))+
  geom_histogram(color="black", fill="white") +
  labs(x = "Gap length", y="Count" )

gap.plot.short <- ggplot(gap.alignments, aes(x=gap.length))+
  geom_histogram(color="black", fill="white") +
  labs(x = "Gap length", y="Count" ) +
  scale_x_continuous(limits=c(0, 40)) + 
  labs(x = "Gap length", y="Count" )

gap.plot.long <- ggplot(gap.alignments, aes(x=gap.length))+
  geom_histogram(color="black", fill="white") +
  labs(x = "Gap length", y="Count" ) +
  scale_x_continuous(limits=c(40, 80)) + 
  labs(x = "Gap length", y="Count" )

# Plot of gaps vs. alignment length
gap.dot.plot <- ggplot(gap.alignments, aes(x=length, y=gap.length))+
  geom_point(alpha = 1/10, color="black", fill="white") +
  labs(x = "Alignment length (bp)", y="Gap length (bp)" )
####################################################################################################################


####################################################################################################################
# Section 11: FORMAT DATA WITH NO OVERLAPS INTO NETWORK WITH  INTRA/INTER LABELS AND EDGEWEIGHTS
####################################################################################################################
almost.full.net.dat <- full.align.set1 %>%
  # Coverting layer and node id's to integers as needed by
  # infomapecology
  mutate(layer_from = as.integer(layer_from),
         node_from = as.integer(node_from),
         layer_to = as.integer(layer_to),
         node_to = as.integer(node_to)) %>%
  # Reorder with layers before nodes, as needed by infomapecology:
  select(layer_from, node_from, layer_to, node_to, everything()) %>%
  mutate(edge_type=ifelse(layer_from==layer_to,'Intra','Inter'))
# 104,094 alignments
####################################################################################################################


####################################################################################################################
# Section 12: REMOVE REVERSED ALIGNMENTS
####################################################################################################################
# To avoid inflating the number of edges, eliminate any alignments that are duplicated but reversed:
# E.g. L1 N2 -> L3 N4 = L3 N4 -> L1 N2
rev.elim <- almost.full.net.dat %>% # From Section 11, above
  rowwise() %>%
  mutate(fr_grp = paste("L",layer_from, "N",node_from,sep = "")) %>%
  mutate(to_grp = paste("L",layer_to, "N",node_to, sep = "")) %>%
  mutate(uniq_pair = paste(sort(c(fr_grp, to_grp)), collapse = "_")) %>%
  select(-fr_grp, -to_grp) %>%
  select(uniq_pair, everything()) %>%
  group_by(uniq_pair, length, pident) %>%
  slice(1) %>%
  ungroup()
# Leaves 80,109 alignments

# Eliminate nested alignments when considering reversed alignments:
rev.elim2 <- rev.elim %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(uniq_pair) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  select(-nested, -sep.align) %>%
  ungroup()
# Leaves 74,253 alignments

#Rerun to determine if there are any further nested alignments:
rev.elim3 <- rev.elim2 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(uniq_pair) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  select(-nested, -sep.align) %>%
  ungroup()
# Leaves 74,243 alignments

rev.elim4 <- rev.elim3 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(uniq_pair) %>% 
  # Order each alignment pair by length
  arrange(qstart, desc(qend), .by_group = TRUE) %>% 
  # Identify nested using lag --> compares each row to row before it
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(nested=ifelse(qstart > lag(qstart) & qend < lag(qend) |
                         qstart == lag(qstart) & qend < lag(qend) |
                         qstart > lag(qstart) & qend == lag(qend), T, F)) %>%
  mutate(sep.align=ifelse(qend < lag(qstart) | qstart > lag(qend), T, F)) %>%
  # For plasmid pairs with only 1 alignment between them, nested will be NA
  # Therefore keep all rows if nested is false OR NA 
  filter(nested==FALSE | is.na(nested)) %>%
  select(-nested, -sep.align) %>%
  ungroup()
# Leaves 74,243 (Same as previous iteration, so no further runs necessary)

# In some cases, reversed alignments overlap. As explained above, for each plasmid pair, 
# keep only the longest alignment when overlap occurs.
rev.elim5<- rev.elim4 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(uniq_pair) %>% 
  # Order each alignment pair by length
  arrange(desc(length), .by_group = TRUE) %>% 
  mutate(rank.in.group = 1:n()) %>%
  # Identify overlapping alignments
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(sep.align=ifelse(qend < first(qstart) | qstart > first(qend), T, F)) %>%
  mutate(overlap=ifelse(qstart < first(qstart) & qend > first(qstart) &qend < first(qend) 
                        | qstart > first(qstart) & qstart < first(qend) 
                        & qend > first(qend) , T, F)) %>%
  mutate(overlap=ifelse(rank.in.group == min(rank.in.group), "TOP", overlap)) %>%
  # Eliminate the overlapping segments
  filter(overlap == "TOP" | overlap == FALSE) %>%
  select(-sep.align, -overlap)
# Leaves 68,604 alignments

# Repeat
rev.elim6 <- rev.elim5 %>%
  # Group by alignment (e.g. nodes/layers, from/to)
  group_by(uniq_pair) %>% 
  # Order each alignment pair by length
  arrange(desc(length), .by_group = TRUE) %>% 
  mutate(rank.in.group = 1:n()) %>%
  # Identify overlapping alignments
  # Sample with lag: mutate(time.difference = time - lag(time))
  mutate(sep.align=ifelse(qend < first(qstart) | qstart > first(qend), T, F)) %>%
  mutate(overlap=ifelse(qstart < first(qstart) & qend > first(qstart) &qend < first(qend) 
                        | qstart > first(qstart) & qstart < first(qend) 
                        & qend > first(qend) , T, F)) %>%
  mutate(overlap=ifelse(rank.in.group == min(rank.in.group), "TOP", overlap)) %>%
  # Eliminate the overlapping segments
  filter(overlap == "TOP" | overlap == FALSE) %>%
  select(-sep.align, -overlap, -rank.in.group)
# Leaves 68,604 (Stabilized, indicating all overlapping alignments have been removed)
####################################################################################################################


####################################################################################################################
# Section 13: Save the network data
####################################################################################################################
full.net.dat <- rev.elim6

# Save this as R data so it's easy to come back to:
save(full.net.dat, file="full.net.dat.Rda")
####################################################################################################################


####################################################################################################################
# Section 14: CREATE AND SAVE FINAL NODE ID, METADATA, STATE NODE TABLES WITH COMPLETE DATA
####################################################################################################################
# Create final state and physical node tables to account for the addition of the "masked" plasmids:

# Masked plasmids: 
# The masked plasmids are located in this dataframe, created in Section 9. 
# Includes the layers (cows) and nodes (plasmids) by ID's only (e.g., no names)
head(plas.cov2)

plas.cov3 <- plas.cov2 %>%
  mutate(node_id = as.numeric(node_id),
         layer_id = as.numeric(layer_id))


# Add in all the necessary metadata to the list of masked plasmids and the layers they were found in.
masked.metadat <- plas.cov3 %>%
  # Join plasmid length
  left_join(., node.id.table.prelim, by="node_id") %>%
  # Join cow names:
  left_join(., layer.id.table, by="layer_id") %>%
  # Join to plasmid_metadat to get plasmid names
  left_join(., plasmid_metadat, by = c("layer_id", "node_id", "plasmid_length", "cow_name")) %>%
  # Reorder according to plasmid_metadat table
  select(node_id, plasmid_name, layer_id, cow_name, plasmid_length)

# Bind these rows to the plasmid_metadat file:
plasmid.metadat.updated <- bind_rows(plasmid_metadat, masked.metadat) %>%
  # Now filter out any repeats:
  distinct() %>%
  # Because the plasmid name can sometimes we missing, group by the other variables:
  arrange(node_id)

# Save plasmid metadata file as Rdata:
save(plasmid.metadat.updated, file="plasmid.metadat.updated.Rda")

# State nodes (all information)
state.nodes <- plasmid.metadat.updated

# Node id table (each node id occurs once)
# Is also the physical node table
node.id.table <- plasmid.metadat.updated %>%
  select(node_id, plasmid_length) %>%
  distinct()

# Save as Rdata:
save(node.id.table, file="node.id.table.Rda")

# Layer id table
layer.id.table <- plasmid.metadat.updated  %>%
  select(layer_id, cow_name) %>%
  arrange(layer_id) %>%
  distinct()

# Save as Rdata:
save(layer.id.table, file="layer.id.table.Rda")
####################################################################################################################


####################################################################################################################
# Section 15: EDIT FASTA FILE
####################################################################################################################
# This section is for consistency in the fasta file but is not necessary for performing downstream analysis.

# Keep only one sequence from each group of identical plasmids

# Start with the identical plasmids from Section 5:
head(id.plas.list7)

# Filter to include only one distinct plasmid name for each plasmid id:
keep.in.fasta.id <- id.plas.list7 %>%
  distinct(plasmid_id, .keep_all = T)

# Make a list of the rest of the plasmids from each group that should be removed from the fasta
remove.fr.fasta <- anti_join(id.plas.list7, keep.in.fasta.id)

# Join keep.in.fasta.id to all the unique / non-identical plasmids, plas.no.ident1 object from Section 6:
keep.in.fasta <- rbind(keep.in.fasta.id, plas.no.ident1)

# Keep only plasmids in the keep.in.fasta list from the fasta file (keeps only one sequence per plasmid id number, 
# e.g. identical plasmid groups)
fasta.remove.dups <- plas.fast[c(which(names(plas.fast) %in% keep.in.fasta$plasmid_name))]

# Outwrite a new fasta file with the duplicates removed:
write.fasta(fasta.remove.dups,names = names(fasta.remove.dups), 
            file.out= "fasta.dups.removed.fasta")

# Rename the plasmids in the fasta file using the plasmid id numbers (from plasmid_cow_matrix_reconstruction script)

# First reformat the plasmid.idnums so the original name is the first
# column (as required by the rename.fasta command)
ref.table <- plasmid.idnums %>%
  select(plasmid_name, plasmid_id)%>%
  rename(old_name = plasmid_name,
         new_name=plasmid_id)

plas.fast.idnums <- rename.fasta(infile="fasta.dups.removed.fasta", ref_table=ref.table,
                                 outfile="plasmid.idnums.fasta")
####################################################################################################################