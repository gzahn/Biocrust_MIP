#!/bin/bash

# this will download and build the custom databases
# The Eukaryome database, supplemented with maarjaam VXT assignments for glomeromycota taxa

# Run this from the directory where you will keep your databases for your project...
# depends on:
# wget
# R v 4.0+
# tidyverse R package suite
# ShortRead R package (will install automatically to .libPaths() )
# mycobank R package (will install automatically to .libPaths() )

# to run this on a cluster, you must have that software installed and loaded (via modules, for example)
# to run, make sure the format_maarjAM.R script is in the same directory as this script,
# then use ./build_taxonomy_database.sh
# and things SHOULD work fine...



# The script:



# get Eukaryome databases
wget https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_SSU_v1.8.zip
unzip General_EUK_SSU_v1.8.zip && rm General_EUK_SSU_v1.8.zip
wget https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v1.8.zip
unzip General_EUK_ITS_v1.8.zip && rm General_EUK_ITS_v1.8.zip


# reformat Eukaryome to play nicely with dada2 assignTaxonomy() function
cat General_EUK_SSU_v1.8.fasta | sed 's/>\([^;]*;\)/>/' | sed 's/>\([^|]*|\)/>/' | sed 's/|/;/g' > Eukaryome_General_SSU_v1.8_reformatted.fasta
cat General_EUK_ITS_v1.8.fasta | sed 's/>\([^;]*;\)/>/' | sed 's/>\([^|]*|\)/>/' | sed 's/|/;/g' > Eukaryome_General_ITS_v1.8_reformatted.fasta


# get maarjAM databases (fasta)
wget https://maarjam.ut.ee/resources/maarjam_database_SSU_TYPE.fasta.2021.zip -O maarjam_database_SSU_TYPE_2021.fasta.gz
wget https://maarjam.ut.ee/resources/maarjam_database_SSU.fasta.2021.zip -O maarjam_database_SSU_2021.fasta.gz
wget https://maarjam.ut.ee/resources/maarjam_database_full_ITS.fasta.2021.zip -O maarjam_database_full_ITS_2021.fasta.gz

# decompress
gunzip maarjam_database_*gz

# Run custom R script to add taxonomy to maarjam databases
# This uses R to work on those files and get the headers into shape
# It's just much nicer than trying to do it all in bash.
Rscript ./format_maarjAM.R

# The R script built some newly formatted files...
# combine the databases
cat maarjam_database_SSU*newheaders.fasta Eukaryome_General_SSU_v1.8_reformatted.fasta > Eukaryome_General_SSU_v1.8_reformatted_VTX.fasta
cat maarjam_database_full_ITS_2021_newheaders.fasta Eukaryome_General_ITS_v1.8_reformatted.fasta > Eukaryome_General_ITS_v1.8_reformatted_maarjam.fasta

# compress the final files
gzip *.fasta

