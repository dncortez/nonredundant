nonredundant v1.0

This script creates a non redundant database of proteins from individual proteomes.
To run, place all the proteome files in a single folder and execute the following command

python nonredundant.py folder_name --check

For additional options see help.txt

This script creates a table that links the proteins from the separate files to the protein 
in the non redundant database. Its columns are as follows:

proteome	The filename of the proteome that the protein belongs to
n_protein	The position of the protein in said proteome, starting in 0
protein_id	The id of the protein in said proteome, defined by the fasta (faa) format
nr_id		The position of the protein in the non redundant database

The first 3 columns are related to the proteins in the original files. Their concatenation
is unique in the table. The last column can repeat, as multiple proteins can have the same
non-redundant id.

This script is a rudimentary attempt to optimize all vs all string comparison using python
datatypes and direct comparison, while keeping track of which protein ends up in which
position. Most likely, this can be highly improved.

Created by Diego Cortez diegonahuel8@gmail.com
2019