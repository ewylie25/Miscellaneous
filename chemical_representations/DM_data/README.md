Data Matrix Metadata
--------------------
--------------------

Directory for storing meta data and helpful scripts related to the Distance Matrices.

#### legend.json

File contains hash map (dictionary) for matrices. The key is the atom hash value and the value is the corresponding row/column in matrix. This file corresponds to what I have labelled the "Ha_hybd_atom_aromat_1_Leg_hash_rings_1" version in Elasticsearch dm_final database.


#### generate_id_files.py  

Script to generate a set of text files formatted as input into `../DM_v_1_generation.py`. This is a helper script which generates the directory `source_files`. The idea is that multiple processes can be run (as in `../run.sh`) dividing the work load if the elasticsearch cluster is big enough and the hosting computer(s) has/have enough resources. 
