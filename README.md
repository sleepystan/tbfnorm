# tbfnorm
Script designed for normalizing b-factors of PDB structures for analyzing dynamics of experimentally derived protein structures

## To run the script, please install following conda environment:

## Following commands can be used to run b-factor normalization with different inputs:
### single PDB
python normalize_bfactors.py -i 6CAX -t pdb -o ./examples/outputs/single_pdb_6CAX_output/

### List of PDBs (txt file)
python normalize_bfactors.py -i ./examples/pdb_list_example.txt -t pdb_list -o ./examples/outputs/pdb_list_output/

### directory with files
python normalize_bfactors.py -i ./examples/dir_example/ -t file_dir -o ./examples/outputs/file_directory_output/

Example outputs can be found in the examples/example_outputs/ directory