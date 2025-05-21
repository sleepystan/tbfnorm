# tbfnorm
### Normalizing b-factors for analyzing dynamics of experimentally derived protein structures

## To run the script, we suggest to use Miniconda to install the environment via following command. Then activate the environment:
```
conda env create -f tbfnorm_env.yml -n tbfnorm
conda activate tbfnorm
```
To install Miniconda, please refer to [Anaconda: Installing Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install)
## Following commands can be used to run b-factor normalization with different inputs:
### Single PDB entry
```
python normalize_bfactors.py -i 6CAX -t pdb -o ./examples/outputs/single_pdb_6CAX_output/
```
### List of PDBs (txt file)
```
python normalize_bfactors.py -i ./examples/pdb_list_example.txt -t pdb_list -o ./examples/outputs/pdb_list_output/
```
### directory with files
```
python normalize_bfactors.py -i ./examples/dir_example/ -t file_dir -o ./examples/outputs/file_directory_output/
```
Example outputs can be found in the examples/example_outputs/ directory

## Credit

This repository is part of following work:
> [Placeholder](link)