# Cognate Model
Project for the evaluation of a model for cognate data implemented in RAxML-NG

## Requirements:
1. Setup and activate the conda environment: 
```
conda env create -f environment.yml 
conda acitvate cognate-model
```
2. Adapted RAxML-NG binaries:
```
git clone --recursive https://github.com/luisevonderwiese/raxml-ng-cognate.git
cd raxml-ng-cognate
git checkout COG
mkdir build && cd build
cmake ..
make
cd ..
rm -r build/
git checkout COGs
mkdir build && cd build
cmake ..
make
```
Copy the binaries `raxml-ng-COG` and `raxml-ng-COGs` from `raxml-ng-cognate/bin` to `cognate-model/bin`

3. Generate Lexibench Data
Clone the [glottolog repo](https://github.com/glottolog/glottolog) to a directory of your choice, then run:
```
lexibench --repos data/lexibench download --upgrade
lexibench --repos data/lexibench lingpy_wordlists --plotstyle 2
lexibench --repos data/lexibench character_matrices --glottolog <your_glottolog_path> --formats bin.phy bin_part_2.phy bin_part_3.phy bin_part_4.phy bin_part_5.phy bin_part_6.phy bv_part_2.phy bv_part_3.phy bv_part_4.phy bv_part_5.phy bv_part_6.phy --plotstyle 2
lexibench --repos data/lexibench cross_validation_data
```

## Execution:
Data Properties:
```
python concept_language_ratio.py
python symbol_counts.py
python kappa_subset_sizes.py
python tabels.py
```
Tree Searches and Evaluation:
```
python inferences.py
```
Cross Validation:
```
python cross_validation.py
python bin_cross_validation.py
``` 
