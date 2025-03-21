# Cognate Model
Project for the evaluation of a model for cognate data implemented in RAxML-NG

## Requirements:
1. Setup and activate the conda environment: 
```
conda env create -f environment.yml python=3.10
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

3. Data
From lexibench, branch `for_cognate_data`

## Execution:
Data Properties:
```
python concept_language_ratio.py
python symbol_counts.py
python polymorphism_analysis.py
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
