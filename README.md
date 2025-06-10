# Background
Scripts for the assembly modeling analysis for the paleo-microbime project.

# Setup

## Create conda environment

``` bash
git clone git@github.com:hhollandmoritz/PaleoMicrobiomeAssemblyAnalysis.git
cd PaleoMicrobiomeAssemblyAnalysis
./install_dependencies.sh
conda activate PaleoMicrobiome_vX
```

If you require packages not installed in the above environment, add them to the .yml file and update the version number.

## Create your analysis file

Source `setup.R` for all inputs.

## Output results from your analysis into results folder

Separate analysis subfolders for each type.
