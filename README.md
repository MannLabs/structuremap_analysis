# Analysis of PTMs in their structural context using the StructureMap and AlphaMap Python packages

This repository contains all data and python notebooks to reproduce the analyses presented in: 
["The structural context of PTMs at a proteome wide scale" by Bludau et al. (2022)](https://doi.org/10.1101/2022.02.23.481596).

To set up the analysis, the first step is to run the commands in [setup_environment.sh](https://github.com/MannLabs/structuremap_analysis/blob/master/setup_environment.sh).

The actual analysis is performed in multiple jupyter notebooks:
* PTM data processing and import: notebooks in [ptm_data_import](https://github.com/MannLabs/structuremap_analysis/blob/master/ptm_data_import)
* Benchmark for IDR prediction: [IDR_benchmark.ipynb](https://github.com/MannLabs/structuremap_analysis/blob/master/IDR_benchmark.ipynb)
* Main structuremap analysis: [data_analysis_structuremap.ipynb](https://github.com/MannLabs/structuremap_analysis/blob/master/data_analysis_structuremap.ipynb)
* Proteasome inhibitor usage in ubiquitination datasets on PhosphoSitePlus: [PSP_Ubi.ipynb](https://github.com/MannLabs/structuremap_analysis/blob/master/PSP_Ubi.ipynb)
