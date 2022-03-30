# TULIP - TUmor CLassIfication Predictor <img src = "images/tulip.svg" width = "50" height = "50">

## Description

## Software Setup

To set up the Python environment for TULIP:
1. Install [conda](https://docs.conda.io/en/latest/) package manager. 
2. Clone this repository. 
3. Create the environment as shown below.

```bash
   conda env create -f environment.yml -n tulip
   conda activate tulip
```

## Downloading Model Weights

To download the model weights needed for running TULIP:
1. Create an account on the Model and Data Clearinghouse [MoDaC](https://modac.cancer.gov). 
2. Follow the instructions in the [Running TULIP](https://github.com/CBIIT/TULIP/blob/main/README.md#running-tulip) section below.
3. When prompted, enter your MoDaC credentials.

## Data Setup



## Running TULIP

```bash
python tulip.py -i example_data/all_genes_htseq_fpkm_uq.csv
```

```bash
python tulip.py -i example_data/all_genes_htseq_fpkm_uq.csv -t 17 -g all -o example_results/
```

## Acknowledgments

TULIP is based on [NCI-DOE-Collab-Pilot1-Tumor-Classifier](https://github.com/CBIIT/NCI-DOE-Collab-Pilot1-Tumor-Classifier).

This work has been supported in part by the Joint Design of Advanced Computing Solutions for Cancer (JDACS4C) program established by the U.S. Department of Energy (DOE) and the National Cancer Institute (NCI) of the National Institutes of Health.
