# TULIP - TUmor CLassIfication Predictor <img src = "images/tulip.svg" alt = "tulip" width = "50" height = "50">

TULIP (the TUmor cLassIfication Predictor) classifies RNA-Seq samples into tumor types.

## Description

TULIP is a 1D convolutional neural network for classifying RNA-Seq data with [60K genes](https://github.com/CBIIT/TULIP/blob/main/gene_lists/all_genes.txt) or [19K protein coding genes](https://github.com/CBIIT/TULIP/blob/main/gene_lists/protein_coding_genes.txt). TULIP can classify either list into [17](https://github.com/CBIIT/TULIP/blob/main/labels/17_tumors.csv) or [32](https://github.com/CBIIT/TULIP/blob/main/labels/32_tumors.csv) tumor types. 

The resource transfer team trained and validated the models in TULIP on over 9,000 TCGA RNA-Seq files from the [Genomic Data Commons](https://portal.gdc.cancer.gov/) (GDC) in February 2022. To use TULIP, the user must provide a file of RNA-Seq data expressed as FPKM-UQ (fragments per kilobase of transcript per million mapped reads upper quartile) for one or more samples. TULIP then converts FPKM-UQ values to TPM (transcript per million), performs log10 normalization, and reformats the data into the correct dimensions before applying the selected model. TULIP generates two files, as described in the following table:

| Example File for 17 Tumor Types and All 60K Genes | Example File for 32 Tumor Types and Protein Coding Genes | File Contents |
| ------------- | ------------- | ------------- |
| [predictions_17_all.csv](https://github.com/CBIIT/TULIP/blob/main/example_results/predictions_17_all.csv)  | [predictions_32_pc.csv](https://github.com/CBIIT/TULIP/blob/main/example_results/predictions_32_pc.csv)  | Only the predicted primary tumor types and their probability scores. |
| [predictions_full_17_all.csv](https://github.com/CBIIT/TULIP/blob/main/example_results/predictions_full_17_all.csv)  | [predictions_full_32_pc.csv](https://github.com/CBIIT/TULIP/blob/main/example_results/predictions_full_32_pc.csv)  | The probabilitity scores for each tumor type for a full reference. |

*The naming convention of the output files reflect the model used and the file contents.*

## Software Setup

To set up the Python environment for TULIP:
1. Install [conda](https://docs.conda.io/en/latest/) package manager. 
2. Clone this repository. 

   ```bash
   git clone https://github.com/CBIIT/TULIP.git
   ```

3. Create the environment as shown below.

   ```bash
   conda env create -f environment.yml -n tulip
   conda activate tulip
   ```

These commands generate an initial folder structure with an empty *models* folder. 

## Downloading Model Weights

To download the model weights needed for running TULIP:
1. Create an account on the Model and Data Clearinghouse [MoDaC](https://modac.cancer.gov). 
2. Run the following command.

   ```bash
   python modac_file_download.py
   ```
   
3. When prompted, enter your MoDaC credentials. The script creates a *models* folder and downloads files from the [TULIP](https://modac.cancer.gov/assetDetails?dme_data_id=NCI-DME-MS01-17794660) asset in MoDaC to this folder. 

## Folder Structure

After performing the Software Setup steps and downloading the model weights, the following folder structure is available:

```
.
├── example_data/                 # folder containing example input files
├── example_results/              # folder containing example output files
├── gene_lists/                   # folder containing lists for 19K and 60K genes
├── labels/                       # folder containing lists for 17 and 32 tumor types
├── models/                       # folder containing model weights
│   ├── cnn_17_pc_weights.h5      # model weights for 17 tumor types and 19K protein coding genes
│   ├── cnn_17_weights.h5         # model weights for 17 tumor types and 60K genes
│   ├── cnn_32_pc_weights.h5      # model weights for 32 tumor types and 19K protein coding genes
│   ├── cnn_32_weights.h5         # model weights for 32 tumor types and 60K genes
├── utils                         # Python helper scripts
├── environment.yml               # Python and libraries to run TULIP
├── modac_file_download.py        # Python script to download model weights from MoDaC
├── tulip.py                      # Python script of TULIP
└── ...

```

## Data Setup

TULIP accepts the RNA-seq data expressed as FPKM-UQ, in CSV, XLSX, and TSV file formats. 

Arrange the data with the Ensembl IDs in the first column and the expression values starting from the third column, as shown in the following example:

| gene_id | gene_name | CPT0019990006 | CPT0017440009 | CPT0077290006 |
| --------| ----------|----------|----------|----------|
| ENSG00000000003.13 | TSPAN6 | 157778.5731 | 76515.8557 | 205326.5947 |
| ENSG00000000005.5 | TNMD | 2828.4868 | 3321.4867 | 5517.4428 |
| ENSG00000000419.11 | DPM1 | 508866.9116 | 332778.5383 | 468852.2266 |

The [example_data](https://github.com/CBIIT/TULIP/tree/main/example_data) folder provides example files with 60K genes and 19K protein coding genes. The CPTAC samples (from [GDC](https://portal.gdc.cancer.gov/)) included in these files represent kidney cancer. 

If the data file contains any duplicate Ensembl IDs, TULIP removes them. Additionally, if the data file is missing any Ensembl IDs, TULIP adds them and sets the expression values to 0 for each sample. 

## Running TULIP

To run TULIP, run a command in the following format:

   ```bash
   python tulip.py -i <path/to/file> [options]
   ```

Consider the following example commands:

 * 17 tumor types and 60K genes
   ```bash
   python tulip.py -i example_data/all_genes_htseq_fpkm_uq.csv -t 17 -g all -o example_results/
   ```
 * 32 tumor types and 19K protein coding genes
   ```bash
   python tulip.py -i example_data/protein_coding_genes_htseq_fpkm_uq.csv -t 32 -g pc -o example_results/
   ```
The results of the above example commands can be found in [example_results](https://github.com/CBIIT/TULIP/blob/main/example_results/).

TULIP parameters:

| Required? | Parameter | Description | Default |
| ------------- | ------------- | ------------- | ------------- |
| Yes  | -i, --input | The full path of the gene expression matrix file (FPKM-UQ) in the required format.  | (None) |
| No  | -t, --types | The number of tumor types, 17 or 32, to use for classification.  | 32 |
| No  | -g, --genes | Indicate 'all' to use all 60K genes or 'pc' to use only the 19K protein coding genes. | pc |
| No  | -o, --output_dir | The full path to the output directory. | (Current directory) |
| No  | -m, --min_score | The minimum probability score (0.0 to 1.0) for keeping the predicted primary tumor type. | (None) |

## Acknowledgments

TULIP is a newer version of [NCI-DOE-Collab-Pilot1-Tumor-Classifier](https://github.com/CBIIT/NCI-DOE-Collab-Pilot1-Tumor-Classifier).

This work has been supported in part by the Joint Design of Advanced Computing Solutions for Cancer (JDACS4C) program established by the U.S. Department of Energy (DOE) and the National Cancer Institute (NCI) of the National Institutes of Health.
