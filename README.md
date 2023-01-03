# DiSiR: a Software framework to identify ligand-receptor interactions at subunit level from single-cell RNA-sequencing data

<img src="https://github.com/miladrafiee/disir_package/blob/main/Data/ReadMe_images/disir_Method.png" width="1000">

**DiSiR** (DiSiR: Dimer Signal Receptor Analysis) is an ultrafast and easy-to-use software framework to investigate how individual cells are interacting with each other by analyzing signaling pathways of multi-subunit ligand-activated receptors from scRNA-seq data, even for the genes that are not listed in available databases. DiSiR has been built upon the framework from existing expression permutation-based methods such as CellPhoneDB (Vento-Tormo, Efremova et al. 2018, Efremova, Vento-Tormo et al. 2020) and has been specifically optimized to analyze and also visualize complex pathways involving multi-subunit receptor complexes with multi-subunit ligands for specific individual user-defined ligand-receptor input.

This repository contains the codes used for the implementation and evaluation of our method and one case study in applying it. An alpha version of a user-friendly internet-based tool for DiSiR is now publicly available [here](https://share.streamlit.io/akurlovs/disir_app/main/disir_app.py).

The primary implementation is in Python 3. To see usage example of DiSiR keep reading, and also look in the Python-module and Analysis directories. Also, the scripts to plot graph representation of resulting ligand-receptor interactions at cell type level are located in  plotting folder. In order to identify ligand-receptor interactions, 'DiSiR_main.py' script in the Python-module directory needs to be run. Below are usage examples for analyzing IL6 signaling pathway in rheumatoid arthritis (RA) Synovium scRNA-seq data and more complex IL1RAP signaling pathway in lung COVID-19 scRNA-seq data. 

## Installation instructions 
1. Install <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank">Miniconda</a> if not already available.

2. Clone this repository:
```bash
  git clone https://github.com/miladrafiee/disir_package.git.
```
3. Navigate to `disir_package` directory:
```bash
  cd disir_package
```
4. Create a conda environment with the required dependencies:
```bash
  conda env create -f environment.yml
```
5. Activate the `disir_package` environment you just created:
```bash
  conda activate disir_package
``` 
6. Install DiSiR by executing:
```bash
  pip install .
``` 

## File formatting (required inputs)
DiSiR requires 4 files as input. All files should be provided in comma-separated values input format (saved as .csv) with no double quotations. Further formatting details for each input file are specified below. 

<img src="https://github.com/miladrafiee/disir_package/blob/main/Data/ReadMe_images/disir_inputs.png" width="1000">

1. __A scRNA-seq gene expression file:__
- The matrix must be genes (rows) by cells (columns) without single cell IDs and gene names (they will be provided as another input). 
- The matrix must not have a header and row names.
- The gene expression data should be represented as non-normalized counts. 
<img src="https://github.com/miladrafiee/disir_package/blob/main/Data/ReadMe_images/scRNAfile.png" width="800"> 

2. __A cell type label file:__
- Cell type labels corresponding to the single cells (columns) in the scRNA-seq gene expression matrix with the same order. 
- The table should contain only one column which contains the cell types corresponding to the single cells in the scRNA-seq gene expression matrix with the same order.
- The columns must not have a header and row names. 
<img src="https://github.com/miladrafiee/disir_package/blob/main/Data/ReadMe_images/celllabelfile.png" width="125"> 

3. __A gene name file:__
- Genes names corresponing to the features (rows) in the scRNA-seq gene expression matrix with the same order.
- The table should contain only one column that contains the gene names corresponding to the features (rows) in the scRNA-seq gene expression matrix with the same order.
- The columns must not have a header and row names.
<img src="https://github.com/miladrafiee/disir_package/blob/main/Data/ReadMe_images/genesfile.png" width="125"> 

4. __A ligand-receptor interaction file:__
- User also needs to provide a comma-separated value (csv) file" that contains one interaction per line, for example: IL6 | IL6, IL6 | IL6ST (in this example, the IL6R and IL6ST are two receptor subunits of the ligand IL6).
- The columns must not have a header and row names.
<img src="https://github.com/miladrafiee/disir_package/blob/main/Data/ReadMe_images/interactionsfile.png" width="200"> 

Also, a name to the "output directory", in order to save the DiSiR results, needs to be provided by the user without double quotations.

## Other input parameters

1. `-threshold-fraction (-tf)` - Threshold on the fraction of cells expressed each ligand or receptor per cell type (default = 0)
 
2. `-threshold-expression (-te)` - Threshold on scaled (max-normalized) average expression of each ligand or receptor within a cell type (default = 0)
  
3. `-threshold-pvalue (-tp)` - Threshold on p-value for filtering non-significant LR interactions (default = 0.05)
   
4. `-iterations-value (-iv)` - Number of iterations for permutating data in statistical test (default = 1000)

To run DiSiR by changing the default parameter values, pass each parameter followed by the correspondong value, for example: `-tf 0.1`, `-te 0.1`, `-tp 0.01` or `-iv 500`.

## Running DiSiR
After activating the `disir_package` conda environment via `conda activate disir_package`, DiSiR can be called from the command line from any folder using `disir_package`. Examples on how to run DiSiR are provided in the section "Example datasets for running DiSiR" below.

A typical DiSiR run with default settings would look like this: 
 ```bash
 disir_package --scRNA-path /path/to/scRNA_gene_expression
    --cell-type-path /path/to/scRNA_cell_labels
    --gene-path /path/to/gene_names
    --subunit-interactions-path /path/to/ligand_receptor_interactions
    --output-directory-path path/to/output_results_directory
```
Or with more condensed parameter names: 
 ```bash
 disir_package --sp /path/to/scRNA_gene_expression
    --ctp /path/to/scRNA_cell_labels
    --gp /path/to/gene_names
    --sip /path/to/ligand_receptor_interactions
    --odp path/to/output_results_directory
```

For full usage details with additional options, see the section "Other input parameters" above. 

## DiSiR outputs
DiSiR will produce the following output files by default.
1. ```expression_info.pdf``` Bubble plot illustrating max-normalized average expressions of the putative signaling pathway components per cell type (color of the circles), and fraction of cells expressing them within its corresponding cell type (size of the circles; adjusted by the total number of cells expressing the respective gene). 
2. ```Links.csv``` Significant ligand-receptor interactions, and their corresponding strengths, identified by DiSiR. These interactions are the links of the resulting inetraction graph.
3. ```Nodes.csv``` Cell type names associated with the identified LR intercations listed in the "Links.csv" file.
4. ```Heatmap.csv``` and ```Heatmap.pdf``` Heatmap plot illustrating significant interactions between different cell types listed in rows and columns of the heatmaps. The color intensity in heatmap representation is associated with the strength of interactions between cell types.
5. ```Heatmap_all_interactions.csv``` and ```Heatmap_all_interactions.pdf``` Heatmap plot illustrating all interactions between different cell types listed in rows and columns of the heatmaps. The color intensity in heatmap representation is associated with the strength of interactions between cell types.

## Example datasets for running DiSiR
For users to test DiSiR, we have included files for two example runs:
1. [AMP consortium’s Phase 1 rheumatoid arthritis (RA) data](https://immunogenomics.io/ampra/) published in [Zhang, Wei et al. 2019](https://www.nature.com/articles/s41590-019-0378-1). Please find the gene expression matrix, genes list and meta data [here](https://drive.google.com/drive/folders/14hHjAMff3i_nFr4zfkptre-1g9i-s2WF?usp=share_link).
2. [Lung COVID-19 scRNA-seq data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926) published in [Liao, Liu et al. 2020](https://www.nature.com/articles/s41590-019-0378-1). You can find the gene expression matrix, genes list and meta data [here](https://drive.google.com/drive/folders/1UiFeTFB7My7UdsRmuZGHRD6yAM8O-0zV?usp=share_link).

Two of LR subunit interaction files (for IL6/IL11 and IL1RAP signaling pathways) are also provided [here](https://drive.google.com/drive/folders/1yLoV3ai51kI2KWE27rh9beN5Dg3R3Cir?usp=share_link) as an example.

### Commands for running example analyses:
Once the example files are downloaded, the commands below can be run from the folders where the example datasets are located (for IL6/IL11 pathway):
1. RA synovium data:
```bash
  disir_package -sp scRNA.csv -gp genes.csv -sip Subunit_interactions_IL6_IL11.csv -ctp cell_type_labels.csv -odp results_RA
```
2. COVID data: 
Since the input scRNA-seq data is very big in this case, we have saved the input data as a .mtx matrix (scRNA.mtx) which used for big sparse matrices. DiSiR can handle this file format by setting the parameter `sparse-matrix (-sm)` as below:
```bash
  disir_package -sm -sp scRNA.mtx -gp genes.csv -sip Subunit_interactions_IL6_IL11.csv -ctp cell_type_labels.csv -odp results_COVID 
```

## DiSiR outputs for example RA synovium scRNA-seq data and IL6 signaling pathway
The main outputs from a DiSiR run are the files named `Links.csv` and `Nodes.csv`, which provide the significant ligand-receptor interactions (and their corresponding strengths) identified by DiSiR and cell type names associated with them, respectively.

The DiSiR results are visualized in heatmaps saved as `Heatmap_all.pdf` and `Heatmap.pdf` illustrating all and significant interactions between different cell types listed in rows and columns of the heatmaps. The color intensity in heatmap representation is associated with the strength of interactions between cell types. Below are the heatmaps produced for the example RA data.

<p align="center">
  <img src="https://github.com/miladrafiee/disir_package/blob/main/Data/ReadMe_images/heatmaps.png" width="1000"> 
</p>

DiSiR also visualizes output cell-cell interactions in two ways: graph representation and chord diagram representation. The outputs of running "DiSiR_main.py" are the links and nodes of the resulting interaction interaction graph at cell type level, and heatmap plots for all and significant interactions. Using the DiSiR outputs "links.csv" and "Nodes.csv" (along with the putative LR interactions list) files as the inputs of the "Graph_representation.R" script, which is located in the "Plotting" directory, users can generate a directed graph in which nodes are associated with the cell types present in the input data set and each edge corresponds to a ligand–receptor interaction between two cell types (from ligand-expressing cells to receptor-expressing cells). For a given interaction, if both ligand and receptor are present in the same cell type, then there is a self-loop in the graph on the corresponding node. We use the “visNetwork version 2.1.0” package in R version 4.0.0 with an interactive environment.

Also, using the DiSiR output "Heatmap.csv" file as the input of the "Chord_diagram_representation.R" script, which is located in the "Plotting" directory, users can generate a chord diagram representation of the resulting interactions at cell type level. We use the “circlize version 0.4.15” package in R for this purpose.

<p align="center">
  <img src="https://github.com/miladrafiee/disir_package/blob/main/Data/ReadMe_images/graph_chord_representation.png" width="1000"> 
</p>
  
## Authors
DiSiR was developed by

* Milad R. Vahid
* Andre Kurlovs 
* Virginia Savova 
* Franck Rapaport
* Tommaso Andreani
* Franck Augé

## License
DiSiR is licensed under the GNU GPL, version 3 or (at your option) any
later version.
DiSiR is Copyright (2022-) by the authors.

## Citation
If you use DiSiR, please cite:  
  
*DiSiR: a software framework to identify ligand-receptor interactions at subunit level from single-cell RNA-sequencing data* (bioxRiv 2022)  
Milad R. Vahid, Andre Kurlovs, Franck Auge, Reza Olfati-Saber, Emanuele de Rinaldis, Franck Rapaport, Virginia Savova . 
doi: https://doi.org/10.1101/2022.03.25.485741
