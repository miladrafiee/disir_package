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

For full usage details with additional options, see the section "Other input parameters" below. 

## DiSiR outputs
DiSiR will produce five output files by default.
1. ```plot_cell_type_locations.pdf``` Heatmaps of cell type assignments within the ST sample. Along with a plot showing the total number of cells mapped to each spot, these show the spatial distribution of cell type assignments. Color bars indicate the number of cells of the respective cell type inferred per spot.
2. ```assigned_locations.csv``` This file will provide the assigned locations of each single cell mapped to ST spots. As some cells may be mapped to multiple locations depending on the size of the input scRNA-seq set, new cell IDs (`UniqueCID`) are assigned to each cell and given in the first column. The second column includes original cell IDs (`OriginalCID`); the third includes assigned spot IDs (`SpotID`); the fourth and fifth columns respectively include  `row` and `column` indices of the corresponding spots; and then optionally, the sixth and seventh columns include `coord_x` and `coord_y` of the corresponding spots if these details were provided in the initial Coordinates file.
3. ```cell_type_assignments_by_spot.csv``` This file gives the raw number of cells of each cell type per spot by `SpotID` as well as the total number of cells assigned to that spot.
4. ```fractional_abundances_by_spot.csv``` This file gives the fractional abundance of cell types assigned to each spot by `SpotID`.
5. ```log.txt``` This file contains a log of CytoSPACE run parameters and running time.

## Example datasets for running CytoSPACE
For users to test CytoSPACE, we have included files for two example runs:
1. A HER2+ breast cancer scRNA-seq atlas by Wu et al. (<a href="https://www.nature.com/articles/s41588-021-00911-1" target="_blank">Nature Genetics, 2021</a>) and a HER2+ breast cancer FFPE specimen profiled by the Visium platform (<a href="https://www.10xgenomics.com/resources/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0" target="_blank">10x Genomics</a>). Default parameters were selected with Visium samples in mind and are appropriate here.
2. A melanoma scRNA-seq atlas by Tirosh et al. (<a href="https://www.science.org/doi/10.1126/science.aad0501?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed" target="_blank">Science, 2016</a>), and a melanoma specimen profiled by the legacy ST platform (Thrane et al, <a href="https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables" target="_blank">Cancer Research, 2018</a>). As the legacy ST platform has larger spot sizes, we recommend mapping an average of 20 cells per spot, i.e. passing the argument `-mcn 20` to your `cytospace` call.

### Download example datasets
Zip files containing the example datasets can be downloaded from the following links:
1. <a href="https://drive.google.com/file/d/1ODN7Ux2x9XNl1K8cOLl4JxIZpjJL2GwG/view?usp=sharing" target="_blank">Breast cancer</a>
2. <a href="https://drive.google.com/file/d/1oHe4UP2K0kQS9gNFVvZtjqtJyJK_VeBp/view?usp=sharing" target="_blank">Melanoma</a> 

To download from the command line using `gdown`:
1. Breast cancer
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1ODN7Ux2x9XNl1K8cOLl4JxIZpjJL2GwG/view?usp=sharing
   unzip CytoSPACE_example_breast_cancer.zip
   ```
2. Melanoma
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1oHe4UP2K0kQS9gNFVvZtjqtJyJK_VeBp/view?usp=sharing
   unzip CytoSPACE_example_melanoma.zip
   ```
   
### Commands for running example analyses:
Once the example files are downloaded, the commands below can be run from the folders where the example datasets are located:
```bash
  cytospace -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp brca_STdata_GEP.txt -cp brca_STdata_coordinates.txt -ctfep brca_cell_fraction_estimates.txt -o cytospace_results_brca -sm lap_CSPR
```

```bash
  cytospace -sp melanoma_scRNA_GEP.txt -ctp melanoma_scRNA_celllabels.txt -stp melanoma_STdata_slide1_GEP.txt -cp melanoma_STdata_slide1_coordinates.txt -ctfep melanoma_cell_fraction_estimates.txt -mcn 20 -o cytospace_results_melanoma -nr 5 -nc 3 -ss 1100 -pm s -nv -sm lap_CSPR 
```
Please note that here we use the `lap_CSPR` solver for compatibility. If your system supports AVX2 intrinsics, you can run the same commands without the final argument to use the `lapjv` solver instead. __These CytoSPACE runs should take around 5 minutes each.__

### CytoSPACE output files for example breast cancer data
The main output from a CytoSPACE run is the file named `assigned_locations.csv`, which provides the ST spots to which the single cells have been assigned. 

<img width="600" src="https://github.com/digitalcytometry/cytospace/blob/main/images/assigned_locations.png">


The CytoSPACE results are visualized in heatmaps saved as `plot_cell_type_locations.pdf` showing the distribution of single cells across ST spots for each cell type. Color bars indicate the number of cells of the respective cell type inferred per spot. Below are the heatmaps produced for the example BRCA data.

<p align="center">
  <img width="800" src="https://github.com/digitalcytometry/cytospace/blob/main/images/BRCA_plot_cell_type_locations.png">
</p>

For comparison, consider the pathologist annotations of this ST sample as provided by 10x:

<p align="center">
  <img width="800" src="https://github.com/digitalcytometry/cytospace/blob/main/images/Visium_FFPE_Human_Breast_Cancer_Pathologist_Annotations.png">
</p>

The number of cells per spot by cell type as well as in total are provided in the file `cell_type_assignments_by_spot.csv`. Fractional abundances of each cell type are returned in the file `fractional_abundances_by_spot.csv`. A log file recording CytoSPACE inputs and running times is output in the file `log.txt`.

Zip files of expected CytoSPACE outputs (with `lap_CSPR` solver) are available to download at the following links:
1. <a href="https://drive.google.com/file/d/1CLfy4Txez8ThID8YzIH04hlvrBRCQ4Rh/view?usp=sharing" target="_blank">Breast cancer results</a>
2. <a href="https://drive.google.com/file/d/1X4jMwctRNmqCRIJcop2hL3jRdlhxnnlc/view?usp=sharing" target="_blank">Melanoma results</a> 

To download from the command line using `gdown`:
1. Breast cancer
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1CLfy4Txez8ThID8YzIH04hlvrBRCQ4Rh/view?usp=sharing
   unzip CytoSPACE_example_breast_cancer_results.zip
   ```
2. Melanoma
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/1X4jMwctRNmqCRIJcop2hL3jRdlhxnnlc/view?usp=sharing
   unzip CytoSPACE_example_melanoma_results.zip
   ```

## Customizing plotting outputs
CytoSPACE provides two mechanisms for plotting output heatmaps. First, CytoSPACE generates heatmaps by default within the main function call.  By default, the geometry and plotting parameters have been optimized for standard 10x Visium spatial slides. To plot non-Visium spatial data, the `plot-nonvisium` flag should be set to `True`. As an example, to change the plotting parameters, e.g., `-nr` (number of rows of heatmaps per page), `-nc` (number of columns of heatmaps per page), `-ss` (spot size), -pm (plot marker which indicates the shape of the spots, e.g., `-pm s` plots squares and `-pm h` plots hexagons) and `-nv` (plot non-Visium), from command line, for the melanoma dataset:
```bash
  cytospace -sp melanoma_scRNA_GEP.txt -ctp melanoma_scRNA_celllabels.txt -stp melanoma_STdata_slide1_GEP.txt -cp melanoma_STdata_slide1_coordinates.txt -ctfep melanoma_cell_fraction_estimates.txt -nr 5 -nc 3 -ss 1100 -pm s -nv
```
`cytospace-plot` can be used to generate heatmaps from CytoSPACE outputs. There are three paths that need to be set by users: `-alp` (assigned_locations_path), `-cp` (coordinates_path) and `-o` (output_filename). After CytoSPACE outputs have been generated, `cytospace-plot` can be run from the command line for the sample breast cancer data as:
```bash
  cytospace-plot -alp cytospace_results_brca/assigned_locations.csv -cp  brca_STdata_coordinates.txt -o brca_results.pdf
```
and for the sample melanoma data as: 
```bash
  cytospace-plot -alp cytospace_results_melanoma/assigned_locations.csv -cp  melanoma_STdata_slide1_coordinates.txt -o melanoma_results.pdf -nr 5 -nc 3 -ss 1100 -pm s -nv
```

## Advanced options
While default options are recommended for most use cases, we do provide additional advanced options.

### User-provided estimates of number of cells per spot
Rather than using the internal mechanism of CytoSPACE for estimating the number of cells per spot, users can provide their own estimates (from image segmentation, for example) in a two-column file with header, in which the first column contains spot IDs and the second contains the number of cells predicted per spot:

<p align="center">
  <img width="300" src="https://github.com/digitalcytometry/cytospace/blob/main/images/n_cells_per_spot.PNG">
</p>

To run CytoSPACE with this option, pass the flag `-ncpsp` or `--n-cells-per-spot-path` followed by the file location.

### Alternative distance metric
By default, CytoSPACE uses Pearson correlation to compare cell and spot transcriptomes. Users can choose to use Spearman correlation or Euclidean distance instead by passing `-dm Spearman_correlation` or `-dm Euclidean` respectively with the function call. 

### Setting a new random seed
While the CytoSPACE algorithm is mostly deterministic, the initial step of sampling cells to be mapped is done at random. To provide an alternative random seed resulting in a different random sampling of cells, users can pass `-se` followed by the desired (integer) seed with the function call. The default random seed for CytoSPACE is 1.

### Alternative handling of sampling
CytoSPACE starts by creating a pool of cells that matches what is expected within the ST data. By default, this is done by resampling single cells to achieve the overall cell type fractions and total cell numbers estimated in the tissue. We recommend that CytoSPACE be run with this default setting for all real data analyses. However, we provide an additional option to generate new "place-holder" cells by sampling from the distribution of gene counts within each cell type instead, and used this option for ensuring uniqueness of mapped cells for benchmarking on simulated data. To run CytoSPACE with this alternative mode, users can pass `-sam place_holders` with the function call. 

### Method extension: mapping quality
While CytoSPACE's formulation as a linear assignment problem guarantees an optimal solution given its cost function, there is no underlying probabilistic framework for estimating mapping uncertainty. One possibility is to determine whether a given cell type belongs to a given spot after mapping - that is, whether a spot contains at least one cell of the same cell type. Notably, this does not distinguish between cells of the same cell type for quality of fit. As such a protocol provides some measure of mapping quality, albeit incomplete, we provide a helper script that implements this via a support vector machine that produces and trains on pseudo-bulks generated from the input scRNA-seq data. This script, `uncertainty_quantification.R`, takes as input the path to the ST dataset count matrix file, the scRNA-seq count matrix file, and the CytoSPACE output file `assigned_locations.csv`, and returns an appended output file with confidence scores in `assigned_locationswConfidenceScores.csv`. The command to run this script following a completed CytoSPACE run is as follows: 
 ```bash
 Rscript uncertainty_quantification.R /path/to/ST_geneexpression /path/to/scRNA_geneexpression /path/to/assigned_locations.csv
```
For interpreting confidence scores, we recommend a cutoff of 0.1, with higher scores indicating increased confidence that a spot contains at least one cell of the same cell type.

Please note that `uncertainty_quantification.R` requires separate dependencies from those included in the provided `environment.yml` file for the `cytospace` conda environment. This script should be run in a separate environment with the following R packages installed: `Seurat` (must be v4; tested with v4.0.1), `data.table` (tested with v1.14.0), and `e1071` (tested with v1.7.8).

### Method extension: single cell ST data
While designed for Visium-type data in which most spots contain RNA from multiple cells, CytoSPACE can also be used with single-cell resolution spatial data such as <a href="https://vizgen.com/resources/meet-the-merscope-platform/" target="_blank">Vizgen's MERSCOPE platform</a>. We expect this extension to be useful for reducing noise and expanding transcriptome coverage of each cell in the ST data. For this single-cell resolution mode, CytoSPACE partitions the ST data into smaller chunks and utilizes multiple CPU cores to assign down-sampled versions of the reference scRNA-seq data to these regions.

To run CytoSPACE with single-cell resolution spatial data:
 ```bash
 cytospace --single-cell 
    --scRNA-path /path/to/scRNA_geneexpression
    --cell-type-path /path/to/scRNA_celllabels
    --st-path /path/to/ST_geneexpression
    --coordinates-path /path/to/ST_coordinates
    --cell-type-fraction-estimation-path path/to/cellfracestimates
    --number-of-processors NUMBER_OF_PROCESSORS
    --number-of-selected-cells  NUMBER_OF_SELECTED_CELLS
    --number-of-selected-spots NUMBER_OF_SELECTED_SPOTS
```
Or with more condensed parameter names: 
 ```bash
 cytospace -sc
    -sp /path/to/scRNA_geneexpression
    -ctp /path/to/scRNA_celllabels
    -stp /path/to/ST_geneexpression
    -cp /path/to/ST_coordinates
    -ctfep path/to/cellfracestimates
    -nop NUMBER_OF_PROCESSORS
    -nosc NUMBER_OF_SELECTED_CELLS
    -noss NUMBER_OF_SELECTED_SPOTS
```
where `NUMBER_OF_PROCESSORS` denotes the number of cores to use, `NUMBER_OF_SELECTED_CELLS` denotes the number of cells in each dowsampled version, and `NUMBER_OF_SELECTED_SPOTS` denotes the size of each ST region. We generally recommend `-nosc 10000 -noss 10000`.

A zip file of example single cell inputs is available to download from Google Drive <a href="https://drive.google.com/file/d/10fhxjCn-VfPPurrI-RE8lbs6NCPqfGXY/view?usp=sharing" target="_blank">here</a>. To download from the command line using `gdown`:
   ```bash
   gdown --fuzzy https://drive.google.com/file/d/10fhxjCn-VfPPurrI-RE8lbs6NCPqfGXY/view?usp=sharing
   unzip single_cell_example_data.zip
   ```

To run CytoSPACE with this example dataset, run the following command from the location of the unzipped inputs and with your CytoSPACE conda environment active:
 ```bash
 cytospace -sc -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp single_cell_ST_example_GEP.txt -cp single_cell_ST_example_coors.txt -ctfep single_cell_ST_example_fractions.txt -nop 1 -nosc 10000 -noss 10000
```

## Extended usage details
```
usage: cytospace [-h] -sp SCRNA_PATH -ctp CELL_TYPE_PATH -stp ST_PATH -cp
                 COORDINATES_PATH -ctfep CELL_TYPE_FRACTION_ESTIMATION_PATH
                 [-ncpsp N_CELLS_PER_SPOT_PATH] [-o OUTPUT_FOLDER]
                 [-op OUTPUT_PREFIX] [-d DELIMITER]
                 [-sm {lapjv,lapjv_compat,lap_CSPR}]
                 [-sam {duplicates,place_holders}] [-mcn MEAN_CELL_NUMBERS]
                 [-se SEED]
                 [-dm {Pearson_correlation,Spearman_correlation,Euclidean}]
                 [-nosc NUMBER_OF_SELECTED_CELLS]
                 [-noss NUMBER_OF_SELECTED_SPOTS] [-nop NUMBER_OF_PROCESSORS]
                 [-sc] [-p] [-nr NUM_ROW] [-nc NUM_COLUMN] [-r] [-nv]
                 [-rd ROTATION_DEGREES] [-ss SPOT_SIZE] [-pm PLOT_MARKER]

CytoSPACE is a computational strategy for assigning single-cell transcriptomes
to in situ spatial transcriptomics (ST) data. Our method solves single
cell/spot assignment by minimizing a correlation-based cost function through a
linear programming-based optimization routine.

optional arguments:
  -h, --help            show this help message and exit
  -ncpsp N_CELLS_PER_SPOT_PATH, --n-cells-per-spot-path N_CELLS_PER_SPOT_PATH
                        Path to number of cells per ST spot file
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Relative path to the output folder
  -op OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix of results stored in the 'output_folder'
  -d DELIMITER, --delimiter DELIMITER
                        Set delimiter of the input files, default '\t'
  -sm {lapjv,lapjv_compat,lap_CSPR}, --solver-method {lapjv,lapjv_compat,lap_CSPR}
                        Which solver to use for the linear assignment problem,
                        default 'lapjv'
  -sam {duplicates,place_holders}, --sampling-method {duplicates,place_holders}
                        Which unerlying method to use for dealing with
                        duplicated cells, default 'duplicates'
  -mcn MEAN_CELL_NUMBERS, --mean-cell-numbers MEAN_CELL_NUMBERS
                        Mean number of cells per spot, default 5 (appropriate
                        for Visium). If analyzing legacy spatial
                        transcriptomics data, set to 20
  -se SEED, --seed SEED
                        Set seed for random generators, default 1
  -dm {Pearson_correlation,Spearman_correlation,Euclidean}, --distance-metric {Pearson_correlation,Spearman_correlation,Euclidean}
                        Which distance metric to use for the cost matrix,
                        default 'Pearson_correlation'
  -nosc NUMBER_OF_SELECTED_CELLS, --number-of-selected-cells NUMBER_OF_SELECTED_CELLS
                        Number of selected cells from scRNA-seq data used in
                        each iteration
  -noss NUMBER_OF_SELECTED_SPOTS, --number-of-selected-spots NUMBER_OF_SELECTED_SPOTS
                        Number of selected spots from ST data used in each
                        iteration
  -nop NUMBER_OF_PROCESSORS, --number-of-processors NUMBER_OF_PROCESSORS
                        Number of processors used for the analysis
  -sc, --single-cell    Use single-cell spatial approach or not
  -p, --plot-off        Turn create plots on/off
  -nr NUM_ROW, --num-row NUM_ROW
                        Number of rows in pdf figure
  -nc NUM_COLUMN, --num-column NUM_COLUMN
                        Number of columns in pdf figure
  -r, --rotation-flag   Rotate plot
  -nv, --plot-nonvisium
                        Plot with custom slide dimensions
  -rd ROTATION_DEGREES, --rotation-degrees ROTATION_DEGREES
                        Rotation on plot
  -ss SPOT_SIZE, --spot-size SPOT_SIZE
                        Set size of ST spots
  -pm PLOT_MARKER, --plot-marker PLOT_MARKER
                        Shape of ST spots

Required arguments:
  -sp SCRNA_PATH, --scRNA-path SCRNA_PATH
                        Path to scRNA-Seq data
  -ctp CELL_TYPE_PATH, --cell-type-path CELL_TYPE_PATH
                        Path to cell type labels
  -stp ST_PATH, --st-path ST_PATH
                        Path to spatial transcriptomics data (expressions)
  -cp COORDINATES_PATH, --coordinates-path COORDINATES_PATH
                        Path to transcriptomics data (coordinates)
  -ctfep CELL_TYPE_FRACTION_ESTIMATION_PATH, --cell-type-fraction-estimation-path CELL_TYPE_FRACTION_ESTIMATION_PATH
                        Path to cell type fraction file
```

You can see this list of variables and default values for running CytoSPACE from the commmand line as well at any time by calling `cytospace` along with the `-h` or 
`--help` flag, i.e., `cytospace -h`.

## CytoSPACE Solver options
1. `lapjv` __(Recommended for most systems)__    By default, CytoSPACE calls the `lapjv` solver from package `lapjv`. This solver is a fast implementation of the Jonker-Volgenant shortest augmenting path assignment algorithm and returns a globally optimal solution given the objective function as defined in our paper [cite]. As noted above, however, this package is not supported on all systems as it achieves its speedup through use of AVX2 instructions. This solver will be selected by default and can be specified explicitly by passing arguments `--solver-method lapjv` or `-sm lapjv` to `cytospace`.
2. `lap_CSPR` __(Recommended for systems not supporting `lapjv`)__    A second solver option is the `linear_assignment` method from the `ortools` package. This solver uses a different method than the first and third options, an assignment algorithm called the cost scaling push relabel method. This algorithm approximates assignment costs to integer values and loses some numerical precision in doing so. Therefore, while it returns a globally optimal solution __after approximation__ given the objective function defined in the paper, it will return similar but generally not identical results to the first two methods. This solver has a similar running time to the first option and is a good option for systems not supporting the `lapjv` package. This solver can be selected by passing arguments `--solver-method lap_CSPR` or `-sm lap_CSPR` to `cytospace`.
3. `lapjv_compat`   A third solver option implements the `lapjv` solver from package `lap`. Like the first option `lapjv`, this solver also implements the Jonker-Volgenant shortest augmenting path assignment algorithm to return the same globally optimal solution given the objective function defined in the paper. Furthermore, it is broadly supported and should work on all standard operating systems. However, it takes 3-4 times as long to run as the first solver option, the `lapjv` solver from the `lapjv` package, so we only recommend it for systems that do not support the first option. This solver can be selected by passing arguments `--solver-method lapjv_compat` or `-sm lapjv_compat` to `cytospace`.

## Updating local installations
To update your local installation of CytoSPACE following updates of this GitHub repository, navigate to your `cytospace` directory and execute the following commands:
```bash
  git pull
  conda env update --name cytospace --file environment.yml
  conda activate cytospace
  pip install .
```
If you have made local updates to your version of the CytoSPACE source code, you should execute 
```bash
  pip install .
``` 
before running. 

## Authors
CytoSPACE was developed in the <a href="https://anlab.stanford.edu/" target="_blank">Newman Lab</a> by

* Milad R. Vahid (miladrv)
* Erin L. Brown (erinlbrown)
* Chloé B. Steen (cbsteen)
* Aaron M. Newman (aaronmnewman)

## Contact
If you have any questions, please contact the CytoSPACE team at cytospaceteam@gmail.com.

## License
CytoSPACE is licensed under the GNU GPL, version 3 or (at your option) any
later version.
CytoSPACE is Copyright (2022-) by the authors.

## Citation
If you use CytoSPACE, please cite:  
  
*Robust alignment of single-cell and spatial transcriptomes with CytoSPACE* (bioxRiv 2022)  
Milad R. Vahid*, Erin L. Brown*,  Chloé B. Steen*,  Minji Kang,  Andrew J. Gentles,  Aaron M. Newman.  
doi: https://doi.org/10.1101/2022.05.20.488356
