# Lipidomics Data Analysis Pipeline

This repository contains an R script for normalizing and analyzing lipidomics data. The pipeline includes normalization, log transformation, differential analysis using LIMMA, and visualization of results through RLE plots, volcano plots, and heatmaps.

---

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Input Files](#input-files)
5. [Output Files](#output-files)
6. [Dependencies](#dependencies)
7. [Contact](#contact)

---

## Overview

The pipeline performs the following steps:
1. **Normalization**: Normalizes lipidomics data using a user-specified method (TIC, Median, or PQN).
2. **Log Transformation**: Applies log2 transformation to the normalized data.
3. **Differential Analysis**: Uses LIMMA to perform differential analysis based on user-provided contrasts.
4. **Visualization**: Generates RLE plots (before and after normalization), volcano plots, and heatmaps.

---

## Installation

### 1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/lipidomics-analysis.git
   cd lipidomics-analysis
   ```
### 2. Install the required R packages:
   ```bash
   install.packages(c("ggplot2", "dplyr", "tidyr", "limma", "pheatmap", "RColorBrewer"))
   ```
---

## Usage

### 1. Run the script from the command line with the following arguments:
```bash
Rscript lipidomics_analysis.R <input_data.csv> <contrasts_file.csv> <output_dir> <normalization_method>
```

### Arguments

- ```<input_data.csv>```: Path to the raw lipidomics data (CSV file).
- ```<contrasts_file.csv>```: Path to the contrasts file (CSV file with a column named ```Contrast```).
- ```<output_dir>```: Path to the output directory.
- ```<normalization_method>```: Normalization method (e.g., ```TIC```, ```Median```, ```PQN```).

**Example**:
```bash
Rscript lipidomics_analysis.R data/input_data.csv data/contrasts.csv results TIC
```

---

## Input Files

**Raw Lipidomics Data** (```input_data.csv```)
- **Format**: CSV file with samples as rows and lipids as columns.
- **Required Columns**:
   - ```Sample Name```: Unique identifier for each sample.
   - ```Sample Group```: Group or condition for each sample.
   - Additional columns: Lipid intensities.

**Example:**
| Sample Name | Sample Group | Lipid1 | Lipid2 | Lipid3 |
|-------------|--------------|--------|--------|--------|
| Sample1     | Group1       | 100    | 200    | 300    |
| Sample2     | Group2       | 150    | 250    | 350    |
| Sample3     | Group1       | 120    | 220    | 320    |

**Contrasts File** (```contrasts_file.csv```)
- **Format**: CSV file with a single column named ```Contrast```.
- Each row specifies a contrast for differential analysis.

**Example**:
| Contrast | 
|-------------|
| Group1-Group2     | 
| Group3-Group4     |
---

## Output Files
```
output_dir/
├── Normalization_Results/
│   ├── RLE_Before_Normalization.pdf
│   ├── RLE_After_Normalization.pdf
│   └── Normalized_Data.csv
├── Differential_Analysis_Results/
│   ├── Volcano_Plots/
│   │   ├── Group1-Group2_Volcano.pdf
│   │   └── Group3-Group4_Volcano.pdf
│   ├── Heatmaps/
│   │   ├── Group1-Group2_Heatmap.pdf
│   │   └── Group3-Group4_Heatmap.pdf
│   └── Group1-Group2_Results.csv
│   └── Group3-Group4_Results.csv
```
---

## Dependencies
The script requires the following R packages:

- ggplot2
- dplyr
- tidyr
- limma
- pheatmap
- RColorBrewer

Install them using:
```R
install.packages(c("ggplot2", "dplyr", "tidyr", "limma", "pheatmap", "RColorBrewer"))
```
---

## Contact
For questions or feedback, please contact Fayrouz Hammal at fayrouz.hammal@petermac.org.
