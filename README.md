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
7. [License](#license)

---

## Overview

The pipeline performs the following steps:
1. **Normalization**: Normalizes lipidomics data using a user-specified method (TIC, Median, or PQN).
2. **Log Transformation**: Applies log2 transformation to the normalized data.
3. **Differential Analysis**: Uses LIMMA to perform differential analysis based on user-provided contrasts.
4. **Visualization**: Generates RLE plots (before and after normalization), volcano plots, and heatmaps.

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/lipidomics-analysis.git
   cd lipidomics-analysis
```
