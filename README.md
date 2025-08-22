# UThPb-ZirChron

A toolkit for U-Th(-Pb) zircon geochronology analysis, including Bayesian eruption age estimation, data reduction schemes for LA-ICP-MS data to be implemented in iolite, and a code for synthetic data generation.

## Overview

This repository contains Python tools and workflows for processing and analyzing U-Th(-Pb) zircon geochronological data. The project focuses on:

- **Bayesian Eruption Age Estimation**: Statistical methods for determining eruption ages from zircon age datasets including uncertainties
- **Iolite Data Reduction Schemes**: Custom data reduction schemes for the iolite software
- **Synthetic Data Generation**: Tool for creating synthetic zircon age datasets for testing and validation

## Repository Structure

```
├── Bayesian eruption estimate/          # Bayesian analysis tools
│   ├── bayesian_age_calculation.py      # Core Bayesian calculation functions
│   ├── Bayesian_input_file.ipynb        # Interactive workflow notebook
│   ├── bayesian_input_file.py           # Python script version
│   └── example_data*.xlsx               # Example datasets and results
├── iolite DRS/                          # Iolite data reduction schemes
│   ├── U-Pb reduction (young zircons).py # U-Pb reduction for young zircons
│   ├── U-Th reduction.py                # U-Th reduction scheme
│   └── Fill in Labels.py                # Label filling utilities
└── Synthetic Data/                      # Synthetic data generation
    ├── definitions.py                   # Core functions for data generation
    ├── generating_synthetic_data.ipynb  # Interactive data generation notebook
    └── stored_variables.npz             # Pre-computed variables
```

## Features

### Bayesian Eruption Estimate
- **MCMC Simulation**: Markov Chain Monte Carlo methods for eruption age estimation
- **Uncertainty Quantification**: Comprehensive error propagation and uncertainty analysis
- **Visualization Tools**: Plotting functions for results interpretation

### Iolite Integration
- **Custom Data Reduction Schemes**: Specialized workflows for young zircon U-Pb dating
- **U-Th Dating Support**: Dedicated reduction scheme for U-Th geochronology
- **Quality Control Tools**: Automated label filling to avoid empty label names

### Synthetic Data Generation
- **Flexible Distribution Modeling**: Generate synthetic zircon ages from various statistical distributions
- **Realistic Uncertainty Simulation**: Apply realistic analytical uncertainties to synthetic data
- **Visualization**: Displays the synthetic data und uncertainties transparently

## Requirements

### Python Dependencies
```
numpy
pandas
scipy
matplotlib
tqdm
jupyter
openpyxl  # for Excel file handling
```

### Software Requirements
- **Python 3.7+**
- **Iolite 4** (for data reduction schemes)
- **Jupyter Notebook** (for interactive workflows)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/moserzoe/UThPb-ZirChron.git
cd UThPb-ZirChron
```

2. Install Python dependencies:
```bash
pip install numpy pandas scipy matplotlib tqdm jupyter openpyxl
```

## Usage

### Bayesian Eruption Age Analysis

1. **Prepare your data**: Ensure your zircon age data is in the correct format: see `example_data.xlsx` (needed columns: Unit, model_age, model_age_1sigma)

- `example_data.xlsx`: Sample zircon age dataset
- `example_data_bayesian_results.xlsx`: Detailed analysis results
- `example_data_bayesian_results_summary.xlsx`: Summary statistics

2. **Run the interactive notebook**: 
```bash
jupyter notebook "Bayesian eruption estimate/Bayesian_input_file.ipynb"
```

3. **Or use the Python script**:
```bash
python "Bayesian eruption estimate/bayesian_input_file.py"
```

### Iolite Data Reduction

1. Copy the desired `.py` and `.pyc` files from the `iolite DRS/` folder to your iolite directory: `iolite/Plugins/Data reduction schemes/`
2. Open iolite, data reduction schemes should appear in the DRS list
3. Process your LA-ICP-MS data according to the scheme parameters (additional information can be found in the supplementary file)

The 'U-Th reduction' data reduction scheme can not only be used for zircon, but also for groundmass glass, ilmenite, garnet, and other mineral phases.

### Synthetic Data Generation

1. **Interactive approach**:
```bash
jupyter notebook "Synthetic Data/generating_synthetic_data.ipynb"
```

## Contributing

Contributions are welcome! If something is not working for you, send me a mail and I can look into it.

## Citation

If you use this code in your research, please cite:

[Add your publication details here when available]

## Contact

**Author**: Zoe Moser  
**Institution**: ETH Zurich  
**Email**: moserz@erdw.ethz.ch

## Acknowledgments

This work was developed as part of PhD research at ETH Zurich. Thanks to my collaborators, particularly Marcel Guillong and Chetan Nathwani for their contributions to develop these codes. Special thanks to the iolite team for providing me with a free student licence.