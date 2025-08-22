# UThPb-ZirChron

A toolkit for U-Th(-Pb) zircon geochronology analysis, including Bayesian eruption age estimation, iolite data reduction schemes for LA-ICP-MS data, and synthetic U-Th data generation code.

## Overview

This repository contains Python tools and workflows for processing and analyzing U-Th(-Pb) zircon geochronological data. The project focuses on:

- **Bayesian Eruption Age Estimation**: Statistical methods for determining eruption ages from zircon age datasets including uncertainties
- **Iolite Data Reduction Schemes**: Custom data reduction schemes for the iolite software
- **Synthetic Data Generation**: Tool for creating synthetic U-Th zircon age datasets for testing and validation

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
- **Bayesian eruption estimation function**: Methods for eruption age estimation following Keller et al. (2018).
- **Prior distributions**: Many different prior distribution can be applied.
- **Visualization Tools**: Plotting functions for results interpretation.

### Iolite Integration
- **Custom Data Reduction Schemes**: Specialized workflows for data processing of young zircon.
- **U-Th Dating**: Dedicated reduction scheme for U-Th geochronology following the method by Guillong et al. (2016).
- **U-Pb Dating**: Dedicated reduction scheme for young zircon U-Pb geochronology.
- **Quality Control Tools**: Independently usable automated label filling scheme to avoid empty label names.

### Synthetic Data Generation
- **Flexible Distribution Modeling**: Generate synthetic U-Th zircon ages from various underlying distributions.
- **Realistic Uncertainty Simulation**: Apply realistic analytical uncertainties to synthetic data.
- **Visualization**: Display the synthetic data und uncertainties transparently.

## Requirements

### Python Dependencies
```
numpy
pandas
scipy
matplotlib
tqdm
jupyter
openpyxl 
```

### Software Requirements
- **Python 3.7+**
- **Iolite 4** (for data reduction schemes)

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

2. **Run the interactive notebook**: 
```bash
jupyter notebook "Bayesian eruption estimate/Bayesian_input_file.ipynb"
```

3. **Or use the Python script**:
```bash
python "Bayesian eruption estimate/bayesian_input_file.py"
```

Example files include
- `example_data.xlsx`: Sample zircon age dataset
- `example_data_bayesian_results.xlsx`: Detailed analysis results of the Bayesian eruption age estimation
- `example_data_bayesian_results_summary.xlsx`: Summary of the eruption estimation results

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

We welcome contributions to improve and expand the UThPb-ZirChron toolkit. Whether you are reporting bugs, suggesting enhancements, or contributing code, your input is valuable. Please do not hesitate to reach out.

## Citation

If you use this code in your research, please cite:

[Will be added after submission]

## Contact

**Author**: Zoe Moser  
**Institution**: ETH Zurich  
**Email**: moserz@erdw.ethz.ch

## References
Guillong, M., Sliwinski, J. T., Schmitt, A., Forni, F., & Bachmann, O. (2016). U‐Th zircon dating by laser ablation single collector inductively coupled plasma‐mass spectrometry (LA‐ICP‐MS). Geostandards and Geoanalytical Research, 40(3), 377-387.
Keller, C. B., Schoene, B., & Samperton, K. M. (2018). A stochastic sampling approach to zircon eruption age interpretation. Geochemical Perspectives Letters (Online), 8(LLNL-JRNL-738859).

## Acknowledgments

This work was developed as part of PhD research at ETH Zurich. Thanks to my collaborators, particularly Marcel Guillong and Chetan Nathwani for their contributions to develop these codes. Special thanks to the iolite team for providing me with a free student licence.