
## Introduction

SAMMPy (Sensitivity Analysis for Multiple Models using Python) is an open-source python package for performing process sensitivity analysis to identify controlling processes under process model uncertainty that a process may be represented by multiple models. This is to advance our understanding of the impacts of process representations on model outputs. SAMMPy can generate a range of process sensitivty indeices including:

* The first-order process sensitivity index to determine which process(es) to be prioritized for achieving the greatest reduction in predictive uncertainty [(Dai et al., 2017)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019715);
* The total-effect process sensitivity index to determine which process(es) to be fixed and excluded from further study (Yang et al., 2020, under preparation);
* The mean and variacne of the output difference of the porcesses by using Multi-model difference-based sensitivity (MMDS) method to screen non-influential process(es) (Yang et al., 2020, under preparation).

SAMMPy is a modular modeling code that can simply and efficiently vary model structure (process representation), allowing for the generation and running of large model ensembles that vary in process representation, parameters, parameter values, and environmental conditions. The repository file contains:

```bash
├── analyze                                         
│   ├── mmds.py                                     # Code for performing mutil-model difference-based process sensitivity analysis
│   ├── vbsa.py                                     # Code for performing variance-based process sensitivity analysis
├── examples         
│   ├── workflow_1D_groundwater_flow.ipynb          # Example workflow of 1D groundwater flow model
│   ├── workflow_Sobol_G_star.ipynb                 # Example workflow of Sobol-G* function
├── plotting
│   ├── bar.py                                      # Code for generating bar plots for first-order and total effect process sensitivity indices    
│   ├── dotty.py                                    # Code for generating scatter plots for mean vs. variance of the output difference  
│   ├── hist.py                                     # Code for generating histograms for sampled parameter values 
├── util
│   ├── results.py                                  # Code for postprocessing the sensitivity results into a dictionary
├── __init__.py                                     # Common script used in the regular package, in which a model class is defined
├── LICENCE     
├── README.md         
└── setup.py                                        # Centre script of and installing this package
```

## Installation

SAMMPy requires **Python** 3.7 (or higher)

The easiest way to install is via `pip`:

To install SAMMPy type:

    pip install sammpy

To update SAMMPy type:

    pip install sammpy --upgrade

To uninstall SAMMPy type:

    pip uninstall sammpy
    
## Requirements:

- [NumPy](https://www.numpy.org)
- [Numba](http://numba.pydata.org)
- [Matplotlib](https://www.scipy.org/scipylib)
- [SALib](https://salib.readthedocs.io/en/latest/)
    
## Getting started with the package

We recommend to start by executing the workflow scripts that apply the different process sensitivity analysis methods to test cases. 

## License

SAMMPy is distributed under the MIT License. See the [LICENSE](https://github.com/jyangfsu/SAMMPy/LICENSE) file for details.

Copyright (c) 2020 Jing Yang and Ming Ye.

## Contributing to SAMMPy

Users are welcome to submit bug reports, feature requests, and code contributions to this project through GitHub or mail to us at mye@fsu.edu.
