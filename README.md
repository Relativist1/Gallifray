<p align="left"><img width=50% src="images/gallifray_logo.png"></p>

![Python](https://img.shields.io/badge/python-v3.6+-blue.svg)
[![Build Status](https://travis-ci.org/anfederico/Clairvoyant.svg?branch=master)](https://travis-ci.org/Relativist1/Gallifray)
<!--![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)-->
[![GitHub Issues](https://img.shields.io/github/issues/Relativist1/Gallifray?color=yellow)](https://github.com/Relativist1/Gallifray/issues)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
## About

**Gallifray** is a geometric modelling and parameter estimation framework for Black hole images from EHT including VLBI. Gallifray has been developed for the purpose of utilizing the Event Horizon Telescope results to perform analysis. Its uses can be extended to general radio interferometry. This acts as a bridge between the user-defined models and model fitting via a template provided to do the necessary changes allowing one to simply use their model (phenomenological or geometrical) or the ones pre-defined in the library itself to directly fit with the VLBI datasets (synthetic/observational).
The framework utilizes eht-imaging library for radio interferometric synthetic data manipulation.

This version of the library is an early release, so please feel free to open an issue or submit a pull request.
<br>

## Requirements
### Non-Standard python packages

Gallifray also the following Python packages:

* NumPy
* Astropy
* Matplotlib
* SciPy

### External packages

[eht-imaging](https://github.com/achael/eht-imaging):  python library for simulating and manipulating VLBI data (To be installed manually)

[emcee](https://github.com/dfm/emcee/): the python ensemble sampling toolkit for affine-invariant MCMC

[seaborn](https://github.com/mwaskom/seaborn): for fancy plotting


## Basic installation

The most recent and stable version can be installed manually, by cloning the repository and installing using pip:

```sh
$ git clone https://github.com/Relativist1/Gallifray.git
$ cd Gallifray
$ pip install -r requirements.txt
$ pip install .
```
## Documentation
Coming soon....

## Models
Included Geometric Models

| Model      | Parameters                            |
|------------|---------------------------------------|
| sym_gauss  | I0, A                                 |
| asym_gauss | I0, A, S, phi                         |
| disk       | I0, R                                 |
| crescent   | I0, R_p, psi, tau, phi                |
| xsring     | I0, R_p, R_n, ecn, f, phi             |
| xsringauss | I0, R_p, R_n, ecn, f,gax, aq, gq, phi |


## Basic usage
Example scripts have been provided in the example_scripts/ folder. It contains examples for fitting a geometric model and a raytracing model with ipole.
```sh
$ python example_scripts/example_gauss.py
```
## Example usage of Tardis module 

```python
# samples are either imported or directly used after mcmc sampler
import gallifray as gr

m_true = -0.9594
b_true = 4.294
f_true = 0.534

truths = [m_true, b_true, f_true]
labels = [r"$m_{true}$", r"$b_{true}$",r"$f_{true}$"]

# if Emcee ensemble.sampler is used
samples = sampler.get_chain(flat=True)

gr.Tardis(samples, truths=truths, labels =labels, 
       savefig='new1.png', diag_shade_color='red',
       shade=True, truth1d=True, truth2d = False)
```
<p align="left"><img width=70% src="images/web.png"></p>

## Example plotting scheme using Tardis module 

<p align="left"><img width=80% src="images/tardis_example.png"></p>

## Features coming soon
- New Models with parameterizations
- User input models
- Inclusion of more sample datasets

## Features to be included
- Space-VLBI simulations
- Addition of more models

## Cite
If you use Gallifray in your publication, please cite _`Saurabh. et. al 2022 <https://ui.adsabs.harvard.edu/abs/2022arXiv221206827S/abstract>`_.
       
## Contact
If you have any issues or want any help, feel free to drop an email at sbhkmr1999@gmail.com

## License
Gallifray is licensed under GPLv3. See LICENSE.txt for more details.
