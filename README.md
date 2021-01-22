<p align="left"><img width=50% src="images/gallifray_logo.png"></p>

![Python](https://img.shields.io/badge/python-v3.6+-blue.svg)
[![Build Status](https://travis-ci.org/anfederico/Clairvoyant.svg?branch=master)](https://travis-ci.org/Relativist1/Gallifray)
<!--![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)-->
[![GitHub Issues](https://img.shields.io/github/issues/Relativist1/Gallifray?color=yellow)](https://github.com/Relativist1/Gallifray/issues)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
## About

**Gallifray** is a geometric modelling and parameter estimation framework for Black hole images from EHT including VLBI. Gallifray has been developed for the purpose of utilizing the Event Horizon Telescope results to perform analysis. Its uses can be extended to general radio interferometry. The framework utilizes emcee library for bayesian analysis and MCMC sampling, along with eht-imaging library for radio interferometric synthetic data manipulation.

This version of the library is an early release, so please feel free to make an issue or submit a pull request.
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

The the most recent and stable version can be installed manually, by cloning the repository and installing using pip:

```sh
$ git clone https://github.com/Relativist1/Gallifray.git
$ cd Gallifray
$ pip install -r requirements.txt
$ python3 setup.py install
```
## Documentation
Coming soon....

## Models
Included Models

| Model      | Parameters                            |
|------------|---------------------------------------|
| sym_gauss  | I0, A                                 |
| asym_gauss | I0, A, S, phi                         |
| disk       | I0, R                                 |
| crescent   | I0, R_p, psi, tau, phi                |
| xsring     | I0, R_p, R_n, ecn, f, phi             |
| xsringauss | I0, R_p, R_n, ecn, f,gax, aq, gq, phi |


## Basic usage
To run a sample script with MCMC sampling on observational  or simulated dataset 
```sh
$ python example_scripts/example_script_run.py
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
- General Relativistic Raytracing for generating static and spherically symmetric Black hole images surrounded by accretion disk.
- New Models with parameterizations
- User input models
- Inclusion of more sample datasets

## Features to be included
- Reading and writing uvfits file independently
- Support for general relativistic ray-tracing for physically motivated Black hole images.
- Space-VLBI simulations
- Addition of more models

## Contact
If you have any issue or want any help, feel free to drop an email at sbhkmr1999@gmail.com

## License
Gallifray is licensed under GPLv3. See LICENSE.txt for more details.
