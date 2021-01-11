<p align="left"><img width=50% src="gallifray_logo.png"></p>

![Python](https://img.shields.io/badge/python-v3.6+-blue.svg)
[![Build Status](https://travis-ci.org/anfederico/Clairvoyant.svg?branch=master)](https://travis-ci.org/Relativist1/Gallifray)
<!--![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)-->
[![GitHub Issues](https://img.shields.io/github/issues/Relativist1/Gallifray?color=yellow)](https://github.com/Relativist1/Gallifray/issues)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)
[![License](https://img.shields.io/github/license/Relativist1/Gallifray?color=blue)](https://github.com/Relativist1/Gallifray/blob/master/LICENSE)
## About

**Gallifray** is a geometric modelling and parameter estimation framework for Black hole images including VLBI developed for using Event Horizon Telescope results or in general for radio interferometry. The framework utilizes the usage emcee for bayesian analysis and MCMC sampling and eht-imaging library for radio interferometric syntheic data manipulation.

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


## Contributing


## Features to be included
- Reading and writing uvfits file independently
- Support general relativistic ray-tracing for physically motivated Black hole images.
- Addition of more models

## Contact
If you have any issue or want any help, feel free to drop an email at sbhkmr1999@gmail.com

## License
Gallifray is licensed under GPLv3. See LICENSE.txt for more details.
