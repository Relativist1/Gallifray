import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if __name__ == "__main__":
    setup(name="gallifray",

          version = "1.0",
          author = "Saurabh",
          author_email = "sbhkmr1999@gmail.com",
          description = "Geometric modelling and parameter estimation framework for Black hole images.",
          license = "GPLv3",
          keywords = "astronomy imaging EHT VLBI MCMC",
          url = "https://github.com/relativist1/gallifray",
          download_url = "https://github.com/achael/eht-imaging/archive/v1.2.1.tar.gz",
          packages = ["gallifray",
                      "gallifray.models",
                      "utils" ],
          long_description=read('README.rst'),
          install_requires=["ehtim",  # https://github.com/achael/eht-imaging.git
                            "emcee",
                            "future",
                            "matplotlib",
                            "numpy",
                            "scipy",
                            "astropy",
                            "corner",
                            "ehtplot"] # git+https://github.com/liamedeiros/ehtplot.git
         )

