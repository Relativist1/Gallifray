import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if __name__ == "__main__":
    setup(name="gallifray",

          version = "0.1",
          author = "Saurabh",
          author_email = "sbhkmr1999@gmail.com",
          description = "Geometric modelling and parameter estimation framework for Black hole images.",
          license = "GPLv3",
          keywords = "astronomy imaging EHT VLBI MCMC",
          url = "https://github.com/Relativist1/Gallifray",
          packages = ["gallifray",
                      "gallifra.const",
                      "gallifray.models",
                      "gallifray.likelihood",
                      "gallifray.priors",
                      "gallifray.mcmc",
                      "gallifray.tardis",
                      "gallifray.utilities"],
          long_description=read('README.md'),
          install_requires=["ehtim",  # https://github.com/achael/eht-imaging.git
                            "emcee",
                            "seaborn",
                            "future",
                            "matplotlib",
                            "numpy",
                            "scipy",
                            "astropy"]
         )
