# Copyright (C) 2022-2022 Shujia Huang <huangshujia9@gmail.com>
import os
from argparse import Namespace


DESCRIPTION = "Genotools: A python utility package for NGS data analysis."
meta = Namespace(
    __DISTNAME__     = "genotools",
    __AUTHOR__       = "Shujia Huang",
    __AUTHOR_EMAIL__ = "huangshujia9@gmail.com",
    __URL__          = "https://github.com/ShujiaHuang/genotools",
    __LICENSE__      = "BSD (3-clause)",
    __DOWNLOAD_URL__ = "https://github.com/ShujiaHuang/genotools",
    __VERSION__      = "0.0.0",
)

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages


if __name__ == "__main__":

    THIS_PATH = os.path.abspath(os.path.dirname(__file__))
    long_description = os.path.join(THIS_PATH, "README.md")

    setup(name=meta.__DISTNAME__,
          author=meta.__AUTHOR__,
          author_email=meta.__AUTHOR_EMAIL__,
          maintainer=meta.__AUTHOR__,
          maintainer_email=meta.__AUTHOR_EMAIL__,
          description=DESCRIPTION,
          long_description=(open(long_description).read()),
          long_description_content_type="text/markdown",
          license=meta.__LICENSE__,
          url=meta.__URL__,
          version=meta.__VERSION__,
          download_url=meta.__DOWNLOAD_URL__,
          install_requires=[],
          packages=find_packages(),
          classifiers=[
             "Intended Audience :: Science/Research",
             "Programming Language :: Python :: 3.7",
             "Programming Language :: Python :: 3.8",
             "Programming Language :: Python :: 3.9",
             "License :: OSI Approved :: BSD License",
             "Topic :: Scientific/Engineering :: Bio-Informatics",
             "Topic :: Scientific/Engineering :: Visualization",
             "Topic :: Multimedia :: Graphics",
             "Operating System :: POSIX",
             "Operating System :: Unix",
             "Operating System :: MacOS"],
          )

