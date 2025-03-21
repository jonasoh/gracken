from setuptools import setup, find_packages
import os

# Read the long description from README.md if present
this_directory = os.path.abspath(os.path.dirname(__file__))
try:
    with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "See project page for more details."

setup(
    name="gracken",
    version="0.1.1",
    packages=find_packages(),
    install_requires=[
        "ete3==3.1.3",
        "pandas==2.2.3",
    ],
    description="Create phylogenetic trees from metagenomic reports.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={"console_scripts": ["gracken=gracken.gracken:main"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Environment :: Console",
    ],
    python_requires=">=3.9, <3.13",  # ete3 requires <3.13
)
