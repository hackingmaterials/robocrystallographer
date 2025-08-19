""""
robocrystallographer: Automatic generation of crystal structure descriptions.

"""

from __future__ import annotations

from os.path import join as path_join

from setuptools import find_packages, setup

with open("README.md") as file:
    long_description = file.read()

with open("robocrys/_version.py") as file:
    version = file.readlines()[-1].split()[-1].strip("\"'")

setup(
    name="robocrys",
    version=version,
    description="Automatic generation of crystal structure descriptions",
    url="https://github.com/hackingmaterials/robocrystallographer",
    author="Alex Ganose",
    author_email="aganose@lbl.gov",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="modified BSD",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Information Technology",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering",
        "Topic :: Other/Nonlisted Topic",
        "Operating System :: OS Independent",
    ],
    keywords="crystal-structure crystallography materials-science",
    test_suite="nose.collector",
    packages=find_packages(),
    install_requires=[
        "spglib>=2.5.0",
        "numpy",
        "scipy",
        "pymatgen>=2024.1.26",
        "inflect",
        "networkx",
        "matminer>=0.9.2",
        "monty",
        "pubchempy",
        "pybtex",
    ],
    extras_require={
        "docs": [
            "sphinx>=5.3.0",
            "sphinx-argparse==0.5.2",
            "sphinx_rtd_theme==3.0.2",
            "sphinx-autodoc-typehints==2.3.0",
        ],
        "dev": ["tqdm", "pybel", "pebble", "maggma"],
        "tests": ["pytest==8.3.3", "pytest-cov==6.1.1"],
        "cli": ["mp-api"],
        "lint": [
            "coverage==7.10.4",
            "codacy-coverage==1.3.11",
            "pycodestyle==2.13.0",
            "mypy==1.15.0",
            "pydocstyle==6.1.1",
            "flake8==7.1.1",
            "pylint==3.3.7",
            "black==24.8.0",
        ],
    },
    package_data={
        "robocrys": [
            path_join("condense", "mineral_db.json.gz"),
            path_join("condense", "molecule_db.json.gz"),
            path_join("condense", "formula_db.json.gz"),
        ]
    },
    python_requires=">=3.10",
    data_files=["LICENSE", "CONTRIBUTING.rst"],
    entry_points={"console_scripts": ["robocrys = robocrys.cli:main"]},
)
