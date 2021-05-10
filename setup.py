""""
robocrystallographer: Automatic generation of crystal structure descriptions.
"""

from os.path import join as path_join

from setuptools import find_packages, setup

with open("README.md", "r") as file:
    long_description = file.read()

version = open("robocrys/_version.py").readlines()[-1].split()[-1].strip("\"'")

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
        "spglib",
        "numpy",
        "scipy",
        "pymatgen>=2020.10.20",
        "inflect",
        "networkx",
        "matminer",
        "monty",
        "pubchempy",
        "pybtex",
        "ruamel.yaml",
    ],
    extras_require={
        "docs": [
            "sphinx==3.5.4",
            "sphinx-argparse==0.2.5",
            "sphinx_rtd_theme==0.5.2",
            "sphinx-autodoc-typehints==1.12.0",
            "m2r2==0.2.7",
        ],
        "dev": ["tqdm", "pybel", "pebble", "maggma"],
        "tests": ["pytest==6.2.3", "pytest-cov==2.11.1"],
        "lint": [
            "coverage==5.5",
            "codacy-coverage==1.3.11",
            "pycodestyle==2.7.0",
            "mypy==0.812",
            "pydocstyle==6.0.0",
            "flake8==3.9.2",
            "pylint==2.8.2",
            "black==21.4b2",
        ],
    },
    package_data={
        "robocrys": [
            path_join("condense", "mineral_db.json.gz"),
            path_join("condense", "molecule_db.json.gz"),
            path_join("condense", "formula_db.json.gz"),
        ]
    },
    data_files=["LICENSE", "CONTRIBUTING.rst"],
    entry_points={"console_scripts": ["robocrys = robocrys.cli:main"]},
)
