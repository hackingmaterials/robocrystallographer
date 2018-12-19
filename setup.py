""""
robocrystallographer: Automatic generation of crystal structure descriptions.
"""

from setuptools import setup, find_packages
from os.path import join as path_join


with open('README.rst', 'r') as file:
    long_description = file.read()

setup(
    name='robocrys',
    version='0.1.0',
    description='Automatic generation of crystal structure descriptions',
    url='https://github.com/hackingmaterials/robocrystallographer',
    author='Alex Ganose',
    author_email='aganose@lbl.gov',
    long_description=long_description,
    license='modified BSD',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering',
        'Topic :: Other/Nonlisted Topic',
        'Operating System :: OS Independent',
        ],
    keywords='crystal-structure crystallography materials-science',
    test_suite='nose.collector',
    packages=find_packages(),
    install_requires=['spglib', 'numpy', 'scipy', 'pymatgen>=2017.12.30',
                      'inflect', 'networkx', 'matminer', 'monty', 'pubchempy'],
    extras_require={'docs': ['sphinx', 'sphinx-argparse',
                             'sphinx-autodoc-types>=2.1.1'],
                    'dev': ['tqdm', 'pybel', 'pebble', 'maggma']},
    package_data={'robocrys': [path_join('condense', 'mineral_db.json.gz'),
                               path_join('condense', 'molecule_db.json.gz'),
                               path_join('condense', 'formula_db.json.gz')]},
    data_files=['LICENSE', 'CONTRIBUTING.rst'],
    entry_points={'console_scripts': ['robocrys = robocrys.cli:main']}
    )
