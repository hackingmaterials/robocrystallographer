""""
robocrystallographer: Automated structure descriptions.
"""

from setuptools import setup, find_packages


with open('README.rst', 'r') as file:
    long_description = file.read()

setup(
    name='robocrys',
    version='0.0.1',
    description='Automated structure descriptions',
    url='https://github.com/hackingmaterials/robocrystallographer',
    author='Alex Ganose',
    author_email='aganose@lbl.gov',
    long_description=long_description,
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics'
        'Topic :: Scientific/Engineering'
        'Topic :: Other/Nonlisted Topic',
        'Operating System :: OS Independent',
        ],
    keywords='crystal structure crystallography descriptions',
    test_suite='nose.collector',
    packages=find_packages(),
    install_requires=['spglib', 'numpy', 'scipy', 'pymatgen>=2017.12.30',
                      'inflect', 'networkx', 'matminer'],
    extras_require={'docs': ['sphinx', 'sphinx-argparse']},
    package_data={'robocrys': ['mineral_db.json.gz']},
    data_files=['LICENSE', 'requirements_rtd.txt'],
    entry_points={'console_scripts': ['robocrys = robocrys.cli:main']}
    )
