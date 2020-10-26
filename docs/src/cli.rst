``robocrys`` program
====================

Robocrystallographer can be used on the command-line through the ``robocrys``
command. Most of the options provided in the Python API are also accessible
on the command-line. This page details the basic usage of the program.


.. contents:: Table of Contents
   :local:
   :backlinks: None

Usage
-----

The full range of options supported by ``robocrys`` are detailed in the
`Command-line interface`_ section, and be can be listed using::

    robocrys -h

The package integrates with the `Materials Project
<https://materialsproject.org>`_ to for allow generation of structure
descriptions directly from Materials Project ids. For example, to generate the
description of SnO​​ :sub:`2` (`mp-856
<https://materialsproject.org/materials/mp-856/>`_), one can simply run::

    robocrys mp-856

Alternatively, a structure file can be specified in place of a Materials Project
id. Robocrystallographer supports the same file formats as
`pymatgen <http://pymatgen.org>`_, including the Crystallographic Information
Format (CIF), common electronic structure package formats such as POSCAR
files, and JSON files containing pymatgen Structure objects. Gzipped files are
also supported. For example, to generate the description of a CIF file::

    robocrys MyStructure.cif

Basic options
~~~~~~~~~~~~~

Generating a structure description is split into two steps, the first is to
condense the structure into an intermediate JSON representation, the second is
to generate a description from the intermediate representation. The ``robocrys``
program contains options for controlling both steps of this process.

Condenser options
~~~~~~~~~~~~~~~~~

By default, robocrystallographer determines whether sites are equivalent using
the site information (element, geometry, order parameter), nearest neighbor
information (number and types of neighbors, bond lengths), and next nearest
neighbor information (number and types of neighbors, connectivity, bond angles).
To instead use symmetry to determine which sites are equivalent, the
``--symmetry`` option can be used. The symmetry tolerance is controlled using
the ``--symprec`` option. Note that the symmetry information will always be used
to determine the symmetry labels for the sites, even if the equivalent sites are
determined through structural properties.

Robocrystallographer will produce a description for the structure provided as
is. However, generally best results will be obtained when describing the
conventional cell structure. The ``--conventional`` optional will automatically
convert the structure into the standardized conventional cell representation.

For example, to generate the description for the conventional cell, with
symmetry used to determine site inequivalence, and a symmetry tolerance of
0.001, the command would be::

    robocrys MyStructure.cif --conventional --symmetry --symprec 0.001

Describer options
~~~~~~~~~~~~~~~~~

The ``robocrys`` program has many options to control the verbosity of the
description. By default, oxidation states and symmetry labels
are included in the description. These can be disabled using the ``--no-oxi``
and ``--no-symmetry--labels`` options, respectively.

By default, only the connectivity of cation-polyhedra is included. To also
describe the connectivity anion polyhedra, the ``--anion-polyhedra`` option
can be used.

Robocrystallographer also supports generating LaTeX ready descriptions through
the ``--latexify`` option.

Command-line interface
----------------------

.. argparse::
   :module: robocrys.cli
   :func: _get_parser
   :prog: robocrys
