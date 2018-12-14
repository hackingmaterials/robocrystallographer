"""
This module implments a class to match molecule graphs to molecule names.

Some functionality relies on having a working internet connection.
"""

from typing import Optional, Tuple

from pkg_resources import resource_filename
from pubchempy import get_compounds, BadRequestError

from pymatgen import loadfn
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.io.babel import BabelMolAdaptor


class MoleculeNamer(object):
    name_sources = ('traditional', 'iupac')

    def __init__(self, use_online_pubchem: bool = True,
                 name_preference: Tuple[str] = name_sources):
        """Class to match molecule graphs to known molecule names.

        Args:
            use_online_pubchem: Whether to try using the Pubchem website for
                matching molecules if a match is not found in the offline
                database. Defaults to ``True``. Requires a working internet
                connection and for the ``pubchempy`` package to be installed.
            name_preference: The order of preference for determining compound
                names. Options are "traditional", and "iupac". If the
                first option is not available, the subsequent options will be
                used. Should be provided as a tuple of options, from 1st choice
                to last.
        """

        db_file = resource_filename('robocrys.condense', 'molecule_db.json.gz')
        self.molecule_db = loadfn(db_file)
        self.matched_molecules = {}
        self.use_online_pubchem = use_online_pubchem

        # append the sources list to the end in case the user only supplies
        # a single preference
        self.name_preference = tuple(list(name_preference) +
                                     list(self.name_sources))

    def get_name_from_molecule_graph(self, molecule_graph: MoleculeGraph,
                                     ) -> Optional[str]:
        """Gets the name of a molecule from a molecule graph object.

        Args:
            molecule_graph: A molecule graph.

        Returns:
            The molecule name if a match is found else ``None``.
        """
        smiles = self.molecule_graph_to_smiles(molecule_graph)

        match = None
        if smiles in self.matched_molecules:
            match = self.matched_molecules[smiles]

        elif smiles in self.molecule_db:

            # we should use the first preference for which there is a match
            for source in self.name_preference:
                if (source in self.molecule_db[smiles] and
                        self.molecule_db[source][smiles]):
                    match = self.molecule_db[smiles][source]
                    break

        if not match and self.use_online_pubchem:
            match = self.get_name_from_pubchem(smiles)

        return self._process_match(smiles, match)

    def get_name_from_pubchem(self, smiles: str) -> Optional[str]:
        """Tries to get the name of a molecule from the Pubchem website.

        Args:
            smiles: A SMILES string.

        Returns:
            The molecule name if a match is found else ``None``.
        """
        try:
            comp = get_compounds(smiles, namespace="smiles")[0]

        except (BadRequestError, IndexError):
            return None

        traditional = comp.synonyms[0] if comp.synonyms else None
        names = {'traditional': traditional,
                 'iupac': comp.iupac_name}

        match = None
        for source in self.name_preference:
            if source in names and names[source]:
                match = names[source]
                break

        if isinstance(match, str):
            match = match.lower()

        return self._process_match(smiles, match)

    @staticmethod
    def molecule_graph_to_smiles(molecule_graph: MoleculeGraph
                                 ) -> Optional[str]:
        """Converts a molecule graph to SMILES string.

        Args:
            molecule_graph: A molecule graph.

        Returns:
            The SMILES representation of the molecule.
        """
        bma = BabelMolAdaptor.from_molecule_graph(molecule_graph)
        pbmol = bma.pybel_mol
        return pbmol.write(str("smi")).split()[0]

    def _process_match(self, smiles: str, match: Optional[str]
                       ) -> Optional[str]:
        """Utility function to store and process match."""
        if isinstance(match, str):
            match = match.lower()
            self.matched_molecules[smiles] = match
            return self.matched_molecules[smiles]
        else:
            return match
