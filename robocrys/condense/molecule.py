"""This module implements a class to match molecule graphs to molecule names.

Some functionality relies on having a working internet connection.
"""

from __future__ import annotations

import warnings

from monty.serialization import loadfn
from importlib.resources import files as import_resource_file
from pubchempy import BadRequestError, get_compounds  # type:ignore[import-untyped]
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.io.babel import BabelMolAdaptor


class MoleculeNamer:
    name_sources = ("traditional", "iupac")

    def __init__(
        self,
        use_online_pubchem: bool = True,
        name_preference: tuple[str, ...] = name_sources,
    ):
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
        db_file = import_resource_file("robocrys.condense") / "molecule_db.json.gz"
        self.molecule_db: dict[str, dict[str, str]] = loadfn(str(db_file))
        self.matched_molecules: dict[str, str] = {}
        self.use_online_pubchem = use_online_pubchem

        # append the sources list to the end in case the user only supplies
        # a single preference
        self.name_preference = tuple(list(name_preference) + list(self.name_sources))

    def get_name_from_molecule_graph(
        self,
        molecule_graph: MoleculeGraph,
    ) -> str | None:
        """Gets the name of a molecule from a molecule graph object.

        Args:
            molecule_graph: A molecule graph.

        Returns:
            The molecule name if a match is found else ``None``.
        """
        smiles = self.molecule_graph_to_smiles(molecule_graph)

        if not smiles:
            return None

        match = None
        if smiles in self.matched_molecules:
            match = self.matched_molecules[smiles]

        elif smiles in self.molecule_db:
            # we should use the first preference for which there is a match
            for source in self.name_preference:
                if (
                    source in self.molecule_db[smiles]
                    and self.molecule_db[source][smiles]
                ):
                    match = self.molecule_db[smiles][source]
                    break

        if not match and self.use_online_pubchem:
            match = self.get_name_from_pubchem(smiles)

        return self._process_match(smiles, match)

    def get_name_from_pubchem(self, smiles: str) -> str | None:
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
        names = {"traditional": traditional, "iupac": comp.iupac_name}

        match = None
        for source in self.name_preference:
            if source in names and names[source]:
                match = names[source]
                break

        if isinstance(match, str):
            match = match.lower()

        return self._process_match(smiles, match)

    @staticmethod
    def molecule_graph_to_smiles(molecule_graph: MoleculeGraph) -> str | None:
        """Converts a molecule graph to SMILES string.

        Args:
            molecule_graph: A molecule graph.

        Returns:
            The SMILES representation of the molecule.
        """
        try:
            bma = BabelMolAdaptor.from_molecule_graph(molecule_graph)
            pbmol = bma.pybel_mol
            return pbmol.write("smi").split()[0]  # type: ignore[attr-defined]
        except RuntimeError:
            warnings.warn(
                "Molecule naming requires openbabel to be installed "
                "with Python bindings. Please get it at http://openbabel.org.",
                RuntimeWarning,
            )
            return None

    def _process_match(self, smiles: str, match: str | None) -> str | None:
        """Utility function to store and process match."""
        if isinstance(match, str):
            match = match.lower()
            self.matched_molecules[smiles] = match
            return self.matched_molecules[smiles]
        return match
