from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.analysis.graphs import MoleculeGraph
from pubchempy import get_compounds


class MoleculeMatcher(object):

    def __init__(self):
        self.matched_molecules = {}

    def get_name_from_molecule_graph(self, molecule_graph: MoleculeGraph):
        smiles = self.molecule_graph_to_smiles(molecule_graph)

        if smiles in self.matched_molecules:
            return self.matched_molecules[smiles]

        self.matched_molecules[smiles] = self.smiles_to_molecule_name(smiles)
        return self.matched_molecules[smiles]

    @staticmethod
    def molecule_graph_to_smiles(molecule_graph: MoleculeGraph):
        bma = BabelMolAdaptor.from_molecule_graph(molecule_graph)
        pbmol = bma.pybel_mol
        return pbmol.write(str("smi")).split()[0]

    @staticmethod
    def smiles_to_molecule_name(smiles):
        """"""
        comp = get_compounds(smiles, namespace="smiles")[0]
        return comp.iupac_name
