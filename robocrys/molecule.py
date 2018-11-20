from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.io.babel import BabelMolAdaptor
from pubchempy import get_compounds


class MoleculeMatcher(object):

    def __init__(self):
        self.matched_molecules = {}

    def get_molecule_name_from_molecule_graph(self,
                                              molecule_graph: MoleculeGraph):
        bma = BabelMolAdaptor.from_molecule_graph(molecule_graph)
        pbmol = bma.pybel_mol
        smiles = pbmol.write(str("smi")).split()[0]

        if smiles in self.matched_molecules:
            return self.matched_molecules[smiles]

        comp = get_compounds(smiles, namespace="smiles")[0]
        self.matched_molecules[smiles] = comp.iupac_name

        return self.matched_molecules[smiles]
