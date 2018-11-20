from typing import List, Dict, Text, Any

from pymatgen import Structure, PeriodicSite


def get_sym_inequiv_components(equivalent_atoms: List[int],
                               components: List[Dict[Text, Any]]
                               ) -> List[Dict[Text, Any]]:
    """

    Component data has to have been generated with inc_site_ids=True

    Args:
        equivalent_atoms:
        components:

    Returns:

    """
    sym_inequiv_components = {}

    for component in components:
        sym_ids = frozenset(
            equivalent_atoms[x] for x in component['site_ids'])

        if sym_ids in sym_inequiv_components:
            sym_inequiv_components[sym_ids]['count'] += 1
            continue

        component['count'] = 1
        component['inequivalent_ids'] = tuple(sym_ids)
        sym_inequiv_components[sym_ids] = component

    return list(sym_inequiv_components.values())


def filter_molecular_components(components):
    molecular_components = [c for c in components if c['dimensionality'] == 0]
    other_components = [c for c in components if c['dimensionality'] != 0]

    return molecular_components, other_components


def get_reconstructed_structure(components: List[Dict[Text, Any]],
                                simplify_molecules: bool=True
                                ) -> Structure:
    if simplify_molecules:
        mol_comps, comps = filter_molecular_components(components)

        if mol_comps:
            lattice = mol_comps[0]['structure'].lattice

            mol_centers = [PeriodicSite(c['structure'][0].specie,
                                        c['molecule_graph'].molecule.center_of_mass,
                                        lattice, coords_are_cartesian=True)
                           for c in mol_comps]
        else:
            mol_centers = []

        sites = [site for c in comps for site in c['structure']] + mol_centers
    else:
        sites = [site for c in components for site in c['structure']]

    return Structure.from_sites(sites)
