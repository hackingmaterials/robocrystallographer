from typing import List, Dict, Text, Any


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
