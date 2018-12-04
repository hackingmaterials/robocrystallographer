from monty.serialization import loadfn

from pkg_resources import resource_filename

from robocrys.mineral import MineralMatcher
from robocrys.site import SiteAnalyzer
from robocrys.structure import StructureCondenser


formula_mapping = loadfn(resource_filename('robocrys', 'formula_db.json.gz'))
