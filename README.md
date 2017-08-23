Pasta is a Package for Automatic Structure Transformation and Analysis.

## Installation
```bash
cd pasta
python setup.py install
```

## Usage

```python
from pasta.structure import Structure
# import structure from vasp type file
s = Structure.import_from_vasp('POSCAR')
# export structure to pwmat type file
s.export_to_pwmat('atom.config')
# create supercell
s.create_supercell([[1,1,0],[0,2,0],[0,0,1]])
# get crystal system
s.get_crystal_system()
# find radial distribution function
s.find_rdf()
```

RAmen!

