 from math import sqrt
 from ase import Atoms
 import numpy as np
 from ase.build import fcc111,add_adsorbate
 from ase.build import cut
 from ase.spacegroup import crystal
 from ase.visualize import view
 from ase.io import write
 from ase.lattice.surface import fcc111, add_adsorbate, molecule
 from ase.build import fcc111, root_surface
 from ase.build import bulk
 
## use this command from the main terminal to expand the slab size
 # ase gui -r 3,3,2 slab.xyz 
 
 #### Cu slab cases
 
# Oxygen on top of Cu
cu = bulk ('Cu', a =3.6)
view(cu)

slab = fcc111('Cu', size=(3,3,5)) # 5­layers of 3x3 Cu (111) surface
# add O atom 2.5 Å above the surface in the 'bridge' site
add_adsorbate(slab, 'O', 2.5, position='bridge') # bridge position
view(slab)
--------------------------------------------------------------

# H on top of Cu
slab = fcc111('Cu', size=(2,2,3), vacuum=10.0) # 3layers of 2x2 Cu (111) surface
# add O atom 1.5 Å above the surface in the 'bridge' site
add_adsorbate(slab, 'H', 1.5, 'ontop') # position is mentioned like above
slab.center(vacuum=10.0, axis=2)
view(slab)

#position='fcc', position='bridge', position='hcp'
---------------------------------------------------------------------
# orthogonal unit cell
slab = fcc111('Cu', size=[2, 4, 3], a=3.55, orthogonal=True)
view(slab)
-----------------------------------------------------------------------
# N2 adsorption on FCC Cu (111) surface
slab = fcc111 ('Cu', size =(4 , 4 , 4) , a=4.0 , vacuum =6.0)
add_adsorbate (slab , molecule ('N2') , height =3.0 , offset =(2 , 2) ,
               position = 'ontop') # offset value move the adsorbate here from edge to middle
view(slab)
print(slab) # print the cell size
-----------------------------------------------------------------------
# Cu-442 (N2 on top of Slab)
h = 1.85
d = 1.10   # N-N bond length

slab = fcc111('Cu', size=(4, 4, 2), vacuum=10.0)
view(slab)
molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)]) #define a N2 molecule by directly specifying the position of two nitrogen atoms
add_adsorbate(slab, molecule, h, 'ontop')
view(slab)
----------------------------------------------------------------------
# Making Manual Cu-Slab [without taking any library help]
a = 3.5
atoms = Atoms('Cu4',
cell=[sqrt(2) * a, sqrt(2) * a, 1.0, 90, 90, 120],
pbc=(1, 1, 0),
scaled_positions=[(0, 0, 0),
(0.5, 0, 0),
(0, 0.5, 0),
(0.5, 0.5, 0)])
atoms.center(vacuum=10.0, axis=2)
atoms.cell
atoms.positions
atoms[0]
atoms.write('slab.xyz')
view(atoms)
----------------------------------------------------------------------
# Manually create an copper (111) slab with three layers
#
# First an unit cell of Cu
a = 4.05
copper = crystal('Cu', [(0,0,0)], spacegroup=225,
                    cellpar=[a, a, a, 90, 90, 90])
view(copper)
# Then cut out the slab
cu111 = cut(copper, (1,-1,0), (0,1,-1), nlayers=3)
view(cu111)

# Add add a Ag aton on 1.9 angstrom above of the Cu slab
 h = 1.9
 relative = (1 / 6, 1 / 6, 0.5)
 absolute = np.dot(relative, atoms.cell) + (0, 0, h)
 atoms.append('Ag')
 atoms.positions[-1] = absolute
 view(atoms)
----------------------------------------------------------------------
# Repeat Cu 111 Primitive cell along X & Y direction
a = 4.6
c = 2.95
atoms = fcc111('Cu', (1, 1, 3),vacuum=30.0)
atoms = root_surface(atoms, 27) # important:
view(atoms)

## cubic 8­atom Cu unit cell
atoms = bulk('Cu', cubic=True)  # As per element input (here Cu), atom numbers dependent
# Generating the 256 atom supercell of Cu
supercell = atoms.repeat((4, 4, 4)) 
view(supercell)
print(supercell)
del supercell[0]  # remove first atom, e.g. create a vacancy
view(supercell)
----------------------------------------------------------------------



## Other than Cu slab 
# Only N2 molecule
atoms = Atoms('N2', positions=[[0, 0, -1], [0, 0, 1]])
atoms.center(vacuum=3.0)
view(atoms)

## Making Rutile TiO2 slab
a = 4.6
c = 2.95

atoms = crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
view(atoms)
----------------------------------------------------------------------
## Create a skutterudite unit cell of Co & Sb
a = 9.04
skutterudite = crystal(
    ('Co', 'Sb'),
    basis=[(0.25,0.25,0.25), (0.0, 0.335, 0.158)],
    spacegroup=204,
    cellpar=[a, a, a, 90, 90, 90])

# Then use *origo* to put 'Co' at the corners and *extend* to
# include all corner and edge atoms.
s = cut(skutterudite, origo=(0.25, 0.25, 0.25), extend=1.01)
view(s)  
----------------------------------------------------------------------
# FCC gold:

a = 4.05  # Gold lattice constant
b = a / 2
fcc = Atoms('Au',
            cell=[(0, b, b), (b, 0, b), (b, b, 0)],
            pbc=True)
view(fcc)
----------------------------------------------------------------------
# Hydrogen wire:

d = 0.9  # H-H distance
h = Atoms('H', positions=[(0, 0, 0)],
          cell=(d, 0, 0),
          pbc=(1, 0, 0))
view(h)

# Graphite or Diamond??
from ase.lattice.hexagonal import Graphite
atoms = Graphite('C', latticeconstant={'a': 2.4612,
'c': 6.7079})
view(atoms)
----------------------------------------------------------------------