 from math import sqrt
 from ase import Atoms, Atom 
 import numpy as np
 from ase.build import fcc111,add_adsorbate
 from ase.build import cut
 from ase.spacegroup import crystal
 from ase.visualize import view
 from ase.io import read, write
 from ase.lattice.surface import fcc111, add_adsorbate
 from ase.build import bulk
 from ase.structure import molecule 
 from numpy import pi
 from ase.calculators.vasp import Vasp
 from ase.lattice.cubic import FaceCenteredCubic


# Making CO and save it as image file in local computer
#C O molecule with the C at the origin of unit cell
# define an Atoms object: type and position of each atom
atoms = Atoms([Atom('C', [0., 0., 0.]),
               Atom('O', [1.1, 0., 0.])],
              cell=(10, 10, 10))
print('V = {0:1.0f} Angstrom^3'.format(atoms.get_volume())) 
write('/Users/dsahu/Desktop/simple-cubic-cell.png', atoms, show_unit_cell=2)

----------------------------------------------------------------------

#centering the atoms of CO in the unit cell
# we have guessed values for until the CO molecules are on average 10 Å apart.
b = 7.1 
atoms = Atoms([Atom('C', [0., 0., 0.]), 
               Atom('O', [1.1, 0., 0.])],
              cell=[[b, b, 0.], # face-centered cubic lattice
                    [b, 0., b], 
                    [0., b, b]])
print('V = {0:1.0f} Ang^3'.format(atoms.get_volume())) # Note the final volume is only about 715 Å , which is smaller than the cube
atoms.center() # translate atoms to center of unit cell 
view(atoms)
# Write POSCAR in cartesian coordinate
write('/Users/dsahu/Desktop/POSCAR',images=atoms,format='vasp',vasp5=True)
# Write POSCAR in fractional coordinate
write('/Users/dsahu/Desktop/POSCAR1',images=atoms,format='vasp',direct=True,vasp5=True)
# save png image file
write('/Users/dsahu/Desktop/fcc-cell.png', atoms, show_unit_cell=2)
os.system('open -a VESTA.app POSCAR')
#compute the distance of a vector as the square root of the sum of squared elements
# get unit cell vectors and their lengths
(a1, a2, a3) = atoms.get_cell() 
print('|a1| = {0:1.2f} Ang'.format(np.sum(a1**2)**0.5)) 
print('|a2| = {0:1.2f} Ang'.format(np.linalg.norm(a2))) 
print('|a3| = {0:1.2f} Ang'.format(np.sum(a3**2)**0.5))
----------------------------------------------------------------------


# Import a coordintae of any molecule (e.g. xyz format) from a local computer and place it in unit cell 
# 
atoms = read('/Users/dsahu/Desktop/molecule.xyz')
view(atoms)
atoms.center(vacuum=5)
write('/Users/dsahu/Desktop/molecule.png', atoms, show_unit_cell=2) # Save the image to see the molecule in cell

----------------------------------------------------------------------
#an acetonitrile molecule from libraty
atoms = molecule('CH3CN')
view (atoms)
atoms.center(vacuum=6) 
print('unit cell') 
print('---------') 
print(atoms.get_cell())
write('/Users/dsahu/Desktop/ch3cn.png', atoms, show_unit_cell=2)
# Rotate the unit cell 45 degree [not atoms] & save image file
write('/Users/dsahu/Desktop/ch3cn-rotated.png', atoms, show_unit_cell=2, rotation='45x,45y,0z')
----------------------------------------------------------------------
#rotates the molecule an angle (in radians) around a vector
atoms = molecule('CH3CN') 
atoms.center(vacuum=6) 
p1 = atoms.get_positions()
# Rotate atoms [not the unit cell]
atoms.rotate('x', pi/4, center='COM', rotate_cell=False)
atoms.rotate('y', pi/4, center='COM', rotate_cell=False)
write('/Users/dsahu/Desktop/ch3cn-rotated-2.png', atoms, show_unit_cell=2)
print('difference in positions after rotating')
print('atomdifference vector')
print('--------------------------------------')
p2 = atoms.get_positions()
diff = p2 - p1 
for i, d in enumerate(diff): print('{0} {1}'.format(i, d))
----------------------------------------------------------------------

atoms1 = molecule('NH3')
atoms2 = molecule('O2') atoms2.translate([3, 0, 0])
bothatoms = atoms1 + atoms2 
bothatoms.center(5)
view(bothatoms)
write('/Users/dsahu/Desktop/bothatoms.png', bothatoms, show_unit_cell=2, rotation='90x')


atoms1 = molecule('NH3')
atoms2 = molecule('O2') 
atoms2.translate([3, 0, 0])
bothatoms = atoms1 + atoms2 
bothatoms.center(5) write('/Users/dsahu/Desktop/bothatoms.png', bothatoms, show_unit_cell=2, rotation='90x')
----------------------------------------------------------------------

from numpy import linalg as LA
a = np.arange(9) - 4
a
LA.norm(a)
b = a.reshape((3, 3))
b

# Read a structure file from the local computer e.g. file.xyz
# Save the structure as png format in the local computer
atoms = read('/Users/dsahu/Desktop/file.xyz')
atoms.center(vacuum=5)
write('/Users/dsahu/Desktop/file-xyz.png', atoms, show_unit_cell=2)
----------------------------------------------------------------------
from ase.data import g2 
keys = g2.data.keys() 
# print in 3 columns
for i in range(len(keys) / 3): 
    print('{0:25s}{1:25s}{2:25s}'.format(*tuple(keys[i * 3: i * 3 + 3])))


atoms = FaceCenteredCubic(directions=[[0, 1, 1], 
                                      [1, 0, 1], 
                                      [1, 1, 0]], 
                          size=(1, 1, 1), 
                          symbol='Cu')
view(atoms)


lp= 4.0
refcell= np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
init_cell= lp*refcell
pos=[(0.,0.,0.),(0.5,0.5,0.),(0.5,0.,0.5),(0.,0.5,0.5)]
elems= ['Al','Al','Al','Al']
atoms= crystal(elems,pos,cell=init_cell)
view(atoms)

atoms=FaceCenteredCubic(symbol='Al',latticeconstant=4.041,size=(1,1,1))
atoms = Atoms(symbols=elems,cell=init_cell,scaled_positions=pos)
atoms = Atoms(symbols=elems,cell=init_cell,positions=pos)



atoms = FaceCenteredCubic(symbol="Al",
                          directions=[[2,3,0],[-3,1,0],[0,0,1]],
                          latticeconstant=4.04,
                          size=(1,1,2))
view(atoms)
----------------------------------------------------------------------
def write_vasp(poscar, atoms, label='', direct=False, sort=None,
               symbol_count=None, long_format=True, vasp5=False)    

from ase.lattice.hexagonal import HexagonalClosedPacked
atoms = HexagonalClosedPacked(symbol='Mg',
                              latticeconstant={'a':3.200,'c':5.188},
                              size=(2,2,2))
atoms.append(Atom('Al', (1/6., 1/6., .1))) # Add the Al atom at origin
atoms.pop(0) # Delete the Mg atom at origin
print atoms.cell
view(atoms)



atoms[3].symbol= 'Mg'
# or
atoms[-1].symbol= 'Mg'

atoms[2].symbol= 'Mg'  # <---Replace the atom symbol
write('/Users/dsahu/Desktop/POSCAR',images=atoms,format='vasp',direct=True,vasp5=True)
----------------------------------------------------------------------
from ase.lattice.cubic import FaceCenteredCubicFactory
class CaF2Factory(FaceCenteredCubicFactory):
    bravais_basis=[[0,0,0],
                   [0.25,0.25,0.25],
                   [0.75,0.75,0.75]]
    element_basis=(0,1,1)

C1 = CaF2 = SiMg2 = CaF2Factory()
atoms= C1(symbol=['Si','Mg'],latticeconstant=6.365,size=(1,1,1))
atoms.cell
atoms.get_scaled_positions()
atoms.get_chemical_symbols()
view (atoms)


# Making Cu-slab with repeatative CO adsorbed ***
adsorbate = Atoms('CO')
adsorbate[1].z = 1.1

adsorbate = atoms
a = 3.61
slab = fcc111('Cu', (2, 2, 3), a=a, vacuum=7.0)
add_adsorbate(slab, adsorbate, 1.8, 'ontop')
s = slab * (3, 3, 1)
view(s)

write('/Users/dsahu/Desktop/movie.gif', [bulk(s) for s in ['Cu', 'Ag', 'Au']], interval=500)

view(slab * (3, 3, 1), rotation='10z,-80x')

------------------------------------------------------------------------------
# Making Pt (110) slab
from ase.lattice.surface import fcc110 
slab = fcc110('Pt', (2, 1, 7), a=4.0, vacuum=6.0)
view(slab)


from ase.lattice.surface import add_adsorbate 
add_adsorbate(slab, 'H','hollow')


------------------------------------------------------------------------------
# To setup a Au(211) surface with 9 layers and 10 Å of vacuum:

from ase.build import surface
s1 = surface('Au', (2, 1, 1), 9)
s1.center(vacuum=10, axis=2)
#This is the easy way, where you use the experimental lattice constant for gold bulk structure. You can write:

from ase.visualize import view
view(s1)
s1.edit() # or simply : if you want to see and rotate the structure.
------------------------------------------------------------------------------
#Next example is a molybdenum bcc(321) surface 
#where we decide what lattice constant to use: a=3.16

from ase.build import bulk
Mobulk = bulk('Mo', 'bcc', a=3.16, cubic=True)
s2 = surface(Mobulk, (3, 2, 1), 9)
s2.center(vacuum=10, axis=2)
------------------------------------------------------------------------------
# creation of alloy surfaces is also very easily carried out with this module. 
#In this example, two Pt3Rh fcc(211) surfaces will be created:

a = 4.0
from ase import Atoms
Pt3Rh = Atoms('Pt3Rh',
              scaled_positions=[(0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5)],
              cell=[a, a, a],
              pbc=True)
s3 = surface(Pt3Rh, (2, 1, 1), 9)
s3.center(vacuum=10, axis=2)

Pt3Rh.set_chemical_symbols('PtRhPt2')
s4 = surface(Pt3Rh, (2, 1, 1), 9)
s4.center(vacuum=10, axis=2)
view(s4)
------------------------------------------------------------------------------
## Wulf construction [not relevent now...]
from ase.cluster.wulff import wulff_construction
from ase.io import write

atoms = wulff_construction('Pd',
                           surfaces=[(1, 0, 0),
                                     (1, 1, 1),
                                     (1, 1, 0)],
                           energies=[0.1, 0.5, 0.15],
                           size=100,
                           structure='fcc',
                           rounding='below')

view(atoms)


write('images/wulff.png', atoms)
