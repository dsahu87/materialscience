 from ase.lattice.surface import fcc111, add_adsorbate
 #from ase.lattice.surface import molecule
 from math import sqrt
 from ase import Atoms, Atom 
 import numpy as pi
 import numpy as np
 from ase.build import fcc111,add_adsorbate,cut,bulk,root_surface
 from ase.spacegroup import crystal
 from ase.visualize import view
 from ase.io import read, write
 from ase.structure import molecule 
 from ase.calculators.vasp import Vasp
 from ase.lattice.cubic import FaceCenteredCubic

# conda update -n base -c defaults conda #update conda to use the command
 # conda install -c conda-forge ase # Install ASE Package to use the command

# Import a coordintae of any molecule (e.g. xyz/mol/pdb format) from a local computer and place it in unit cell 
# 
# atoms = read('/Users/dsahu/Desktop/Dibromo-Bianthryl-Precursor-Monomer.mol') 
adsorbate = read('/Users/dsahu/Desktop/Br-Br.mol') 
#adsorbate2 = read('/Users/dsahu/Desktop/ASE_Python/CONTCAR.xyz') 
adsorbate.rotate(35, 'z', center='COM')


# a = 3.61
# slab = fcc111('Cu', (12, 12, 3), a=a, vacuum=30.0)
slab = fcc111('Cu', (9, 9, 3), vacuum=10.0)
# s = slab * (3, 3, 1)
add_adsorbate(slab, adsorbate, 5.5, 'ontop', offset =(4.5 , 4))
#add_adsorbate(slab, adsorbate2, 6.5, 'ontop', offset =(7 , 6.5))
#slab.wrap()
view(slab) 
print(slab) # See the cell information
# Write POSCAR in cartesian coordinate
write('/Users/dsahu/Desktop/POSCAR_prf',images=slab,format='vasp',vasp5=True, sort = True)
# Write POSCAR in fractional coordinate
write('/Users/dsahu/Desktop/POSCAR_Br',images=slab,format='vasp',direct=True,vasp5=True,sort=True)

### Increase slab size
a = 3.61
slab = fcc111('Cu', (6, 6, 3), a=a, vacuum=30.0) # It is same slab like above
view(slab)

# Insert 2 molecue of adsorbate on the Cu slab
# Following possibilites:
# Copy & paste from the single import molecule. possible?
# Import 2 molecule together
# Import one by one molecue


# Python programs which read your VASP input/output and then work on them programmatically. In fact, with ASE, it is almost trivial to make supercells. You can do it with 3 lines of Python code:

import ase.io.vasp
cell = ase.io.vasp.read_vasp("POSCAR")
ase.io.vasp.write_vasp("POSCAR.4x4x4",cell*(4,4,4), label='444supercell',direct=True,sort=True)
#The code above reads a POSCAR file in the current working directory, 
# transforms it to a 4x4x4 supercell and writes the results to disk as “POSCAR.4x4x4”.

ase.build.rotate(adsorbate2, a1, a2, b1, b2, rotate_cell=True, center=0, 0, 0)#[source]
# Rotate atoms, such that a1 will be rotated in the direction of a2 and b1 in the direction of b2. 
# The point at center is fixed. Use center=’COM’ to fix the center of mass. If rotate_cell is true,
# the cell will be rotated together with the atoms.
# Note that the 000-corner of the cell is by definition fixed at origo.
# Hence, setting center to something other than (0, 0, 0) will rotate the atoms out of the cell, 
# even if rotate_cell is True.