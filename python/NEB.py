#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 20:17:32 2020

@author: dsahu
"""


# compute initial and final states
from ase import Atoms
import numpy as np 
from ase.constraints import FixAtoms
from ase.structure import molecule 
from ase.calculators.vasp import Vasp
 
 
atoms = molecule('NH3')
constraint = FixAtoms(mask=[atom.symbol == 'N' for atom in atoms])
atoms.set_constraint(constraint)

Npos = atoms.positions[0]

# move N to origin
atoms.translate(-Npos)
atoms.set_cell((10, 10, 10), scale_atoms=False)

atoms2 = atoms.copy()
pos2 = atoms2.positions

for i,atom in enumerate(atoms2): 
    if atom.symbol == 'H': 
# reflect through z 
        pos2[i] *= np.array([1, 1, -1])
        atoms2.positions = pos2
        
#now move N to center of box 
atoms.translate([5, 5, 5]) 
atoms2.translate([5, 5, 5])

calcs = [Vasp('molecules/nh3-initial',
              xc='PBE', 
              encut=350, 
              ibrion=1, 
              nsw=10,
              atoms=atoms),
Vasp('molecules/nh3-final',
     xc='PBE',
     encut=350,
     ibrion=1,
     nsw=10,
     atoms=atoms2)]

 print [c.potential_energy for c in calcs]
 
 
 
 
 
 from ase.build import fcc100, add_adsorbate
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton

# 2x2-Al(001) surface with 3 layers and an
# Au atom adsorbed in a hollow site:
slab = fcc100('Al', size=(2, 2, 3))
add_adsorbate(slab, 'Au', 1.7, 'hollow')
slab.center(axis=2, vacuum=4.0)

# Make sure the structure is correct:
#view(slab)

# Fix second and third layers:
mask = [atom.tag > 1 for atom in slab]
#print(mask)
slab.set_constraint(FixAtoms(mask=mask))

# Use EMT potential:
slab.set_calculator(EMT())

# Initial state:
qn = QuasiNewton(slab, trajectory='initial.traj')
qn.run(fmax=0.05)

# Final state:
slab[-1].x += slab.get_cell()[0, 0] / 2
qn = QuasiNewton(slab, trajectory='final.traj')
qn.run(fmax=0.05)


#Now, do the NEB calculation:

from ase.io import read
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.optimize import BFGS

initial = read('initial.traj')
final = read('final.traj')

constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])

images = [initial]
for i in range(3):
    image = initial.copy()
    image.set_calculator(EMT())
    image.set_constraint(constraint)
    images.append(image)

images.append(final)

neb = NEB(images)
neb.interpolate()
qn = BFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)

#Visualize the results with:

ase gui 
neb.traj@-5:
    
import matplotlib.pyplot as plt
from ase.neb import NEBTools
from ase.io import read

images = read('neb.traj@-5:')

