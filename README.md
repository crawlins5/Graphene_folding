# Graphene_folding
Resource for generating folded sheets of graphene, running a simulation with LAMMPS and extracting information related to friction.
This respository includes the folders:

example/
INPUT: A sample LAMMPS input file

param_g: A sample parameter file for generating folded graphene sheets

ffield: ReaxFF parameter file for C/H/O (references therein)

Graphene_generator/

folded_general.py: Set-up script for reading in param_g file

folded_arm_zig.py: Performs the calculations for generating folded sheets with an offset angle of 0 or 30 degrees and prints co-ordinate file (graphene.dat).

folded_otherangles.py: Similar to folded_arm_zig but for arbitrary angles. Lacks certain speed-ups available for simpler angles.

Friction/
friction_fold_graphene.py: Reads in INPUT and graphene_xyz_v.lammpstrj files and outputs Distance.dat and velocity.dat files.
