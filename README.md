# LAMMPS Triangle Lattice with Dipoles

Simulate magnetic origami/kirigami sheets in LAMMPS. Elasticity is modeled using the triangle lattice model 
implemented via `bond_style harmonic` and `dihedral_style harmonic`. Dipolar particles are embedded at triangle centers.

## Installation

Copy the files in `src/` to your LAMMPS `src/` folder and recompile. This adds the `improper_style triangle/center` command, 
which constrains the positions of the dipolar particles to triangle centers. 
Dipole orientation is constrained using the existing `angle_style dipole` command.

## Example Usage

### MATLAB files (`write_lam/`)

The kirigami_sheet_aba2lam.m code generates a LAMMPS data file for a purely elastic kirigami sheet using a triangular mesh in ABAQUS .inp format.

The magnetic_kirigami_sheet_aba2lam.m code generates a LAMMPS data file for a magnetized kirigami sheet using the same `.inp` file.

An `.xyz` file from a purely elastic tensile test simulation is needed to define the deformed state in which the sheet is magnetized (along z).

Rotations are extracted from the deformed configuration and used to set the initial orientation of the dipoles in the flat state.

### Tensile Test (`tensile-test/`)

Contains example LAMMPS input script and data file to run a tensile test simulation. 
