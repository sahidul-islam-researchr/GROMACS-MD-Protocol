# GROMACS MD Simulation & Analysis Protocol for Protein-Ligand System

This repository contains all commands and workflows to run a molecular dynamics (MD) simulation and subsequent analysis using **GROMACS (CHARMM36 + TIP3P)**.  
The protocol is designed for publication-ready results.

# GROMACS MD Simulation & Analysis – Complete Command List

This file contains all the commands needed to run a protein-ligand MD simulation and perform subsequent analysis using GROMACS (CHARMM36, TIP3P).

## 🚀 Simulation Commands

Run these in order after preparing `fixed_REC.pdb`, `LIG.pdb`, and obtaining `LIG.itp` from SwissParam.

```bash
# Generate protein topology
gmx pdb2gmx -f fixed_REC.pdb -ignh -ff charmm36 -water tip3p

# Convert ligand pdb to gro
gmx editconf -f LIG.pdb -o LIG.gro

# --- MANUAL STEP ---
# Append LIG.gro content (excluding first two lines) into conf.gro.
# Update total atom count in conf.gro.

# Edit topol.top: add '#include "LIG.itp"' and under [ molecules ] add 'LIG    1'

# Edit LIG.itp: change 'lig_gmx2 3' to 'LIG 3' inside [ moleculetype ]

# Define box
gmx editconf -f conf.gro -d 1.0 -bt triclinic -o box.gro

# Solvate
gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o box_sol.gro

# Add ions
gmx grompp -f ions.mdp -c box_sol.gro -p topol.top -o ION.tpr
echo "15" | gmx genion -s ION.tpr -p topol.top -conc 0.1 -neutral -o box_sol_ion.gro

# Energy minimization
gmx grompp -f EM.mdp -c box_sol_ion.gro -p topol.top -o EM.tpr
gmx mdrun -v -deffnm EM

# NVT equilibration
gmx grompp -f NVT.mdp -c EM.gro -r EM.gro -p topol.top -o NVT.tpr
gmx mdrun -deffnm NVT

# NPT equilibration
gmx grompp -f NPT.mdp -c NVT.gro -r NVT.gro -p topol.top -maxwarn 2 -o NPT.tpr
gmx mdrun -deffnm NPT

# Production MD
gmx grompp -f MD.mdp -c NPT.gro -t NPT.cpt -p topol.top -o MD.tpr
gmx mdrun -v -deffnm MD
