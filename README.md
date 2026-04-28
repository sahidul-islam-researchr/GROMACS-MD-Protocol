# GROMACS Protein-Ligand Molecular Dynamics Protocol

[cite_start]This repository contains a standardized workflow for performing and analyzing Molecular Dynamics (MD) simulations for protein-ligand complexes using GROMACS, based on established research protocols[cite: 1].

## 1. System Preparation

### 1.1 Protein Preparation
* [cite_start]Clean and fix the protein structure using `pdbfixer` to handle missing atoms or loops[cite: 3].
* [cite_start]Alternatively, use UCSF Chimera's **Dock Prep** tool to add hydrogens and assign Gasteiger charges[cite: 4, 5].
* [cite_start]Final structures should be saved in `.pdb` format[cite: 6].

### 1.2 Ligand Topology
* [cite_start]Standardize the ligand `.mol2` file by ensuring the `@<TRIPOS>MOLECULE` header is correct and the residue name is set to `LIG`[cite: 9, 10].
* [cite_start]Arrange bond orders using `sort_mol2_bonds.pl` to prevent processing errors[cite: 11].
* [cite_start]Generate force field parameters using **SwissParam** (CHARMM compatible)[cite: 11].

## 2. Simulation Setup

### 2.1 Topology Generation
* Generate protein topology using the CHARMM36 force field and TIP3P water model:
  [cite_start]`gmx pdb2gmx -f REC.pdb -ignh -ff charmm36 -water tip3p`[cite: 14, 15].
* [cite_start]Convert ligand `.pdb` to `.gro` format: `gmx editconf -f LIG.pdb -o LIG.gro`[cite: 17].
* [cite_start]Merge protein and ligand `.gro` files and update `topol.top` to include `LIG.itp` and the correct molecule count[cite: 18, 23, 25].

### 2.2 Solvent and Ions
* Define a triclinic box with 1.0 nm clearance:
  [cite_start]`gmx editconf -f conf.gro -d 1.0 -bt triclinic -o box.gro`[cite: 28, 30].
* [cite_start]Solvate the system: `gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o box_sol.gro`[cite: 33].
* [cite_start]Neutralize with 0.1 M salt concentration: `gmx genion -s ION.tpr -p topol.top -conc 0.1 -neutral -o box_sol_ion.gro`[cite: 36].

## 3. MD Simulation Workflow

1. [cite_start]**Energy Minimization**: Relax the system using `gmx mdrun -v -deffnm EM`[cite: 36].
2. [cite_start]**Restraints**: Generate position restraints for ligand heavy atoms using `gmx genrestr`[cite: 41].
3. **Equilibration**:
   * [cite_start]**NVT**: Thermalize the system[cite: 44].
   * [cite_start]**NPT**: Stabilize the pressure and density[cite: 44].
4. [cite_start]**Production MD**: Run the final simulation trajectory[cite: 45].

## 4. Analysis and Visualization

### 4.1 Trajectory Post-processing
* Correct Periodic Boundary Conditions (PBC) by centering the protein and making the molecules whole:
  [cite_start]`gmx trjconv -s MD.tpr -f MD.xtc -o MD_center.xtc -center -pbc mol -ur compact`[cite: 46, 49].

### 4.2 Stability and Interaction Analysis
* [cite_start]**RMSD**: Calculate for backbone and ligand to assess stability[cite: 47, 50].
* [cite_start]**RMSF**: Analyze per-residue fluctuations[cite: 47, 50].
* [cite_start]**Hydrogen Bonds**: Monitor protein-ligand interactions over time[cite: 47, 51].
* [cite_start]**SASA**: Measure Solvent Accessible Surface Area for both protein and ligand[cite: 51].
* [cite_start]**Radius of Gyration**: Assess the compactness of the complex[cite: 48].

### 4.3 Advanced Analysis
* [cite_start]**MM-GBSA**: Calculate binding free energies using `gmx_MMPBSA`[cite: 53].
* [cite_start]**PCA**: Perform Principal Component Analysis to identify dominant motions[cite: 56].
* [cite_start]**SSE**: Analyze secondary structure elements using DSSP[cite: 58].
* 
